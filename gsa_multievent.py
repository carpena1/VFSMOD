#!/usr/bin/env python3
"""
gsa_multievent.py — Global Sensitivity Analysis for VFSMOD Freundlich vs Linear
Runs N_SAMPLES Latin Hypercube parameter sets, each as a 3-event sequential
simulation, for both linear (IKD=1) and Freundlich (IKD=2) isotherms.

Usage:
  python3 gsa_multievent.py [--samples 128] [--events 3] [--output gsa_results.csv]

Outputs:
  gsa_results.csv  — one row per sample, columns for all inputs and outputs
  gsa_summary.txt  — Sobol first-order and total sensitivity indices
"""

import sys, os, re, shutil, subprocess, argparse, csv, math, time
import numpy as np

BASE_PRJ  = 'sampleP'
VFSM_EXE  = './vfsm'
C_REF     = 1.0   # mg/L reference concentration for Kf=Kd equivalence

# ── Parameter space ─────────────────────────────────────────────────────────
# Each entry: (name, scale, min, max, description)
PARAMS = [
    ('Koc',    'log', 1.0,    1e6,   'Organic carbon partition coeff (L/kg)'),
    ('N',      'lin', 0.70,   1.00,  'Freundlich exponent (-)'),
    ('t_half', 'log', 1.0,    300.0, 'Pesticide half-life at ref conditions (d)'),
    ('VL',     'lin', 2.0,    30.0,  'VFS length (m)'),
    ('DM',     'lin', 0.05,   0.35,  'Initial soil water deficit (-)'),
    ('Ks',     'log', 1e-6,   1e-4,  'Saturated hydraulic conductivity (m/s)'),
]

OCP = 1.0   # %OC — fixed at 1% as in the abstract

# ── Latin Hypercube Sampling ─────────────────────────────────────────────────
def lhs(n_params, n_samples, seed=42):
    rng = np.random.default_rng(seed)
    result = np.zeros((n_samples, n_params))
    for i in range(n_params):
        perm = rng.permutation(n_samples)
        result[:, i] = (perm + rng.random(n_samples)) / n_samples
    return result

def sample_to_params(unit_sample):
    """Convert unit hypercube sample [0,1]^n to physical parameter values."""
    vals = {}
    for i, (name, scale, lo, hi, _) in enumerate(PARAMS):
        u = unit_sample[i]
        if scale == 'log':
            vals[name] = 10**(math.log10(lo) + u*(math.log10(hi)-math.log10(lo)))
        else:
            vals[name] = lo + u*(hi - lo)
    return vals

# ── IWQ file manipulation ────────────────────────────────────────────────────
def read_iwq(path):
    with open(path) as f:
        return f.readlines()

def write_iwq_linear(lines, Koc, OCP, t_half, dgPin, mresn0):
    """Write IWQ for linear isotherm (IKD=1)."""
    out = []
    data_count = 0
    for line in lines:
        s = line.strip()
        if not s or s.startswith('-'):
            out.append(line)
            continue
        nums = [x for x in re.split(r'\s+', s.split(';')[0].strip()) if x]
        if data_count == 0:   # line 1: IWQPRO
            out.append(line)
            data_count += 1
        elif data_count == 1: # line 2: IKD
            out.append(f'1  {Koc:.4f}  {OCP:.1f}                    ; IKD=1 Koc(%OC)\n')
            data_count += 1
        elif data_count == 2: # line 3: CCP
            out.append(line)
            data_count += 1
        elif data_count == 3: # line 4: IDG
            out.append(line)
            data_count += 1
        elif data_count == 4: # line 5: ndgday dgHalf FC dgPin dgML dgLD dgmres0
            parts = nums
            dghalf = math.log(2)/math.log(2)*t_half  # = t_half directly
            out.append(f'{parts[0]}  {t_half:.4f}  {parts[2]}  {dgPin:.6f}  {parts[4]}  {parts[5]}  {mresn0:.6E}  ; ndgday tHalf FC dgPin dgML dgLD dgmres0\n')
            data_count += 1
        else:
            out.append(line)
    return out

def write_iwq_freundlich(lines, Kf, N, t_half, dgPin, mresn0):
    """Write IWQ for Freundlich isotherm (IKD=2). Kf=Kd at C_ref=1 mg/L."""
    out = []
    data_count = 0
    for line in lines:
        s = line.strip()
        if not s or s.startswith('-'):
            out.append(line)
            continue
        nums = [x for x in re.split(r'\s+', s.split(';')[0].strip()) if x]
        if data_count == 0:
            out.append(line)
            data_count += 1
        elif data_count == 1:
            out.append(f'2  {Kf:.6f}  {N:.4f}                ; IKD=2 Freundlich Kf N\n')
            data_count += 1
        elif data_count == 2:
            out.append(line)
            data_count += 1
        elif data_count == 3:
            out.append(line)
            data_count += 1
        elif data_count == 4:
            parts = nums
            out.append(f'{parts[0]}  {t_half:.4f}  {parts[2]}  {dgPin:.6f}  {parts[4]}  {parts[5]}  {mresn0:.6E}  ; ndgday tHalf FC dgPin dgML dgLD dgmres0\n')
            data_count += 1
        else:
            out.append(line)
    return out

def write_iso(src_path, dst_path, Ks, DM, theta_i_override=None):
    """Write ISO file with updated Ks and DM (theta_i = theta_s - DM).

    If theta_i_override is provided it is used directly as theta_i instead of
    theta_s - DM.  Used for events 2+ where OI must equal the final
    inter-event theta from the previous event's OWQ (not the sampled DM).
    """
    with open(src_path) as f:
        lines = f.readlines()
    # ISO format: Ks  Sav  theta_s  theta_i  Sm  schk
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith('_') or s.startswith('-') or s.startswith('K'):
            continue
        nums = re.findall(r'[\d.eE+\-]+', s)
        if len(nums) >= 4:
            theta_s = float(nums[2])
            if theta_i_override is not None:
                theta_i = max(0.025, min(float(theta_i_override), theta_s))
            else:
                theta_i = max(0.01, theta_s - DM)
            # Sav stays as original
            Sav = nums[1]
            Sm  = nums[4] if len(nums)>4 else '0.0'
            schk= nums[5] if len(nums)>5 else '1'
            lines[i] = f'{Ks:.4E}  {Sav}  {theta_s}  {theta_i:.5f}  {Sm}  {schk}\n'
            break
    with open(dst_path, 'w') as f:
        f.writelines(lines)


def parse_owq_final_theta(owq_path):
    """Extract the LAST soil-moisture theta value from the inter-event
    degradation table in an OWQ file (= DGTHETAN used by VFSMOD).
    This must be used as OI for the next event's ISO file.

    Table row format:   day   T(C)   theta(-)
    """
    if not os.path.exists(owq_path):
        return None
    with open(owq_path) as f:
        text = f.read()
    matches = re.findall(
        r'^\s+\d+\s+[\d.]+\s+([\d.]+)\s*$',
        text, re.MULTILINE
    )
    return float(matches[-1]) if matches else None

def write_ikw(src_path, dst_path, VL):
    """Write IKW file with updated filter length VL."""
    with open(src_path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith('t') or s.startswith('-') or s.startswith('f'):
            continue
        nums = re.findall(r'[\d.eE+\-]+', s)
        if len(nums) == 1:  # fwidth line
            continue
        if len(nums) >= 2 and i > 1:  # vl n thetaw ... line
            # replace VL (first number)
            lines[i] = re.sub(r'^(\s*)([\d.eE+\-]+)',
                              lambda m: m.group(1)+f'{VL:.4f}', line, count=1)
            break
    with open(dst_path, 'w') as f:
        f.writelines(lines)

# ── OWQ parsing ──────────────────────────────────────────────────────────────
def parse_owq(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        text = f.read()
    pats = {
        'dP':    r'([\d.]+)\s+%\s+=\s+Pesticide reduction',
        'dQ':    r'([\d.]+)\s+%\s+=\s+Infiltration',
        'dE':    r'([\d.]+)\s+%\s+=\s+Sediment reduction',
        'mresn': r'([\d.E+\-]+)\s+mg/m2=\s+Next event residue',
        'mfml':  r'Mixing layer total mass \(mfml mg\)=\s*([\d.]+)',
        'R':     r'Retardation factor \(R\)=\s*([\d.]+)',
        'C1':    r'Pulse concentration \(C1\)=\s*([\d.]+)',
        'mi':    r'([\d.E+\-]+)\s+mg/m2=\s+Pesticide input',
        'mo':    r'([\d.E+\-]+)\s+mg/m2=\s+Pesticide output',
    }
    res = {}
    for k, p in pats.items():
        m = re.search(p, text)
        res[k] = float(m.group(1)) if m else None
    return res

def parse_osm(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        text = f.read()
    pats = {
        'VIN':  r'Volume from up-field=\s*([\d.E+\-]+)',
        'VOUT': r'Volume from outflow =\s*([\d.E+\-]+)',
        'VF':   r'Volume infiltrated  =\s*([\d.E+\-]+)',
        'TE':   r'Trapping efficiency =\s*([\d.]+)',
    }
    res = {}
    for k, p in pats.items():
        m = re.search(p, text)
        res[k] = float(m.group(1)) if m else None
    return res

# ── Single 3-event run ───────────────────────────────────────────────────────
def run_3events(sample_id, pvals, isotherm, work_dir, base_lines, n_events=3):
    """Run n_events sequential events. Returns dict of results from event n_events."""
    Koc    = pvals['Koc']
    N      = pvals['N']
    t_half = pvals['t_half']
    VL     = pvals['VL']
    DM     = pvals['DM']
    Ks     = pvals['Ks']
    Kd     = Koc * OCP * 0.01
    Kf     = Kd   # equivalent at C_ref=1 mg/L

    run_dir = os.path.join(work_dir, f'gsa_runs', f's{sample_id}_{isotherm}')
    os.makedirs(run_dir, exist_ok=True)

    # Read base project to get file names
    with open(os.path.join(work_dir, f'{BASE_PRJ}.prj')) as f:
        prj_lines = f.readlines()
    file_map = {}
    for line in prj_lines:
        if '=' in line:
            code, path = line.strip().split('=',1)
            file_map[code.strip()] = path.strip()

    mresn       = 0.0
    dgPin_base  = None
    final_theta = None   # last inter-event theta from previous OWQ → OI for next ISO
    results_final = {}

    for ev in range(1, n_events+1):
        ev_dir = os.path.join(run_dir, f'ev{ev}')
        inp_dir = os.path.join(ev_dir, 'inputs')
        out_dir = os.path.join(ev_dir, 'output')
        os.makedirs(inp_dir, exist_ok=True)
        os.makedirs(out_dir, exist_ok=True)

        # Copy and modify input files
        src_iwq = os.path.join(work_dir, file_map.get('iwq','inputs/sampleP.iwq'))
        src_iso = os.path.join(work_dir, file_map.get('iso','inputs/sample.iso'))
        src_ikw = os.path.join(work_dir, file_map.get('ikw','inputs/sampleP.ikw'))

        # Get base dgPin from original IWQ line 5
        if dgPin_base is None:
            with open(src_iwq) as f:
                raw = f.readlines()
            dc=0
            for line in raw:
                s=line.strip()
                if not s or s.startswith('-'): continue
                nums=[x for x in re.split(r'\s+',s.split(';')[0].strip()) if x]
                if dc==4 and len(nums)>=4:
                    dgPin_base=float(nums[3]); break
                if dc<5 and nums: dc+=1
            if dgPin_base is None: dgPin_base=100.0

        dgPin_ev = dgPin_base  # base pesticide input per event

        # Write modified IWQ
        dst_iwq = os.path.join(inp_dir, os.path.basename(src_iwq))
        if isotherm == 'linear':
            new_lines = write_iwq_linear(base_lines, Koc, OCP, t_half, dgPin_ev, mresn)
        else:
            new_lines = write_iwq_freundlich(base_lines, Kf, N, t_half, dgPin_ev, mresn)
        with open(dst_iwq, 'w') as f:
            f.writelines(new_lines)

        # Write modified ISO (Ks, DM).
        # For events 2+, override theta_i with the final inter-event theta from
        # the previous OWQ so OI is consistent with how mresn was computed.
        dst_iso = os.path.join(inp_dir, os.path.basename(src_iso))
        write_iso(src_iso, dst_iso, Ks, DM, theta_i_override=final_theta)

        # Write modified IKW (VL) — keep N nodes proportional
        dst_ikw = os.path.join(inp_dir, os.path.basename(src_ikw))
        write_ikw(src_ikw, dst_ikw, VL)

        # Copy remaining input files unchanged
        for code, path in file_map.items():
            if code in ('owq','og1','og2','ohy','osm','osp'): continue
            src = os.path.join(work_dir, path)
            dst_name = os.path.basename(path)
            dst = os.path.join(inp_dir, dst_name)
            if os.path.exists(src) and not os.path.exists(dst):
                shutil.copy2(src, dst)

        # Write event PRJ
        prj_path = os.path.join(ev_dir, f'ev{ev}.prj')
        with open(prj_path,'w') as f:
            for line in prj_lines:
                if '=' in line:
                    code, path = line.strip().split('=',1)
                    code=code.strip(); fname=os.path.basename(path.strip())
                    if path.strip().startswith('output'):
                        f.write(f'{code}=output/{fname}\n')
                    else:
                        f.write(f'{code}=inputs/{fname}\n')
                else:
                    f.write(line)

        # Run VFSM
        res = subprocess.run(
            [os.path.join(work_dir, 'vfsm'), f'ev{ev}.prj'],
            cwd=ev_dir, capture_output=True, text=True, timeout=120
        )

        # Parse outputs
        owq_name = os.path.basename(file_map.get('owq','output/sampleP.owq'))
        osm_name = os.path.basename(file_map.get('osm','output/sampleP.osm'))
        owq_path = os.path.join(out_dir, owq_name)
        owq_res = parse_owq(owq_path)
        osm_res = parse_osm(os.path.join(out_dir, osm_name))

        mresn = owq_res.get('mresn', 0.0) or 0.0
        # Carry forward the final inter-event theta as OI for the next event
        final_theta = parse_owq_final_theta(owq_path)
        results_final = {**owq_res, **osm_res}

    return results_final

# ── Sobol sensitivity (rank-based approximation) ─────────────────────────────
def sobol_rank_approx(X, y, param_names):
    """
    Simple rank-based sensitivity: Spearman correlation squared
    as a proxy for first-order Sobol indices. Fast and robust.
    """
    from scipy.stats import spearmanr
    n = len(y)
    y_arr = np.array(y, dtype=float)
    mask = np.isfinite(y_arr)
    results = []
    for i, name in enumerate(param_names):
        x = X[mask, i]
        yy = y_arr[mask]
        if len(x) < 5:
            results.append((name, 0.0))
            continue
        r, p = spearmanr(x, yy)
        results.append((name, r, r**2, p))
    results.sort(key=lambda x: -abs(x[1]))
    return results

# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples',  type=int, default=64,  help='LHS samples')
    parser.add_argument('--events',   type=int, default=3,   help='Events per run')
    parser.add_argument('--output',   default='gsa_results.csv')
    parser.add_argument('--seed',     type=int, default=42)
    args = parser.parse_args()

    work_dir = os.path.dirname(os.path.abspath(__file__))
    n_params = len(PARAMS)
    param_names = [p[0] for p in PARAMS]

    # Read base IWQ template
    base_prj = os.path.join(work_dir, f'{BASE_PRJ}.prj')
    with open(base_prj) as f:
        for line in f:
            if 'iwq=' in line.lower():
                iwq_path = os.path.join(work_dir, line.strip().split('=',1)[1])
                break
    with open(iwq_path) as f:
        base_lines = f.readlines()

    # LHS sampling
    unit_samples = lhs(n_params, args.samples, seed=args.seed)
    param_sets = [sample_to_params(unit_samples[i]) for i in range(args.samples)]

    # Output containers
    rows = []
    outputs_lin = {'mresn':[], 'dP':[], 'mfml':[], 'R':[]}
    outputs_fre = {'mresn':[], 'dP':[], 'mfml':[], 'R':[]}
    X_lin = np.zeros((args.samples, n_params))
    X_fre = np.zeros((args.samples, n_params))

    print(f'\nVFSMOD GSA — {args.samples} LHS samples × 2 isotherms × {args.events} events')
    print(f'Total VFSMOD runs: {args.samples * 2 * args.events}')
    print(f'{"─"*65}')

    t0 = time.time()
    for i, pvals in enumerate(param_sets):
        Kd = pvals['Koc'] * OCP * 0.01
        print(f'  Sample {i+1:3d}/{args.samples}  Koc={pvals["Koc"]:8.1f}  N={pvals["N"]:.3f}  '
              f't½={pvals["t_half"]:6.1f}d  VL={pvals["VL"]:5.1f}m  ...', end='', flush=True)

        # Linear run
        r_lin = run_3events(i, pvals, 'linear', work_dir, base_lines, args.events)
        # Freundlich run
        r_fre = run_3events(i, pvals, 'freundlich', work_dir, base_lines, args.events)

        mresn_l = r_lin.get('mresn', None)
        mresn_f = r_fre.get('mresn', None)
        ratio = (mresn_f/mresn_l) if (mresn_l and mresn_l>0 and mresn_f) else None
        print(f' mresn lin={mresn_l:.3f}  fre={mresn_f:.3f}  ratio={ratio:.2f}' if ratio else ' ERROR')

        # Store
        for j, name in enumerate(param_names):
            X_lin[i,j] = unit_samples[i,j]
            X_fre[i,j] = unit_samples[i,j]
        for k in outputs_lin:
            outputs_lin[k].append(r_lin.get(k))
            outputs_fre[k].append(r_fre.get(k))

        row = {**{f'p_{n}': pvals[n] for n in param_names},
               'Kd': Kd, 'Kf': Kd,
               **{f'lin_{k}': r_lin.get(k) for k in outputs_lin},
               **{f'fre_{k}': r_fre.get(k) for k in outputs_fre},
               'ratio_mresn': ratio}
        rows.append(row)

    elapsed = time.time()-t0
    print(f'\nCompleted {args.samples} samples in {elapsed:.1f}s ({elapsed/args.samples:.1f}s/sample)')

    # Write CSV
    csv_path = os.path.join(work_dir, args.output)
    if rows:
        with open(csv_path,'w',newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f'Results written to {args.output}')

    # Sensitivity analysis
    try:
        from scipy.stats import spearmanr
        print(f'\n{"═"*65}')
        print('SENSITIVITY ANALYSIS — Spearman rank correlation with mresn')
        print(f'{"═"*65}')
        for label, X, outs in [('Linear', X_lin, outputs_lin),
                                ('Freundlich', X_fre, outputs_fre)]:
            y = np.array(outs['mresn'], dtype=float)
            valid = np.isfinite(y) & (y>0)
            print(f'\n{label} isotherm (n={valid.sum()} valid):')
            print(f'  {"Parameter":<12} {"Spearman r":>12} {"r²":>8} {"p-value":>10}')
            print(f'  {"─"*45}')
            sens = []
            for j, name in enumerate(param_names):
                if name == 'N' and label == 'Linear':
                    continue  # N not used in linear
                r, p = spearmanr(X[valid,j], np.log10(y[valid]+1e-10))
                sens.append((name, r, r**2, p))
            sens.sort(key=lambda x: -abs(x[1]))
            for s in sens:
                bar = '█'*int(abs(s[1])*20)
                sign = '+' if s[1]>0 else '-'
                print(f'  {s[0]:<12} {s[1]:>+12.4f} {s[2]:>8.4f} {s[3]:>10.4f}  {sign}{bar}')
    except ImportError:
        print('scipy not available — skipping sensitivity indices')

    print(f'\n{"═"*65}')
    print('Done.')

if __name__ == '__main__':
    main()
