#!/usr/bin/env python3
"""
run_multievent.py — VFSMOD multi-event sequential simulation
Runs N identical storm events, passing dgmres0 (porewater residue) from
each event's OWQ output into the next event's IWQ input.

Usage:
  python3 run_multievent.py <base_prj> <n_events> [--ikd2 Kf N]

Examples:
  python3 run_multievent.py sampleP 3
  python3 run_multievent.py sampleP 3 --ikd2 8.13 0.85

Output:
  Creates event_1/, event_2/, ... subdirectories with full I/O for each event
  Prints a summary table of key results for all events
"""

import sys, os, re, shutil, subprocess, argparse


def parse_owq_mresn(owq_path):
    """Extract mresn (next event remobilization, mg/m2) from OWQ file."""
    with open(owq_path) as f:
        text = f.read()
    m = re.search(r'([\d.E+\-]+)\s+mg/m2=\s+Next event residue', text)
    if m:
        return float(m.group(1))
    return 0.0


def parse_owq_results(owq_path):
    """Extract key results from OWQ file."""
    results = {}
    with open(owq_path) as f:
        text = f.read()
    patterns = {
        'dP (%)':       r'([\d.]+)\s+%\s+=\s+Pesticide reduction',
        'dQ (%)':       r'([\d.]+)\s+%\s+=\s+Infiltration',
        'dE (%)':       r'([\d.]+)\s+%\s+=\s+Sediment reduction',
        'mi (mg/m2)':   r'([\d.E+\-]+)\s+mg/m2=\s+Pesticide input',
        'mo (mg/m2)':   r'([\d.E+\-]+)\s+mg/m2=\s+Pesticide output',
        'mfml (mg)':    r'Mixing layer total mass \(mfml mg\)=\s*([\d.]+)',
        'mresn(mg/m2)': r'([\d.E+\-]+)\s+mg/m2=\s+Next event residue',
        'R':            r'Retardation factor \(R\)=\s*([\d.]+)',
        'C1 (mg/L)':    r'Pulse concentration \(C1\)=\s*([\d.]+)',
    }
    for key, pat in patterns.items():
        m = re.search(pat, text)
        results[key] = float(m.group(1)) if m else None
    return results



def parse_owq_final_theta(owq_path):
    """Extract the LAST soil-moisture theta value from the inter-event
    degradation table in an OWQ file.
    """
    if not owq_path or not os.path.exists(owq_path):
        return None
    with open(owq_path) as f:
        text = f.read()
    matches = re.findall(
        r'^\s+\d+\s+[\d.]+\s+([\d.]+)\s*$',
        text, re.MULTILINE
    )
    return float(matches[-1]) if matches else None

def update_iso_oi(iso_path, new_oi):
    """Update the Initial soil-water content (OI) in the ISO file."""
    with open(iso_path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith('_') or s.startswith('-') or s.startswith('K'):
            continue
        nums = re.findall(r'[\d.eE+\-]+', s)
        if len(nums) >= 4:
            # Reconstruct the line with the new OI (which is the 4th element)
            # The structure is: Ks  Sav  theta_s  theta_i  Sm  schk
            Ks = nums[0]
            Sav = nums[1]
            theta_s = float(nums[2])
            theta_i = min(float(new_oi), theta_s)
            Sm = nums[4] if len(nums) > 4 else '0.0'
            schk = nums[5] if len(nums) > 5 else '1'
            lines[i] = f'{Ks}  {Sav}  {theta_s}  {theta_i:.5f}  {Sm}  {schk}\n'
            break
    with open(iso_path, 'w') as f:
        f.writelines(lines)

def update_iwq_mresn(iwq_path, new_mresn):
    """Update the dgmres0 value in an IWQ file (last number on line 5)."""
    with open(iwq_path) as f:
        lines = f.readlines()
    new_lines = lines.copy()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith('-') or stripped.startswith(';'):
            continue
        nums = re.findall(r'[\d.E+\-]+', stripped.split(';')[0])
        if len(nums) >= 7:  # line 5: ndgday dgHalf FC dgPin dgML dgLD dgmres0 (7 values)
            new_val = f'{new_mresn:.6E}'
            # Replace last numeric token before any comment
            new_line = re.sub(
                r'([\d.E+\-]+)(\s*(?:;.*)?)\s*$',
                f'{new_val}\\2',
                line.rstrip(), count=1
            ) + '\n'
            new_lines[i] = new_line
            break
    with open(iwq_path, 'w') as f:
        f.writelines(new_lines)


def set_dgpin_zero(iwq_path):
    """Set dgPin (4th number on line 5) to 0.0 — single-application events 2+."""
    with open(iwq_path) as f:
        lines = f.readlines()
    new_lines = lines.copy()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith('-') or stripped.startswith(';'):
            continue
        nums = re.findall(r'[\d.E+\-]+', stripped.split(';')[0])
        if len(nums) >= 7:  # ndgday dgHalf FC dgPin dgML dgLD dgmres0
            count = [0]
            def repl(m, _c=count):
                _c[0] += 1
                return '0.000000E+00' if _c[0] == 4 else m.group(0)
            if ';' in line:
                code, rest = line.rstrip().split(';', 1)
                new_lines[i] = re.sub(r'[\d.E+\-]+', repl, code) + ';' + rest + '\n'
            else:
                new_lines[i] = re.sub(r'[\d.E+\-]+', repl, line.rstrip()) + '\n'
            break
    with open(iwq_path, 'w') as f:
        f.writelines(new_lines)


def update_iwq_ikd2(iwq_path, kf, n):
    """Replace IKD line with Freundlich parameters."""
    with open(iwq_path) as f:
        lines = f.readlines()
    new_lines = lines.copy()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith('-'):
            continue
        # IKD line: first token is 0, 1, or 2 (second data line, after IWQPRO)
        parts = stripped.split()
        if parts[0] in ('0', '1', '2') and i > 0:
            new_lines[i] = (
                f'2  {kf:.4f}  {n:.4f}'
                f'                    ; IKD=2 Freundlich: Kf(L^N/Kg) N(-)\n'
            )
            break
    with open(iwq_path, 'w') as f:
        f.writelines(new_lines)


def run_event(event_num, base_prj, work_dir, vfsm_exe, kf=None, n=None, single_app=False):
    """Set up and run one event in its own subdirectory."""
    event_dir = os.path.join(work_dir, f'event_{event_num}')
    os.makedirs(event_dir, exist_ok=True)
    inp_dir = os.path.join(event_dir, 'inputs')
    out_dir = os.path.join(event_dir, 'output')
    os.makedirs(inp_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    # Read base project file to get input file list
    with open(os.path.join(work_dir, base_prj + '.prj')) as f:
        prj_lines = f.readlines()

    # Copy all input files to event inputs/
    file_map = {}
    for line in prj_lines:
        if '=' in line:
            code, path = line.strip().split('=', 1)
            code = code.strip()
            path = path.strip()
            fname = os.path.basename(path)
            src = os.path.join(work_dir, path)
            if os.path.exists(src) and not path.startswith('output'):
                dst = os.path.join(inp_dir, fname)
                shutil.copy2(src, dst)
                file_map[code] = fname

    # Write event project file
    prj_path = os.path.join(event_dir, f'event_{event_num}.prj')
    with open(prj_path, 'w') as f:
        for line in prj_lines:
            if '=' in line:
                code, path = line.strip().split('=', 1)
                code = code.strip()
                path = path.strip()
                fname = os.path.basename(path)
                if path.startswith('output'):
                    newpath = f'output/{fname}'
                else:
                    newpath = f'inputs/{fname}'
                f.write(f'{code}={newpath}\n')
            else:
                f.write(line)

    # Apply Freundlich if requested
    iwq_file = os.path.join(inp_dir, file_map.get('iwq', ''))
    if kf is not None and os.path.exists(iwq_file):
        update_iwq_ikd2(iwq_file, kf, n)

    # Single-app: events 2+ get dgPin=0 (only residue remobilization)
    if single_app and event_num > 1 and os.path.exists(iwq_file):
        set_dgpin_zero(iwq_file)

    # Run vfsm
    result = subprocess.run(
        [vfsm_exe, f'event_{event_num}.prj'],
        cwd=event_dir,
        capture_output=True, text=True
    )

    # Find OWQ output path
    owq_path = None
    for line in prj_lines:
        if line.strip().startswith('owq='):
            fname = os.path.basename(line.strip().split('=')[1])
            owq_path = os.path.join(out_dir, fname)
            break

    return owq_path, iwq_file


def main():
    parser = argparse.ArgumentParser(description='VFSMOD multi-event runner')
    parser.add_argument('base_prj', help='Base project name (without .prj)')
    parser.add_argument('n_events', type=int, help='Number of events to simulate')
    parser.add_argument('--ikd2', nargs=2, type=float, metavar=('Kf', 'N'),
                        help='Use Freundlich isotherm: Kf (L^N/kg) and N (-)')
    parser.add_argument('--single-app', action='store_true',
                        help='Apply pesticide only in event 1; events 2+ remobilization only')
    args = parser.parse_args()

    work_dir = os.path.dirname(os.path.abspath(__file__))
    vfsm_exe = os.path.join(work_dir, 'vfsm')

    kf = args.ikd2[0] if args.ikd2 else None
    n  = args.ikd2[1] if args.ikd2 else None

    isotherm = f'Freundlich (Kf={kf}, N={n})' if kf else 'Linear (Kd from .iwq)'
    app_mode = 'Single-application (dgPin=0 for events 2+)' if args.single_app else 'Repeated application'
    print(f'\nVFSMOD Multi-Event Simulation')
    print(f'Base project : {args.base_prj}.prj')
    print(f'Events       : {args.n_events}')
    print(f'Isotherm     : {isotherm}')
    print(f'App. mode    : {app_mode}')
    print(f'{"─"*80}')

    mresn = 0.0  # initial residue for event 1
    final_theta = None
    all_results = []

    for ev in range(1, args.n_events + 1):
        print(f'\nRunning Event {ev}  (dgmres0 = {mresn:.6f} mg/m2)...')

        owq_path, iwq_path = run_event(
            ev, args.base_prj, work_dir, vfsm_exe, kf, n,
            single_app=args.single_app
        )

        # For events 2+, update dgmres0 and OI from previous event and re-run
        if ev > 1:
            needs_rerun = False
            if mresn > 0 and os.path.exists(iwq_path):
                update_iwq_mresn(iwq_path, mresn)
                needs_rerun = True
            
            # Find the iso file in the event's inputs dir
            event_dir = os.path.join(work_dir, f'event_{ev}')
            iso_file = None
            for fname in os.listdir(os.path.join(event_dir, 'inputs')):
                if fname.endswith('.iso'):
                    iso_file = os.path.join(event_dir, 'inputs', fname)
                    break
            
            if final_theta is not None and iso_file and os.path.exists(iso_file):
                update_iso_oi(iso_file, final_theta)
                needs_rerun = True

            if needs_rerun:
                subprocess.run(
                    [vfsm_exe, f'event_{ev}.prj'],
                    cwd=event_dir, capture_output=True, text=True
                )

        if owq_path and os.path.exists(owq_path):
            results = parse_owq_results(owq_path)
            mresn = results.get('mresn(mg/m2)', 0.0) or 0.0
            final_theta = parse_owq_final_theta(owq_path)
            all_results.append(results)
            print(f'  dP={results.get("dP (%)"):.2f}%  '
                  f'mfml={results.get("mfml (mg)"):.1f} mg  '
                  f'mresn={mresn:.4f} mg/m2  '
                  f'R={results.get("R"):.2f}')
        else:
            print(f'  ERROR: OWQ output not found at {owq_path}')
            all_results.append({})
            mresn = 0.0

    # Summary table
    print(f'\n{"═"*80}')
    print(f'SUMMARY — {args.base_prj} — {isotherm}')
    print(f'{"═"*80}')
    keys = ['mi (mg/m2)', 'mo (mg/m2)', 'dP (%)', 'dQ (%)', 'dE (%)',
            'R', 'C1 (mg/L)', 'mfml (mg)', 'mresn(mg/m2)']
    print(f'{"Metric":<20}' + ''.join(f'{"Event "+str(i+1):>15}' for i in range(len(all_results))))
    print('─'*80)
    for key in keys:
        row = f'{key:<20}'
        for r in all_results:
            val = r.get(key)
            row += f'{val:>15.4f}' if val is not None else f'{"N/A":>15}'
        print(row)
    print(f'{"═"*80}\n')


if __name__ == '__main__':
    main()
