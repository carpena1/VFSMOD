#!/bin/sh
# check_vfsm.sh - regression check for VFSMOD.
#
# Runs every sample project, then compares each regenerated output file
# against its committed reference (output/<base>_chk.<ext>).  The run
# PASSES when every regenerated file is byte-identical to its _chk
# reference; any differences are reported per file.  Exit status equals
# the number of files that differ (0 = all good).

# -- locate the vfsm executable (handles a stale ./vfsm symlink) -----------
if [ -x ./vfsm ]; then
    VFSM=./vfsm
elif command -v vfsm >/dev/null 2>&1; then
    VFSM=vfsm
elif [ -x executables/vfsm4621_arm64 ]; then
    VFSM=executables/vfsm4621_arm64
elif [ -x executables/vfsm4621.exe ]; then
    VFSM=executables/vfsm4621.exe
else
    echo "ERROR: cannot find a vfsm executable (./vfsm, on PATH, or in executables/)." >&2
    exit 2
fi

projects="sample.prj sampleP.prj sampleP2.prj sampleP3.prj samplePF.prj \
sampleP3nodeg.prj samplePnodeg.prj samplePnoRain.prj sampleWT.prj sampleWTP.prj"

# -- 1. run all projects ---------------------------------------------------
echo "Running ${VFSM} on all sample projects ..."
for prj in $projects; do
    if "$VFSM" "$prj" >/dev/null 2>&1; then
        echo "  ran  $prj"
    else
        echo "  WARN $prj (vfsm returned non-zero)"
    fi
done

# -- 2. compare each regenerated output against its _chk reference ---------
echo
echo "==================== output vs reference (_chk) ===================="
ok=0; differ=0; skip=0
for ref in output/*_chk.o*; do
    [ -e "$ref" ] || { echo "No reference files (output/*_chk.o*) found."; break; }
    base=$(basename "$ref")
    # live counterpart = reference name with the _chk tag removed
    live=$(printf '%s' "$ref" | sed 's/_chk\././')
    if [ ! -f "$live" ]; then
        echo "  SKIP  $base  (no current output - project not run)"
        skip=$((skip + 1))
    elif diff -q "$live" "$ref" >/dev/null 2>&1; then
        echo "  OK    $base"
        ok=$((ok + 1))
    else
        echo "  DIFF  $base"
        diff "$ref" "$live" | sed 's/^/          /'
        differ=$((differ + 1))
    fi
done

# -- 3. summary ------------------------------------------------------------
echo "===================================================================="
echo "Summary:  $ok identical, $differ differing, $skip skipped"
if [ "$differ" -eq 0 ]; then
    echo "RESULT: PASS - all regenerated outputs match their references"
else
    echo "RESULT: FAIL - $differ file(s) differ from reference"
fi

exit "$differ"
