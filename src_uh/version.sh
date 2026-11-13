#!/bin/sh
# version.sh - declare a new UH version in the source banners.
#
# Usage:  sh version.sh vOLD vNEW OLDDATE NEWDATE
# e.g.    sh version.sh version.sh v3.0.9 v3.1.0 04/2025 10/2026
#
# Dates may be plain (05/2026) or backslash-escaped (05\\/2026); the older
# instructions showed the escaped form because the dates used to be handed
# straight to sed.  Both work.
#
# Only BANNER lines are rewritten - the lines that state the version the program
# currently is.  Every other mention of a version is history and must survive,
# or the source ends up claiming that old fixes shipped in the new release.
# The previous version of this script rewrote history too.
#
# It then reports what changed, what it left alone, and any banner missing the
# new date.

set -e

if [ $# -ne 4 ]; then
    sed -n '2,6p' "$0" >&2
    exit 2
fi

OLDV=$1
NEWV=$2
OLDD=`printf '%s' "$3" | tr -d '\\'`
NEWD=`printf '%s' "$4" | tr -d '\\'`

FILES=`ls *.f 2>/dev/null` || true
[ -n "$FILES" ] || { echo "no .f files here - run this from the source directory" >&2; exit 2; }

# A banner states the running program's identity.  Anything else is history.
BANNER='Version:|Last Updated:|\.\.\.FINISHED\.\.\.|File: |@@'

total=0
for f in $FILES; do
    awk -v ov="$OLDV" -v nv="$NEWV" -v od="$OLDD" -v nd="$NEWD" -v ban="$BANNER" '
        function lit(s) { gsub(/[.[\]()*+?^$|\\]/, "\\\\&", s); return s }
        # Replace the old version only where no further version number
        # follows it, so v4.6.3 does not match the start of v4.6.3.1 (and a
        # second run cannot turn v4.6.3.1 into v4.6.3.1.1).
        function subv(s, re, rep,   out, rest) {
            out = ""
            while (match(s, re)) {
                rest = substr(s, RSTART + RLENGTH)
                out = out substr(s, 1, RSTART - 1) (rest ~ /^\.?[0-9]/ ? substr(s, RSTART, RLENGTH) : rep)
                s = rest
            }
            return out s
        }
        BEGIN { ovr = lit(ov); odr = lit(od) }
        {
            before = $0
            if (before ~ ban) { $0 = subv($0, ovr, nv); gsub(odr, nd) }
            if ($0 != before) print FILENAME ":" FNR "  " $0 > "/dev/stderr"
            print
        }
    ' "$f" > "$f.tmp" 2> "$f.log"

    c=`wc -l < "$f.log" | tr -d ' '`
    if [ "$c" -gt 0 ]; then
        mv "$f.tmp" "$f"
        echo "$f: $c banner line(s) updated"
        sed 's/^/    /' "$f.log"
        total=`expr $total + $c`
    else
        rm -f "$f.tmp"
    fi
    rm -f "$f.log"
done

echo ""
echo "$total banner line(s) updated."
echo ""
echo "Left as history, not touched - confirm this is correct:"
# same boundary as the replacement: the old version not followed by more numbers
OLDRE=`printf '%s\n' "$OLDV" | sed 's/[.]/[.]/g'`
grep -nE "${OLDRE}([^.0-9]|[.][^0-9]|[.]?\$)" *.f 2>/dev/null | sed 's/^/    /' || echo "    (none)"
echo ""
echo "Banner lines carrying $NEWV but NOT the new date $NEWD:"
grep -nE "$BANNER" *.f 2>/dev/null | grep -- "$NEWV" | grep -v -- "$NEWD" | sed 's/^/    /' \
    || echo "    (none)"
