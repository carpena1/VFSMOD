# Contributing to VFSMOD

VFSMOD is a modeling system with two programs:

- `vfsm` - the core model, source in `src_vfsm/`.
- `uh` - unit hydrograph utility for runoff/sediment, source in `src_uh/`.

All development happens on the single trunk branch `main`. The old versioned
branches (`4.6.0`, `4.6.1`, `4.6.2`, ...) and `master` are retired.

## Branch model

- `main` is the only permanent branch. Every merge lands here.
- **No versioned branches.** Never create a branch named after a version
  (`v4.7.0`, `4.6.3`, ...). A version is a tag on `main`, not a branch.
- Feature/bugfix work uses short-lived topic branches that name the work, not
  the target version:
  `git switch -c feature/my-feature`
- Cut releases from `main` (see below). *Never* cut one from a topic branch.

## Submitting a change

1. Branch from an up-to-date `main`:
   ```
   git fetch origin && git switch -c <topic> origin/main
   ```
2. Edit.
3. Build locally before pushing:
   ```
   ./setup            # builds vfsm + uh, symlinks both at repo root
   ```
   gfortran override (the `make release` CI path is static):
   ```
   make -C src_vfsm FC=gfortran LD=gfortran
   ```
4. Run regression checks / tests:
   ```
   ./check_vfsm.sh
   ./check_uh.sh
   ./vfsm sample
   ```
5. Push and open a PR against `main`. The GitHub `Build` workflow runs the
   three-platform (Linux/macOS/Windows) build on the PR.
6. Wait for CI success, then merge.

Do not `git commit` your object files or build outputs; `vfsm`/`uh` at the
repo root are symlinks, not tracked binaries.

## Declaring a new version (before tagging)

Versions are declared on `main` in a single commit, then tagged afterwards. All
version strings in the code itself are updated together (sync, or the release is
inconsistent):

1. On an up-to-date `main`, add the new version's entry to
   `src_vfsm/CHANGES.txt` (the authoritative changelog) under a
   `Changes in X.Y.Z` heading; keeping the same style. `CHANGES.md` at the repo
   root is the Markdown mirror; update it to match.
2. Update the version and release date in the source code locations.
   From `src_vfsm/` and `src_uh/`:
   ```
   sh version.sh vOLD vNEW 'MM\/YYYY' 'MM\/YYYY'
   ```
   e.g. `sh version.sh v4.6.2.1 v4.7.0 '05\/2026' '11\/2026'`
3. Update the global `VERSION` file at the repo root to the new number. e.g. `4.6.2.1`
4. Commit all of it on `main` with a message naming the version
   (e.g. `v4.7.0: <name-or-summary>`). No branch, no PR needed.

## Releasing (tagging)

Releases are triggered from `main` by pushing a **version tag**.

1. Make sure `main` is pushed so CI sees the release state:
   ```
   git push origin main
   ```
2. Tag the release commit and push the tag (include the `v`!):
   ```
   git tag -a v4.7.0 -m "VFSMOD 4.7.0"
   git push origin v4.7.0
   ```

   The `Build` workflow re-runs on the tag commit and produces the
   Linux/macOS/Windows packages. The annotated `vX.Y.Z` tag on `main` is
   the release marker - it records exactly which commit was shipped.

3. To re-release / roll back (although, could use the 4th version number slot 
   for small changes), move the tag to the corrected `main` commit and
   force-push it (NOTE: only use if necessary; versions should not be ambiguous 
   once released!):
   ```
   git tag -f -a v4.7.0 -m "VFSMOD 4.7.0 (fixed)"
   git push -f origin v4.7.0
   ```

### Rules

- The tag must point at a commit on `main`, never at a topic branch.
- Bump Version -> Tag (in that order) on the same commit: version strings should be settled
  on `main`, then the tag added, so the built binaries will carry the declared version.
- One permanent branch, one release trigger (a tag push). There is no
  long-lived maintenance branch; patches go to `main` and the release tag now
  points at the last commit for that particular version.
