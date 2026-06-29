# Remove Unused Paper Cleanup Candidates

Generated on 2026-06-26 from `docs/Candidates2Remove.txt`.

## Goal

Remove exactly the tracked files left in `docs/Candidates2Remove.txt`, and nothing more.

The plan intentionally avoids shell glob deletion. The only expanded group is `Data/matlab/*.txt`, and that expansion is performed through `git ls-files` so only tracked files can be removed.

## Preconditions

- Run from the repository root: `/Users/4470246/Repositories/Gemcitabine-model`.
- Review `docs/Candidates2Remove.txt` before running.
- Expect unrelated dirty/untracked files in the worktree; do not stage them.

## Commands

```bash
set -euo pipefail

# 0. Abort if the index already contains unrelated staged work.
git diff --cached --quiet || { echo "Index has pre-existing staged changes; abort."; exit 1; }

# 1. Build an explicit path list from the reviewed candidates.
# B_4/B_5/B_6 are covered by Data/matlab/*.txt, so sort -u removes duplicates.
cleanup_paths="$(mktemp "${TMPDIR:-/tmp}/gemcitabine_cleanup_paths.XXXXXX")"

{
  printf '%s\n' \
    'Code/GemcitabineModel_imagingPipeline_DiagramR.R' \
    'Code/Gemcitabine_Metabolomics_Heatmap/Pathway_enrichment_heatmap_Gemcitabine_2fold.py' \
    'Code/wassersteinFun/ws_distance.m' \
    'Code/Gemcitabine_TrackingResultsAnalysis_V2.R' \
    'Code/cellCyclePseudotime.R' \
    'Code/Manager_combined.m' \
    'Code/combined_cost.m' \
    'Code/combined_ODE.m' \
    'Code/modelDiagram.R' \
    'Code/svmFeatures.RObj' \
    'Code/svm_A2_InterDivision.RObj' \
    'Code/svm_B2_InterDivision.RObj' \
    'Code/svm_C2_InterDivision.RObj' \
    'Code/svm_D2_InterDivision.RObj' \
    'Code/svm_E2_InterDivision.RObj' \
    'Code/svm_F2_InterDivision.RObj' \
    'Code/svm_G2_InterDivision.RObj' \
    'Code/solutions.mat' \
    'Figs/ODE_model.png' \
    'Figs/ccmodel.png'
  git ls-files -- 'Data/matlab/*.txt'
} | sort -u > "$cleanup_paths"

# 2. Review the exact removal set.
cat "$cleanup_paths"
wc -l "$cleanup_paths"
test "$(wc -l < "$cleanup_paths" | tr -d ' ')" -eq 100

# 3. Stage only those tracked removals.
git --literal-pathspecs rm --pathspec-from-file="$cleanup_paths"

# 4. Verify staged changes are only intended deletions.
git diff --cached --name-status
git diff --cached --name-only --diff-filter=D | sort > "${cleanup_paths}.staged"
diff -u "$cleanup_paths" "${cleanup_paths}.staged"
test -z "$(git diff --cached --name-only --diff-filter=ACMRTUXB)"

# 5. Commit the cleanup.
git commit -m "Remove unused paper cleanup candidates"
```

## Expected Removal Count

The command should remove:

- 20 explicitly listed files.
- 80 tracked files expanded from `Data/matlab/*.txt`.

Expected total: 100 staged `D` entries.

One explicit file, `Code/Gemcitabine_Metabolomics_Heatmap/Pathway_enrichment_heatmap_Gemcitabine_2fold.py`, may already be deleted in the worktree before running this plan. It is still tracked, so the expected result remains 100 staged deletions, not necessarily 100 newly removed files from disk.

## Safety Checks

Before committing, `git diff --cached --name-status` should show only `D` entries and only paths from `$cleanup_paths`.

If it shows any `A`, `M`, or unexpected path, stop and unstage before proceeding:

```bash
git restore --source=HEAD --staged --worktree --pathspec-from-file="$cleanup_paths"
```
