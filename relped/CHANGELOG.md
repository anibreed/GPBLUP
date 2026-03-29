# Changelog — relped

## 2026-03-26
- chore: Removed Colleau debug-target API and diagnostic prints from `colleau_mod.f90`.
- chore: Removed temporary tabular comparison and per-target diagnostic dumps from `renped.f90`.
- chore: Removed generated diagnostic artifacts and run logs (cleanup of committed debug outputs).

Notes:
- Sparse Colleau inbreeding implementation validated against the tabular fallback on small test pedigrees before cleanup.
- If further debugging is needed, use `diag_colleau.f90` (kept as a small diagnostic tool) or re-enable verbose output via `set_colleau_verbose(.true.)`.
