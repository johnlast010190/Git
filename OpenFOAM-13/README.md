# OpenFOAM-13

Quick reference for this subproject. See root [README](../README.md) and [AGENTS](../AGENTS.md) for workspace-wide conventions and policies.

## Quick Start
- Source env: `source etc/bashrc`
- Full build: `./Allwmake -j` (repeat if dependencies resolve on second pass)
- Component build: `cd src/<lib> && wmake` or `cd applications/<path> && wmake`
- Clean: `wclean` (in target) or `wclean all`
- Smoke test: run a small case from `tutorials/` with `./Allrun`

## Notes
- Headers `.H` pair with sources `.C`; use `Foam::` namespace
- Add focused tests under `test/` when modifying libraries/utilities
