# OpenFOAM-v2506

Quick reference for this subproject. See root [README](../README.md) and [AGENTS](../AGENTS.md) for workspace-wide conventions and policies.

## Quick Start
- Source env: `source etc/bashrc`
- Full build: `./Allwmake -j -s -q -l`
- Rebuild a component: `wmake` in a library/app dir; clean with `wclean`
- System checks: `foamSystemCheck` (and `foamInstallationTest` if needed)
- Tutorials: `foamTestTutorial -full <case>` or run `tutorials/.../Allrun`

## Notes
- Supports optional `modules/` and `plugins/` alongside standard `src/` and `applications/`
- File pairs: `.H` headers with `.C` sources; prefer small, focused patches
