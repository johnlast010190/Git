# CFD Workspace

Quick overview of this multi-repo workspace and pointers to the shared guidelines.

## Layout
- OpenFOAM-13: Upstream core
  - Subproject: [README](OpenFOAM-13/README.md) • [AGENTS](OpenFOAM-13/AGENTS.md)
- OpenFOAM-v2506: Foundation fork with extras (`modules/`, `plugins/`)
  - Subproject: [README](OpenFOAM-v2506/README.md) • [AGENTS](OpenFOAM-v2506/AGENTS.md)
- HELYXcore-4.4.0: CMake-based GUI/tools and runtime post-processing
  - Subproject: [README](HELYXcore-4.4.0/README.md) • [AGENTS](HELYXcore-4.4.0/AGENTS.md)

## Getting Started
- Source exactly one OpenFOAM environment per shell:
  - `source OpenFOAM-13/etc/bashrc` or `source OpenFOAM-v2506/etc/bashrc`
- Build:
  - `cd OpenFOAM-13 && ./Allwmake`
  - `cd OpenFOAM-v2506 && ./Allwmake` (extras: `./Allwmake-modules`, `./Allwmake-plugins`)
  - `cd HELYXcore-4.4.0 && cmake -S . -B build && cmake --build build -j`

## Guidelines
- Development, style, testing, and PR conventions: see [root AGENTS](AGENTS.md).
- Prefer relative paths and `$WM_PROJECT_DIR`/`$FOAM_*` variables; avoid committing binaries or build artifacts.
