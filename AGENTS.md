# Repository Guidelines

## Project Structure & Module Organization
- OpenFOAM-13: upstream core — `applications/`, `src/`, `wmake/`, `tutorials/`, `test/`.
- OpenFOAM-v2506: foundation fork + extras — standard layout plus `modules/`, `plugins/`.
- HELYXcore-4.4.0: CMake-based GUI/tools — `src/`, `applications/` (and optional `thirdParty/`).
- Naming: headers `.H` pair with sources `.C` (e.g., `MySolver.H`/`MySolver.C` defining `MySolver`).

## Build, Test, and Development Commands
- Environment: source exactly one OpenFOAM env per shell
  - `source OpenFOAM-13/etc/bashrc` OR `source OpenFOAM-v2506/etc/bashrc`
- Build
  - `cd OpenFOAM-13 && ./Allwmake`
  - `cd OpenFOAM-v2506 && ./Allwmake` (extras: `./Allwmake-modules`, `./Allwmake-plugins`)
  - `cd HELYXcore-4.4.0 && cmake -S . -B build && cmake --build build -j`
  - OpenFOAM (Debug): `export WM_COMPILE_OPTION=Debug && ./Allwmake`
- Clean
  - OpenFOAM: `wclean` (in target) or `wclean all`; some dirs provide `Allwclean`
  - HELYXcore: `cmake --build build --target clean`
- Quick tutorial check (example)
  - `run && cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily ./ && (cd pitzDaily && blockMesh && simpleFoam)`

## Coding Style & Naming Conventions
- Language: C++; follow upstream project style; indent 4 spaces, no tabs.
- Files: one primary class per `.H/.C`; minimal headers; sorted includes; use `Foam::` namespace (no `using` in headers).
- Names: Classes CamelCase; methods/vars lowerCamelCase; constants UPPER_CASE.

## Testing Guidelines
- OpenFOAM-13: targeted tests under `test/` (`./Allrun`, clean with `./Allclean`).
- Tutorials as validation: prefer small `tutorials/...` cases demonstrating the change; include a brief `README` and logs.
- HELYXcore: build in `build/`; place tests alongside modules or under `tests/` if present.

## Commit & Pull Request Guidelines
- Commits: present-tense subject (≤72 chars). Example: `OpenFOAM-13: fix CFL check in pisoFoam`.
- Body: motivation, scope, affected paths; reference issues (`Fixes #123`).
- PRs: describe subprojects touched, commands used to build/test, and attach tutorial/test evidence (logs, screenshots where relevant). Keep diffs focused.
- Docs: include dictionary examples when interfaces/keywords change.

## Security & Configuration Tips
- Source only one `etc/bashrc` per shell to avoid ABI/path conflicts; use separate shells for 13 vs v2506.
- Prefer relative paths and `$WM_PROJECT_DIR`, `$FOAM_*` vars; do not commit binaries, generated data, or local `build/` artifacts.
- Keep environment local: set options via `etc/prefs.sh` rather than editing `etc/bashrc`.
