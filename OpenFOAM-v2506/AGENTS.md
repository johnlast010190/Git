# Repository Guidelines

## Contribution Checklist
- Source env: `source etc/bashrc`.
- Build: `./Allwmake -j`; rebuild a component with `wmake`, clean with `wclean`.
- Sanity: run `foamSystemCheck` and, if needed, `foamInstallationTest`.
- Tutorials/tests: `foamTestTutorial -full <case>` or run `tutorials/...` with `Allrun`.
- Scope: keep diffs small; avoid committing local `etc/` or machine-specific changes.
- Style: follow OpenFOAM conventions; headers `.H`, sources `.C`.
- Docs: include dictionary examples when interfaces/keywords change.
- Commit: concise subject; link issues; scope per module.
- PR: include summary, build commands, and tutorial output/logs.

## Project Structure & Module Organization
- Source: `src/` (core libraries, C++ headers `.H` and sources `.C`).
- Applications: `applications/` (solvers, utilities) built with `wmake`.
- Configuration: `etc/` (environment, compiler/toolchain settings).
- Build system: `wmake/` and top-level `Allwmake*` scripts.
- Extensions: `modules/` and `plugins/` (optional add-ons).
- Documentation: `doc/` (build and usage docs), `README.md`.
- Tutorials: `tutorials/` (cases used for validation and smoke tests).

## Build, Test, and Development Commands
- Setup environment: `source <path>/OpenFOAM-v2506/etc/bashrc`.
- Full build: `./Allwmake -j -s -q -l` (parallel, quiet, queued, logged).
- Rebuild a component: `wmake` in a library or application dir; clean with `wclean`.
- System check: `foamSystemCheck`; installation test: `foamInstallationTest`.
- Run a tutorial test: `foamTestTutorial -full incompressible/simpleFoam/pitzDaily`.
- Manual tutorial run (example):
  - `run && cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily ./`
  - `( cd pitzDaily && blockMesh && simpleFoam )`
- Debug build: `export WM_COMPILE_OPTION=Debug && ./Allwmake`.

## Coding Style & Naming Conventions
- Language: C++ (OpenFOAM framework). Indent with 4 spaces; no tabs.
- Files: headers `.H`, sources `.C`; file names match class names.
- Types/classes: `CamelCase` (e.g., `volScalarField`); functions/members lowerCamelCase.
- Dictionaries use standard OpenFOAM layout: `system/`, `constant/`, time directories (`0`, `1`, ...).
- Follow existing patterns in `src/` and `applications/`; prefer small, focused patches.
- Namespace: Use `Foam::`; avoid `using namespace` in headers.

## Testing Guidelines
- Prefer tutorial-based validation in `tutorials/` with `Allrun`/`Allclean`.
- For changes affecting a solver/utility, provide a minimal tutorial case or reference an existing one and include command logs.
- Use `tutorials/Alltest` or `Allrun -test` options for quick sweeps where available.

## Commit & Pull Request Guidelines
- Commit messages: concise subject (<=72 chars), body describing motivation, scope, and affected components (e.g., `src/transportModels: fix viscosity bounds`).
- PRs: include summary, linked issues, build commands used, and tutorial evidence (logs, brief notes, optional screenshots of results).
- Scope PRs per module; avoid unrelated formatting churn. Adhere to `etc/` config conventions and do not commit local environment changes.

## Security & Configuration Tips
- Do not hardcode absolute paths; rely on `$WM_PROJECT_DIR`, `$FOAM_SRC`, `$FOAM_APP`, `$FOAM_TUTORIALS`.
- Keep ThirdParty usage configurable; see `doc/Build.md` for details. Avoid committing binaries or large case data.
- Keep environment local: set options via `etc/prefs.sh` rather than editing `etc/bashrc`.
