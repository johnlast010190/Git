# Repository Guidelines

## Contribution Checklist
- Sync branch: `git pull --rebase` then `git checkout -b feature/<name>`.
- Configure/build: `bin/emake` (refresh cache with `bin/emake -r`).
- Targeted build: `bin/emake <lib|app>`; clean with `bin/emake -clean <target>`.
- Tests: enable via `unitTests_REQUIRED=ON`; run `ctest -j$(nproc)` from the build dir.
- Manual check: run an example or app and verify logs/output.
- Docs: update `doc/` or examples when behavior/CLI changes.
- Commit: scoped prefix (`src:`, `apps:`, `cmake:`, `docs:`); link issues.
- PR: include build/test commands and screenshots for UI changes.

## Project Structure & Module Organization
- `src`: Core C++ libraries (OpenFOAM/HELYX-style `.C`/`.H`).
- `applications`: Executables; targets end with `-exe` (e.g., `simpleFoam-exe`).
- `modules`: Optional/feature modules enabled via CMake.
- `thirdParty`, `VTK`: External dependencies and EVTK/VTK integration used by RTPP.
- `bin`: Developer tools (`emake`, `helyxRun`, utilities).
- `etc`: CMake config, user settings, and templates.
- `examples`, `doc`: Example cases and documentation.
- `platforms`: Generated at build; contains `activeBuild.shrc` and per-build env files.

## Build, Test, and Development Commands
- Build current tree: `bin/emake` (run from repo root or any subdir).
- Targeted build: `bin/emake finiteVolume`, `bin/emake simpleFoam-exe`.
- Refresh CMake cache: `bin/emake -r` (after toolchain/settings changes).
- Clean target: `bin/emake -clean <target>` or `bin/emake --clean-first <target>`.
- Package release: `bin/emake package` (CMake packaging rules).
- Activate runtime env: `source platforms/activeBuild.shrc` (PATH, autocompletion).
- Run tests (if enabled): from build dir (e.g., `cbuild/<HELYX_OPTIONS>_<MPI>`), `ctest -j$(nproc)`.

## Coding Style & Naming Conventions
- Language: Modern C++ within `Foam::` namespace; follow OpenFOAM style.
- Indentation: 4 spaces, no tabs; 80‑column guideline; braces on new lines.
- Files: `.C` implementations, `.H` headers; class/file names use PascalCase.
- Targets: Libraries use their base name (e.g., `finiteVolume`); executables end with `-exe`.
- Namespace: Use `Foam::`; avoid `using namespace` in headers.

## Testing Guidelines
- Framework: CTest via CMake. Enable with `unitTests_REQUIRED=ON` (cache or `etc/userSettings.cmake`).
- Layout: Place tests under `unitTests/`, mirroring module paths; register with `add_test(...)`.
- Execute: From the build directory, run `ctest -j$(nproc)`.

## Commit & Pull Request Guidelines
- Commits: Imperative subject; scoped prefixes like `src:`, `apps:`, `cmake:`, `docs:`; reference issues (e.g., `#123`).
- PRs: Include summary, linked issues, modules/targets touched, build steps (`bin/emake ...`), and test notes; attach screenshots for user-facing changes.

## Security & Configuration Tips
- Use `etc/userSettings.cmake` for local overrides; do not commit `platforms/*.shrc` or build artifacts.
- Prefer `HELYX_PROJECT_DIR` and related vars over hard-coded paths.
- For RTPP/EVTK, set `VTK_ARCH_PATH` as needed in `etc/userSettings.cmake`.
- Never commit credentials or machine-specific paths.
