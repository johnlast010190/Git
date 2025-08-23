# Repository Guidelines

## Contribution Checklist
- Source env: `source etc/bashrc`.
- Build changed targets: `wmake` in `src/<lib>` or `applications/...`; full build with `./Allwmake`.
- Tests: `cd test/<area> && ./Allrun` (clean with `./Allclean`) and verify pass.
- Tutorials: run a minimal relevant case under `tutorials/...` and capture logs.
- Style: 4-space indent; one class per `.H/.C`; use `Foam::`.
- Docs: include dictionary examples when interfaces/keywords change.
- Commit: present-tense subject (<=72 chars); reference issues.
- PR: describe scope, affected components, build logs, and test/tutorial evidence.

## Project Structure & Module Organization
- src: Core C++ libraries (Foam types, utilities). Each class has a `.H` header and `.C` source.
- applications: Solvers and utilities. Build per subdir with `wmake`.
- tutorials: Runnable cases demonstrating solvers/utilities; also used as functional tests.
- test: Targeted tests grouped by domain. `Allrun`/`Allclean` symlink to tutorials helpers.
- wmake: Build system scripts and rules. Do not modify lightly.
- etc, bin, doc: Environment, helper scripts, and documentation.

## Build, Test, and Development Commands
- Source env: `source etc/bashrc` (sets `WM_*/FOAM_*`).
- Full build: `./Allwmake -j$(nproc)` (repeat to resolve order-dependent builds).
- Component build: `cd src/<lib> && wmake` or `cd applications/<path> && wmake`.
- Clean: `wclean` (in a target) or `wclean all`.
- Run tutorial: `cd tutorials/<area>/<case> && ./Allrun`.
- Debug build: `export WM_COMPILE_OPTION=Debug && ./Allwmake`.

## Coding Style & Naming Conventions
- Language: Modern C++ (as used by OpenFOAM). Indent 4 spaces, no tabs.
- Files: One class per pair `.H`/`.C`; include guards in headers; keep headers minimal.
- Names: Types and classes UpperCamelCase; variables lowerCamelCase; constants UPPER_CASE.
- Namespace: Use `Foam::`; avoid `using namespace` in headers.
- Match nearby code; prefer existing utilities over new dependencies.

## Testing Guidelines
- Location: Add focused tests under `test/<area>`; use existing patterns in that folder.
- Tutorials as tests: Prefer verifying behavior by running an appropriate `tutorials` case via `./Allrun`.
- Conventions: Mirror source layout; name cases and dictionaries descriptively.
- Run: `cd test/<area> && ./Allrun` (clean with `./Allclean`).

## Commit & Pull Request Guidelines
- Commits: Present tense, concise subject (≤72 chars), detailed body with rationale.
- Reference issues: `Fixes #123` or `Refs #123` when applicable.
- PRs must: Describe change, link issues, list affected components, include build logs, and show a passing tutorial/test (`test/...` or `tutorials/...`).
- Include performance/accuracy notes and any dictionary changes with examples.

## Security & Configuration Tips
- Never commit credentials or machine-specific paths; use dictionaries under `system/` and `constant/` for case config.
- Keep environment local: set options via `etc/prefs.sh` rather than editing `etc/bashrc`.
