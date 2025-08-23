# HELYXcore 4.4.0

Quick reference for this subproject. See root [README](../README.md) and [AGENTS](../AGENTS.md) for workspace-wide conventions and policies.

## Quick Start
- Configure: `cmake -S . -B build`
- Build: `cmake --build build -j`
- Clean: `cmake --build build --target clean`
- Run tests (if present): consult module `tests/` or project docs

## Notes
- CMake-based project; prefer out-of-source builds in `build/`
- Sources under `src/` and `applications/`; keep generated artifacts out of VCS
