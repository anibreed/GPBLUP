# Packaging and Release

This repository can be packaged into a tar.gz archive while excluding
example or sensitive data files.

## Quick Start

```bash
chmod +x tools/package_release.sh
./tools/package_release.sh
```

The archive is written to dist/ by default.

## Exclusions

The packaging script excludes:
- build artifacts (build/, build-*)
- binaries (bin/)
- object and module files (*.o, *.mod)
- local data files marked as example or sensitive

Adjust exclusions in tools/package_release.sh if you add new data files.
