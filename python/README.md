# Python packages

Each subdirectory contains an independently buildable Python distribution with its own `pyproject.toml`, build configuration, and package sources. Native extensions consume the installed ALPS C++ SDK.

The current package is [pyalps](pyalps/README.md), which provides Python bindings, simulation helpers, and analysis tools. Its tests live in [`tests/pyalps/`](../tests/pyalps/), with packaging checks in [`tests/packaging/`](../tests/packaging/).

Additional packages belong alongside `pyalps/`, with their own distribution and import names so they can coexist during a migration.
