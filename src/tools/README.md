# Command-line tools

Each directory owns its command sources and installation rules. C++ commands
install in the `tools` component; `alps-xml` installs in `xml`. `pconfig` links
`ALPS::utilities`; other C++ commands use `ALPS::alps`.

Parameter tools retain the legacy `alps::Parameters` grammar, job generation and
seeding. CLI integration tests and fixtures live in `tests/cli/`.

`archive/` owns `txt2archive`, `xml2archive` and the optional SQLite `archive`
command. `xml/` provides `alps-xml` for plotting, conversion and extraction;
see the [command examples](../../CONTRIBUTING.md#xml-resources-and-tools).

`xml/` also installs the historical command names as aliases of `alps-xml`,
generated from `xml/legacy-alias.in`. It also installs three deprecated helpers
without an `alps-xml` equivalent: `archivecat`, `xslttransform` and
`use_local_stylesheet`. They warn on standard error; remove them once the
deprecation period ends.
