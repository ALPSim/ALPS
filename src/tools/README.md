# Command-line tools

Each directory owns its command sources and installation rules. C++ commands
install in the `tools` component; `alps-xml` installs in `xml`. `pconfig` links
`ALPS::utilities`; other C++ commands use `ALPS::alps`.

Parameter tools retain the legacy `alps::Parameters` grammar, job generation and
seeding. CLI integration tests and fixtures live in `tests/cli/`.
