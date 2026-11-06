# Contributing to Piazzi

## Environment
```sh
cmake -S cpp -B cpp/build -G Ninja && cmake --build cpp/build
```
Or open the repository in the devcontainer (`.devcontainer/`), which installs everything.

## Tests
```sh
ctest --test-dir cpp/build --output-on-failure
bash cpp/tests/accuracy.sh                    # positions compared against JPL Horizons
```
Algorithm changes must maintain (or improve) the accuracy measured against JPL Horizons.

## Independent implementation
Piazzi uses its own algorithms, based on published literature and described in [`DESIGN.md`](DESIGN.md).
**Do not contribute code, constants or details obtained by decompiling or analyzing the executable of any other
software.** Compatibility with other programs' file formats (for example, reading a `.cfg`) is welcome; copying their
implementation is not. Cite the reference (paper, book, public specification) for every new algorithm.

## Style
- C++: `clang-format -i` with `cpp/.clang-format`; `cpp/.clang-tidy` warnings must be addressed.
- Comments and UI messages in Portuguese.
- External processes (curl, open) always with an argument list, never a shell string.

## Commits
We use [Conventional Commits](https://www.conventionalcommits.org/pt-br/): `feat:`, `fix:`, `refactor:`, `test:`,
`docs:`, `build:`, `ci:`, `chore:`. Versions follow SemVer, with `vX.Y.Z` tags and an entry in `CHANGELOG.md`.

## Pull requests
1. One topic per PR, with tests covering the change.
2. CI (`.github/workflows/ci.yml`) must pass: formatting, build, tests, accuracy and Docker.
3. No third-party files without a redistribution license.
