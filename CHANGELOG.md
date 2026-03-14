# Changelog

## 0.1.1 2026-Mar-11

* Add changelog
* Robustify condition for using parallel implementation in python bindings to handle zero input
* Format python and rust & clear commented blocks
* Fix clippy lints
* Add type stubs for python bindings & add py.typed marker to indicate type hints are present
* Add placeholder file to prevent tests/fig/ folder from being deleted, causing test failures
* Set python package version based on rust crate version
* Parametrize magnetization tests over field method
  * Dipole passing, tet failing for now
* Add rust and python lint/test workflows
* Add dev dep group
* Add doc deps to dev group
* Set minimum numpy dep version to the lowest compatible version
* Use uv install for doc deps
* Suppress interactive plots when running tests
* Use pathlib for test_utils paths
* Add .cargo/config.toml to use x86-64-v3 reference cpu for x86 builds
  * Enables modern vector instructions, among other things
* Use abi3-py310 for bindings, s.t. one binary is compatible with python >=3.10
* Unmask _h_demag_tet4 function
* Add ruff and ty for linting and typechecking
* Add ruff linter configuration in pyproject.toml
* Derive PartialEq on Vec3 and Mat3 structs instead of manually implementing methods
* Use slice `&[T]` instead of vec reference `&Vec<T>`
* Use PI and LN_2 built-in constants instead of baked values
* Only deploy docs when running actions on main branch
* Move python package to `src/thor/` for canonical "src-layout" package format
* Rename `test_utils` module to `testing` to avoid being picked up by pytest
* Update tests to use test functions so that they don't run on import
* Add `test_examples`
* Move examples and python benchmarks to `examples/`
