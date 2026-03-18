# Changelog

## 0.1.0 

* Added changelog 
* Ran `rustfmt` and `clippy` and fixed all issues 
* Using `std::f64::consts::{LN_2, PI}` instead of approximate values (2 places fixed)
* Python project version is now dynamic to match rust project
* Added placeholder for `tests/fig/` to avoid errors in running tests
* Added `_thor.pyi` for typing functions in Rust dll
* Added `.git-blame-ignore-revs` file
* Added github actions workflows for testing python and rust, updated deploy_docs to use uv