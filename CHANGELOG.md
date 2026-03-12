# Changelog

## 0.1.1

* Add changelog
* Robustify condition for using parallel implementation in python bindings to handle zero input
* Format & clear commented blocks
* Add type stubs for python bindings
* Add placeholder file to prevent tests/fig/ folder from being deleted, causing test failures
* Set python package version based on rust crate version
* Parametrize magnetization tests over field method
  * Dipole passing, tet failing for now
* Add rust lint/test workflow
* Add dev dep group
* Suppress interactive plots when running tests
* Use pathlib for test_utils paths
