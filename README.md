# How to build

```bash
cd build
cmake ..
make
```

# How to load in FEBioStudio

In FEBioStudio, click `FEBio -> Manage plugins` and load the dynamic library in `build/lib/`.

# FAQ

* Got `symbol not found in flat namespace` when loading the plugin?
    ** Did you forget `DECLARE_FECORE_CLASS();`?
    ** Did you define an un-implemented virtual function but forget `func() = 0;` ?