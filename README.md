# HartreeFockCppProgram

This is a reimplementation of the HartreeFockPythonProgram from
[NickelAndCopper's](https://youtube.com/playlist?list=PL-hN8vfRaU7jSFHN1ZSAMNe_2nXhwAmzM&si=ANjI8kPn-5v_3Kvs)
YouTube Playlist (also, here is the [GitHub](https://github.com/nickelandcopper/HartreeFockPythonProgram) for the
HatreeFockPythonProgram).

Currently the code is set up with the 6-31G basis set for Hydrogen, based on the values
available from the [Basis Set Exhange](www.basissetexchange.org/basis/6-31g/format/json/?version=1&elements=1) [Ref](https://www.basissetexchange.org/references/6-31g/format/txt/?version=1&elements=1)

## Build Instructions

Dependencies:
- C++ compiler (I used GCC version 15.1.1)
- cmake
- Eigen
- Boost

I recently used this repo to learn how to use vcpkg, so there is now a `vcpkg.json` manifest
available. I was able to install the Eigen package listed in `vcpkg.json` by running `vcpkg install`
in the repo directory. There are other actions necessary for CMake and vcpkg to work together, so
please read the vcpkg docs to ensure these steps are taken. If you already have Eigen installed on
your system, the regular compilation should still work.

Clone the repo:
```bash
git clone https::/github.com/jlheflin/HartreeFockCppProgram.git
```


### CMake Build
Configure the build:
```bash
cd ./HatreeFockCppProgram
cmake -S . -B build
```

Compile the code:
```bash
cmake --build build
```

Run the program:
```bash
./build/program
```
### Nix Build
I recently added a flake.nix to this repo to test out how to use Nix Flakes.
The neat thing is that with the `flake.nix` and `flake.lock` files, you don't
have to go out and grab the packages yourself to build the exectuable. You can
try out this build setup by having the Nix package manager installed and running
the following:

Build the program and have the executable in `./result/bin/`:
```bash
cd ./HartreeFockCppProgram
nix build
```

Run the program directly:
```bash
cd ./HartreeFockCppProgram
nix run
```


The output should be the following:
```bash
Total energy: -1.078522581996274
```
