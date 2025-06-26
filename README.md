# HartreeFockC++Program

This is a reimplementation of the HartreeFockPythonProgram from
[NickelAndCopper's](https://youtube.com/playlist?list=PL-hN8vfRaU7jSFHN1ZSAMNe_2nXhwAmzM&si=ANjI8kPn-5v_3Kvs)
YouTube Playlist (also, here is the [GitHub](https://github.com/nickelandcopper/HartreeFockPythonProgram) for the
HatreeFockPythonProgram).

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
git clone https::/github.com/jlheflin/HartreeFockC++Program.git
```

Configure the build:
```bash
cd ./HatreeFockC++Program
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

The output should be the following:
```bash
Total energy: -1.078522581996274
```
