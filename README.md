# HartreeFockC++Program

This is a reimplementation of the HartreeFockPythonProgram from
[NickelAndCopper's](https://youtube.com/playlist?list=PL-hN8vfRaU7jSFHN1ZSAMNe_2nXhwAmzM&si=ANjI8kPn-5v_3Kvs)
YouTube Playlist (also, here is the [GitHub](https://github.com/nickelandcopper/HartreeFockPythonProgram) for the
HatreeFockPytonProgram).

## Build Instructions

Dependencies:
- C++ compiler (I used GCC version 15.1.1)
- cmake
- Eigen
- Boost

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
Total energy: -1.065999461556561
```
