# HartreeFockCppProgram

This is a reimplementation of the HartreeFockPythonProgram from
[NickelAndCopper's](https://youtube.com/playlist?list=PL-hN8vfRaU7jSFHN1ZSAMNe_2nXhwAmzM&si=ANjI8kPn-5v_3Kvs)
YouTube Playlist (also, here is the
[GitHub](https://github.com/nickelandcopper/HartreeFockPythonProgram)
for the HatreeFockPythonProgram). Additionally, this code now utilizes the [Libint](https://github.com/evaleev/libint)
library under the hood, so it is necessary to install Libint to be used.


This is more of a testing "toy" code, for me to learn how to use certain libraries as well as learn the
concepts of the Hartree Fock equation. It is in no way performant. All I have really done is made sure
that the values I am getting for the energies make sense when compared to the values from NWChem for
the same H2 system and basis set.

Currently the code is set up with the STO-3G, STO-6G, 3-21G, 6-31G, and 6-311G** basis set for Hydrogen,
based on the values available from the Libint library. Whatever basis sets are available in Libint are
available for use in this code.

STO-3G, STO-6G Reference: [Ref](https://www.basissetexchange.org/references/sto-6g/format/txt/?version=1&elements=1)  
3-21G Reference: [Ref](https://www.basissetexchange.org/references/3-21g/format/txt/?version=1&elements=1)  
6-31G Reference: [Ref](https://www.basissetexchange.org/references/6-31g/format/txt/?version=1&elements=1)  
6-311G Reference: [Ref](https://www.basissetexchange.org/references/6-311g/format/txt/?version=0&elements=1)

## Important Notes
Most of this code is "vibe" coded, meaning that I utilized ChatGPT to help understand logic,
process, and diagnose the code as I wrote it. While the non-Libint implementation of the code
only followed the logic in NickelAndCopper's code, the current implementation of the functions
in `functions.cpp` were less understood. I personally need to better understand the `shells2bf`
functionality of the `libint2::BasisSet` class. ChatGPT directly suggested using the example code
in the Libint Wiki, which is concerning for the license of this current code (as well as other more
important concerns). With that in mind, I have made this version of the code GPL 3.0 licensed.

(Not quite up to speed on licensing so if there are any problems please make and issue and get me up to speed!)

Additionally, the code uses a "screening" tool within the `scf_cycle` function within `functions.cpp`. This is
the only way I could get basis sets with l > 0 to successfully run. The impact on adding this code (and why it
is necessary) is beyond my grasp at the moment, so I will come back and update this when I know for sure why
it needs to be added.

## Build Instructions

Dependencies:
- C++ compiler (I used GCC version 15.1.1)
- cmake
- Eigen
- Boost
- Libint

I recently used this repo to learn how to use vcpkg, so there is now a
`vcpkg.json` manifest available. I was able to install the Eigen package
listed in `vcpkg.json` by running `vcpkg install` in the repo directory.
There are other actions necessary for CMake and vcpkg to work together,
so please read the vcpkg docs to ensure these steps are taken. If you
already have Eigen installed on your system, the regular compilation
should still work. NOTE: This will not install Libint.

Clone the repo:

``` bash
git clone https://github.com/jlheflin/HartreeFockCppProgram.git
```

### CMake Build

Configure the build:

``` bash
cd ./HatreeFockCppProgram
cmake -S . -B build -DCMAKE_PREFIX_PATH=<path-to-libint-install>
```

Compile the code:

``` bash
cmake --build build
```

Run the program:

``` bash
./build/program
```

### Nix Build

I recently added a flake.nix to this repo to test out how to use Nix
Flakes. The neat thing is that with the `flake.nix` and `flake.lock`
files, you don't have to go out and grab the packages yourself to build
the exectuable. You can try out this build setup by having the Nix
package manager installed and running the following:

Build the program and have the executable in `./result/bin/`:

``` bash
cd ./HartreeFockCppProgram
nix build
./result/bin/program
```

Run the program directly:

``` bash
cd ./HartreeFockCppProgram
nix run
```

## Build Result

The output should be the following:

``` bash
Basis: sto-3g
Total energy: -1.0661086493089351
# NWChem E:   -1.066108669518

Basis: sto-6g
Total energy: -1.0735829307648752
# NWChem E:   -1.073582951298

Basis: 3-21g
Total energy: -1.0913860702729483
# NWChem E:   -1.091386084615


Basis: 6-31g
Total energy: -1.0948079614680448
# NWChem E:   -1.094807976031


Basis: 6-311gss
Total energy: -1.1015899868969445
# NWChem E:   -1.101590002337
```

# References
Pritchard, Benjamin P., Doaa Altarawy, Brett Didier, Tara D. Gibson, and
Theresa L. Windus. 2019. <span>“New Basis Set Exchange: An Open,
up-to-Date Resource for the Molecular Sciences Community.”</span>
<em>Journal of Chemical Information and Modeling</em> 59 (11): 4814–20.
<a
href="https://doi.org/10.1021/acs.jcim.9b00725">https://doi.org/10.1021/acs.jcim.9b00725</a>.
