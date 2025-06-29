# HartreeFockCppProgram

This is a reimplementation of the HartreeFockPythonProgram from
[NickelAndCopper's](https://youtube.com/playlist?list=PL-hN8vfRaU7jSFHN1ZSAMNe_2nXhwAmzM&si=ANjI8kPn-5v_3Kvs)
YouTube Playlist (also, here is the
[GitHub](https://github.com/nickelandcopper/HartreeFockPythonProgram)
for the HatreeFockPythonProgram). Additionally, this code now utilizes the [Libint](https://github.com/evaleev/libint)
library under the hood, so it is necessary to install Libint to be used.


This is more of a testing "toy" code, for me to learn how to use certain libraries as well as learn the
concepts of the Hartree Fock equation. It is in no way performant. All I have really done is made sure
that the values I am getting for the energies make sense when compared to the values from
[NWChem](https://github.com/nwchemgit/nwchem) for the same H2 system and basis set.

Currently the code is set up to use the STO-3G, STO-6G, 3-21G, 6-31G, and 6-311G** basis sets,
based on the values available from the Libint library. Whatever basis sets are available in Libint are
available for use in this code.

## Important Notes
Most of this code is "vibe" coded, meaning that I utilized ChatGPT to help understand logic,
process, and diagnose the code as I wrote it. While the non-Libint implementation of the code
only followed the logic in NickelAndCopper's code, the current implementation of the functions
in `functions.cpp` were less understood. I personally need to better understand the `shells2bf`
functionality of the `libint2::BasisSet` class. ChatGPT directly suggested using the example code
in the Libint Wiki, which is concerning for the license of this current code (as well as other more
important concerns). With that in mind, I have made this version of the code GPL 3.0 licensed.

(Not quite up to speed on licensing so if there are any problems please make and issue and get me up to speed!)

~~Additionally, the code uses a "screening" tool within the `scf_cycle` function within `functions.cpp`. This is
the only way I could get basis sets with l > 0 to successfully run. The impact on adding this code (and why it
is necessary) is beyond my grasp at the moment, so I will come back and update this when I know for sure why
it needs to be added.~~

The screening option was needed due to a bug in the indexing of the one-electron integrals in my code. As far as
I can tell this has been fixed.

The code is not memory friendly, I ran the cc-pVQZ basis set on H2 and saw that the memory usage was ~8GB so be
careful on how big of a basis set you use.

I have also added in the DIIS algorithm (I will not lie, completely copied from ChatGPT. You'll find the
`struct` called `DIISManager` in the `scf.cpp` file.) Does work though when compared with NWChem. The
expected outcomes of using DIIS is less iterations taken to converge the SCF, and making the convergence
process smoother. I highly suggest you try out a larger basis set (6-31g, 6-311gss) with H2O and `-d true` vs
`-d false` to see the difference.

## Build Instructions

Dependencies:
- C++ compiler (I used GCC version 15.1.1)
- cmake
- Eigen
- Boost
- Libint
- spdlog
- cxxopts

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

Available options:
```bash
❯ ./build/program -h                                            
Hartree Fock Toy Code
Usage:
  program [OPTION...]

  -l, --log-level arg       Set logging level [debug, info] (default: info)
  -f, --xyz-file arg        XYZ File to use (default: h2.xyz)
  -b, --basis arg           basis set to use (default: sto-3g)
  -c, --charge arg          charge of system (default: 0)
  -i, --max-iterations arg  Max number of SCF iterations (default: 20)
  -t, --tolerance arg       Tolerance for dE during SCF (default: 1e-6)
  -d, --diis arg            Whether DIIS is enabled or not [true, false] 
                            (default: true)
  -h, --help                Print usage
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

If you want to use `nix run` with the available command line arguments, here is how:
```bash
nix run . -- [-h, --help, -l [debug, info], --log-level [debug, info], etc.]
# Example: nix run . -- --log-level debug
```

## Build Result

The output should be the following:

``` bash
./build/program
[2025-06-29 16:06:38.864] [ info] Using XYZ File: h2.xyz
[2025-06-29 16:06:38.864] [ info] Basis: sto-3g
[2025-06-29 16:06:38.868] [ info] Total energy: -1.0661086493089351
# NWChem Energy (default scf settings):         -1.066108669518
```

# References
Pritchard, Benjamin P., Doaa Altarawy, Brett Didier, Tara D. Gibson, and
Theresa L. Windus. 2019. <span>“New Basis Set Exchange: An Open,
up-to-Date Resource for the Molecular Sciences Community.”</span>
<em>Journal of Chemical Information and Modeling</em> 59 (11): 4814–20.
<a
href="https://doi.org/10.1021/acs.jcim.9b00725">https://doi.org/10.1021/acs.jcim.9b00725</a>.

STO-3G Reference: [Ref](https://www.basissetexchange.org/references/sto-3g/format/txt/?version=1&elements=1)
