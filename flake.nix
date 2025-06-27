{
  description = "Hartree-Fock C++ program";

  inputs.nixpkgs.url = "github:NixOs/nixpkgs/nixos-unstable";

  outputs = {self, nixpkgs, ...}:
    let
      forAllSystems = nixpkgs.lib.genAttrs [ "x86_64-linux" "aarch64-linux" ];
    in
    {
      packages = forAllSystems (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in {
        default = pkgs.stdenv.mkDerivation {
          pname = "program";
          version = "1.0.0";

          src = ./.;
          nativeBuildInputs = [ pkgs.cmake ];
          buildInputs = [ pkgs.boost pkgs.eigen pkgs.nlohmann_json ];

          cmakeFlags = [ "-DCMAKE_INSTALL_PREFIX=${placeholder "out"}" ];

          installPhase = ''
          runHook preInstall
          cmake --install .
          runHook postInstall
          '';
        };
      });
      
      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
        in {
          default = pkgs.mkShell {
            packages = [
              pkgs.boost
              pkgs.eigen
              pkgs.cmake
            ];
            shellHook = ''
              echo "Welcome to dev shell!"
            '';
          };
        });
    };
}
