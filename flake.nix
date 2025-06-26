{
  inputs.nixpkgs.url = "github.com:nixos/nixpkgs";
  outputs = {self, nixpkgs, ...}:
    let
      system = builtins.currentSystem;
      pkgs = import nixpkgs { inherit system; };
    in
    {
      devShells.${system}.default = import ./shell.nix { inherit pkgs; };
    };
}
