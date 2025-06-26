{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  packages = [
    pkgs.boost
    pkgs.eigen3
    pkgs.cmake
  ];
  shellHook = ''
  echo "Welcome to dev hell!"
  '';
}
