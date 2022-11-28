{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.sqlite
  ];

  LD_LIBRARY_PATH = "${pkgs.sqlite}/lib";
}
