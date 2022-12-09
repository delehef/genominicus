{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.cargo pkgs.rust-analyzer pkgs.rustc pkgs.rustfmt pkgs.clippy
    pkgs.sqlite
  ];

  LD_LIBRARY_PATH = "${pkgs.sqlite}/lib";
}
