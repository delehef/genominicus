{ pkgs, lib, config, inputs, ... }:

{
  packages = [ pkgs.git pkgs.sqlite pkgs.git-cliff pkgs.iconv ];

  languages.rust.enable = true;

  env.LD_LIBRARY_PATH = "${pkgs.sqlite}/lib";
}
