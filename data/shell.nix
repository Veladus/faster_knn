let
  pkgs = import <nixpkgs> {};

  nativePackages = with pkgs; [
    python39
  ];

  pythonPackages = with pkgs.python39Packages; [
    geopandas
    pandas
    pyosmium
  ];

in pkgs.mkShell {
  packages = nativePackages ++ pythonPackages;
}
