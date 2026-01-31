{
  description = "A standard-compliant Python project using uv and Nix";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    # points to parent dev env repository -- has toolchain/linting built in
    python-toolchain.url = "github:j-schlesinger/nix-templates?dir=development/python";
    mcl-repo.url = "github:j-schlesinger/mcl_flake";
  };

  outputs = { self, nixpkgs, flake-utils, python-toolchain, mcl-repo }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        mkPythonShell = python-toolchain.lib.${system}.mkPythonShell;
      in
      {
        devShells.default = mkPythonShell {
          # Python version as specified by pyproject.toml
          pythonPackage = pkgs.python312;

          # any non-python packages needed as required by this project
          extraPkgs = with pkgs; [
            zlib
            blast-bin
            mcl-repo.packages.${system}.default
          ];

          # synchronize the environment on entry
          extraShellHook = ''
            # Ensure the .venv is up to date with pyproject.toml
            if [ -f "pyproject.toml" ]; then
              echo "Syncing Python dependencies with uv..."
              uv sync
            fi
          '';
        };
      });
}
