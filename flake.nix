{
  description = "A standard-compliant Python project using uv and Nix";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    # points to parent dev env repository
    python-toolchain.url = "github:j-schlesinger/nix-templates?dir=development/python";
  };

  outputs = { self, nixpkgs, flake-utils, python-toolchain }:
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
