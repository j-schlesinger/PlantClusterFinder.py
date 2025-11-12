{
  description = "DevShell for Python";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      rec {
        devShells.default =
          pkgs.mkShell {
            buildInputs = with pkgs; [
              #python-bin.beta.latest.default
              llvmPackages.clang
              llvmPackages.libclang
              pkg-config
              #(python-bin.selectLatestNightlyWith ( toolchain: toolchain.default.override { extensions = [ "python-analyzer" ]; } ))
              blast
            ];

            shellHook = ''
              export LIBCLANG_PATH="{$pkgs.llvmPackages.libclang}/lib";
              fish
              set -x LIBCLANG_PATH "{$pkgs.llvmPackages.libclang}/lib";
            '';
          };
      }
    );
}
