{
  description = "Atomic Potentials";

  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.nur.url = "github:nix-community/NUR";
  inputs.nixpkgs.url = "github:nixos/nixpkgs/master";
  inputs.operon.url = "github:heal-research/operon/cpp20";
  inputs.pratt-parser.url = "github:foolnotion/pratt-parser-calculator?rev=a15528b1a9acfe6adefeb41334bce43bdb8d578c";
  inputs.vstat.url = "github:heal-research/vstat/cpp20-eve";

  outputs = { self, flake-utils, nixpkgs, nur, operon, pratt-parser, vstat }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ nur.overlay ];
        };
        repo = pkgs.nur.repos.foolnotion;
      in rec {
        defaultPackage = pkgs.gcc11Stdenv.mkDerivation {
          name = "atomic-potentials";
          src = self;

          hardeningDisable = [ "all" ];
          impureUseNativeOptimizations = true;
          nativeBuildInputs = with pkgs; [
            bear
            cmake
            clang_13
            clang-tools
            cppcheck
          ];
          buildInputs = with pkgs; [
            # Project dependencies and utils for profiling and debugging
            ceres-solver
            cxxopts
            diff-so-fancy
            doctest
            eigen
            fmt
            gdb
            glog
            hotspot
            hyperfine
            jemalloc
            linuxPackages.perf
            mimalloc
            ninja
            openlibm
            operon.defaultPackage.${system}
            pratt-parser.defaultPackage.${system}
            vstat.defaultPackage.${system}
            pkg-config
            valgrind
            xxHash

            # Some dependencies are provided by a NUR repo
            repo.aria-csv
            repo.cmake-init
            repo.cmaketools
            repo.cpp-sort
            repo.fast_float
            repo.robin-hood-hashing
            repo.scnlib
            repo.taskflow
            repo.eve
          ];

          shellHook = ''
            LD_LIBRARY_PATH=${
              pkgs.lib.makeLibraryPath [ pkgs.gcc11Stdenv.cc.cc.lib ]
            };
          '';
        };

        devShell = pkgs.gcc11Stdenv.mkDerivation {
          name = "atomic-dev";
          hardeningDisable = [ "all" ];
          impureUseNativeOptimizations = true;
          nativeBuildInputs = with pkgs; [
            bear
            cmake
            clang_14
            clang-tools
            cppcheck
            include-what-you-use
          ];
          buildInputs = defaultPackage.buildInputs ++ (with pkgs; [
            gdb
            hotspot
            hyperfine
            valgrind
            jemalloc
            linuxPackages.perf
            gnuplot
          ]);

          shellHook = ''
            LD_LIBRARY_PATH=${
              pkgs.lib.makeLibraryPath [ pkgs.gcc11Stdenv.cc.cc.lib ]
            };
          '';
        };
      });
}
