{
  description = "Atomic Potentials";

  inputs = {
    flake-utils.url = "github:numtide/flake-utils";
    nur.url = "github:foolnotion/nur-pkg";
    nixpkgs.url = "github:nixos/nixpkgs/master";

    operon.url = "github:heal-research/operon/cpp20";
    vstat.url = "github:heal-research/vstat/cpp20-eve";
    pratt-parser.url =
      "github:foolnotion/pratt-parser-calculator?rev=a15528b1a9acfe6adefeb41334bce43bdb8d578c";
    pypi-deps-db.url = "github:DavHau/pypi-deps-db";

    mach-nix = {
      url = "github:DavHau/mach-nix";
      #url = "github:bjornfor/mach-nix/adapt-to-make-binary-wrapper";
      inputs = {
        nixpkgs.follows = "nixpkgs";
        flake-utils.follows = "flake-utils";
        pypi-deps-db.follows = "pypi-deps-db";
      };
    };
  };

  outputs = { self, flake-utils, nixpkgs, nur, operon, pratt-parser, vstat
    , pypi-deps-db, mach-nix }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        repo = pkgs.nur.repos.foolnotion;
        mach = mach-nix.lib.${system};

        python-env = mach.mkPython {
          requirements = ''
            colorama
            coloredlogs
            numpy
            optuna
            pandas
            sympy
          '';

          python = "python39";

          ignoreDataOutdated = true;
        };

      in rec {
        defaultPackage = pkgs.gcc12Stdenv.mkDerivation {
          name = "atomic-potentials";
          src = self;

          hardeningDisable = [ "all" ];
          impureUseNativeOptimizations = true;
          nativeBuildInputs = with pkgs; [
            bear
            cmake
            clang_14
            clang-tools
            cppcheck
          ];

          buildInputs = (with pkgs; [
            cmake
            cxxopts
            doctest
            eigen
            fmt_8
            git
            jemalloc
            openlibm
            pkg-config
            xxHash
            taskflow
            operon.defaultPackage.${system}
            pratt-parser.defaultPackage.${system}
            vstat.defaultPackage.${system}
          ]) ++ (with nur.packages.${system}; [
            # Some dependencies are provided by a private repo 
            aria-csv
            cpp-sort
            eve
            fast_float
            robin-hood-hashing
            scnlib
            seer
          ]);
          shellHook = ''
            LD_LIBRARY_PATH=${
              pkgs.lib.makeLibraryPath [ pkgs.gcc12Stdenv.cc.cc.lib ]
            };
          '';
        };

        devShell = pkgs.gcc12Stdenv.mkDerivation {
          name = "atomic-dev";
          hardeningDisable = [ "all" ];
          impureUseNativeOptimizations = true;
          buildInputs = defaultPackage.buildInputs ++ (with pkgs; [
            gdb
            hotspot
            hyperfine
            valgrind
            linuxPackages.perf
            gnuplot
          ]);

          shellHook = ''
            LD_LIBRARY_PATH=${
              pkgs.lib.makeLibraryPath [ pkgs.gcc12Stdenv.cc.cc.lib ]
            };
          '';
        };
      });
}
