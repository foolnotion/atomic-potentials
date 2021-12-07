#include "lib.hpp"

auto main() -> int
{
  library lib;

  return lib.name == "atomic-potentials" ? 0 : 1;
}
