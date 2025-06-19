#include <classes.hpp>
#include <iostream>

std::ostream &operator<<(std::ostream &os, const primitive_gaussian &pg) {
  os << "primitive_gaussian(alpha=" << pg.alpha << ", coeff=" << pg.coeff
     << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const atomic_orbital &ao) {
  for (const auto &pg : ao.primitives) {
    os << pg << std::endl;
  }
  return os;
}
