#include <classes.hpp>
#include <iostream>

std::ostream &operator<<(std::ostream &os, const primitive_gaussian &pg) {
  os << "primitive_gaussian(alpha=" << pg.alpha << ", coeff=" << pg.coeff
     << ", coords={";

  for (const auto &[x, y, z] : pg.coords) {
    os << "(" << x << ", " << y << ", " << z << ")";
  }
  os << "})";
  return os;
}

std::ostream &operator<<(std::ostream &os, const atomic_orbital &ao) {
  for (const auto &pg : ao.primitives) {
    os << pg << std::endl;
  }
  return os;
}

std::ostream &operator<<(std::ostream &os, const atom &atm) {
  os << "Symbol: " << atm.symbol << "\n";
  os << "Z: " << atm.Z << "\n";
  os << "x: " << atm.x << "\n";
  os << "y: " << atm.y << "\n";
  os << "z: " << atm.z << "\n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const molecule &mol) {
  for (const auto &atm : mol.atoms) {
    os << atm.symbol << " " << atm.x << " " << atm.y << " " << atm.z << "\n";
  }
  return os;
}
