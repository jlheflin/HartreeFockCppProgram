#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using coord_type = std::vector<std::tuple<double, double, double>>;

struct primitive_gaussian {
  double alpha = 0;
  double coeff = 0;
  const double pi = M_PI;
  double A = (2.0 * alpha / pi);
  coord_type coords;

  primitive_gaussian(double alpha, double coeff, coord_type coords)
      : alpha(alpha), coeff(coeff), coords(coords) {}

  coord_type coordinates() const { return coords; }
};

std::ostream &operator<<(std::ostream &os, const primitive_gaussian &pg) {
  os << "primitive_gaussian(alpha=" << pg.alpha << ", coeff=" << pg.coeff
     << ", coords={";

  for (const auto &[x, y, z] : pg.coords) {
    os << "(" << x << ", " << y << ", " << z << ")";
  }
  os << "})";
  return os;
}

struct atomic_orbital {
  std::vector<primitive_gaussian> primitives;

  atomic_orbital(std::initializer_list<primitive_gaussian> prims)
      : primitives(prims) {}

  void push_back(const primitive_gaussian &pg) { primitives.push_back(pg); }

  size_t size() const { return primitives.size(); }

  const primitive_gaussian &operator[](size_t index) const {
    return primitives[index];
  }
};

std::ostream &operator<<(std::ostream &os, const atomic_orbital &ao) {
  for (const auto &pg : ao.primitives) {
    os << pg << std::endl;
  }
  return os;
}

struct atom {
  std::string symbol;
  int Z;
  double x;
  double y;
  double z;

  atom(std::string symbol, int Z, double x, double y, double z)
      : symbol(symbol), Z(Z), x(x), y(y), z(z) {}
};

std::ostream &operator<<(std::ostream &os, const atom &atm) {
  os << "Symbol: " << atm.symbol << "\n";
  os << "Z: " << atm.Z << "\n";
  os << "x: " << atm.x << "\n";
  os << "y: " << atm.y << "\n";
  os << "z: " << atm.z << "\n";
  return os;
}

struct molecule {
  std::vector<atom> atoms;

  molecule(std::initializer_list<atom> atms) : atoms(atms) {}

  void push_back(const atom &atm) { atoms.push_back(atm); }

  size_t size() const { return atoms.size(); }

  const atom &operator[](size_t index) const { return atoms[index]; }
};

std::ostream &operator<<(std::ostream &os, const molecule &mol) {
  for (const auto &atm : mol.atoms) {
    os << atm.symbol << " " << atm.x << " " << atm.y << " " << atm.z << "\n";
  }
  return os;
}

int main() {

  atom h1("H", 1, 0.00234235, 0, 0);
  atom h2("H", 1, 0, 0, 1);

  molecule mol({h1, h2});
  mol.push_back(h1);
  std::cout << mol[0] << std::endl;
  std::cout << mol << std::endl;
}
