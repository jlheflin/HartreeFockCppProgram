#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

struct coord_type {
  double x, y, z;

  coord_type operator-(const coord_type &other) const {
    return {x - other.x, y - other.y, z - other.z};
  }

  coord_type operator+(const coord_type &other) const {
    return {x + other.x, y + other.y, z + other.z};
  }

  coord_type operator/(double &scalar) const {
    return {x / scalar, y / scalar, z / scalar};
  }
};

inline coord_type operator*(double scalar, const coord_type &coord) {
  return {scalar * coord.x, scalar * coord.y, scalar * coord.z};
}

struct primitive_gaussian {
  double alpha = 0;
  double coeff = 0;
  const double pi = M_PI;
  double A = pow((2.0 * alpha / pi), 0.75);

  primitive_gaussian(double alpha, double coeff) : alpha(alpha), coeff(coeff) {}
};

struct atomic_orbital {
  std::vector<primitive_gaussian> primitives;
  coord_type coords;
  double x, y, z;

  atomic_orbital() = default;
  atomic_orbital(std::initializer_list<primitive_gaussian> prims,
                 coord_type coords)
      : primitives(prims), coords(coords) {}

  void push_back(const primitive_gaussian &pg) { primitives.push_back(pg); }

  size_t size() const { return primitives.size(); }

  const primitive_gaussian &operator[](size_t index) const {
    return primitives[index];
  }
};

struct orbital_set {
  std::vector<atomic_orbital> aobs;

  orbital_set() = default;
  orbital_set(std::initializer_list<atomic_orbital> aos) : aobs(aos) {}

  void push_back(const atomic_orbital &ao) { aobs.push_back(ao); }

  size_t size() const { return aobs.size(); }

  const atomic_orbital &operator[](size_t index) const { return aobs[index]; }

  auto begin() const { return aobs.begin(); }
  auto end() const { return aobs.end(); }
};

struct atom {
  std::string symbol;
  int Z;
  double x;
  double y;
  double z;

  atom(std::string symbol, int Z, double x, double y, double z)
      : symbol(symbol), Z(Z), x(x), y(y), z(z) {}
};

struct molecule {
  std::vector<atomic_orbital> aos;
  std::vector<int> Z_list;
  std::vector<coord_type> coord_list;

  molecule(std::initializer_list<atomic_orbital> aos) : aos(aos) {}

  void push_back(const atomic_orbital &ao) { aos.push_back(ao); }

  size_t size() const { return aos.size(); }

  const atomic_orbital &operator[](size_t index) const { return aos[index]; }
};

struct molecular_system {
  using mol_type = std::tuple<molecule, orbital_set>;
  std::vector<mol_type> mol_sys;

  molecular_system() = default;
  molecular_system(std::initializer_list<mol_type> mol_ao) : mol_sys(mol_ao) {}

  void push_back(const mol_type &mols) { mol_sys.push_back(mols); }

  size_t size() const { return mol_sys.size(); }

  const mol_type &operator[](size_t index) const { return mol_sys[index]; }
};

std::ostream &operator<<(std::ostream &os, const primitive_gaussian &pg);
std::ostream &operator<<(std::ostream &os, const atomic_orbital &ao);
std::ostream &operator<<(std::ostream &os, const atom &atm);
std::ostream &operator<<(std::ostream &os, const molecule &mol);
