#pragma once
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
  double A = pow((2.0 * alpha / pi), 0.75);
  coord_type coords;

  primitive_gaussian(double alpha, double coeff, coord_type coords)
      : alpha(alpha), coeff(coeff), coords(coords) {}

  coord_type coordinates() const { return coords; }
};

struct atomic_orbital {
  std::vector<primitive_gaussian> primitives;

  atomic_orbital() = default;
  atomic_orbital(std::initializer_list<primitive_gaussian> prims)
      : primitives(prims) {}

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
  std::vector<atom> atoms;

  molecule(std::initializer_list<atom> atms) : atoms(atms) {}

  void push_back(const atom &atm) { atoms.push_back(atm); }

  size_t size() const { return atoms.size(); }

  const atom &operator[](size_t index) const { return atoms[index]; }
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
