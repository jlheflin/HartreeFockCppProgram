#pragma once
#include <cmath>
#include <iostream>
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

struct molecule {
  std::vector<atomic_orbital> aos;
  std::vector<int> Z_list;
  std::vector<coord_type> coord_list;

  molecule() = default;
  molecule(std::initializer_list<atomic_orbital> aos) : aos(aos) {}

  void push_back(const atomic_orbital &ao) { aos.push_back(ao); }

  size_t size() const { return aos.size(); }

  const atomic_orbital &operator[](size_t index) const { return aos[index]; }
};

std::ostream &operator<<(std::ostream &os, const primitive_gaussian &pg);
std::ostream &operator<<(std::ostream &os, const atomic_orbital &ao);
std::ostream &operator<<(std::ostream &os, const molecule &mol);
