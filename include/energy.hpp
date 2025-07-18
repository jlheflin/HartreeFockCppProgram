#pragma once
#include <Eigen/Dense>
#include <libint2.hpp>
#include <vector>

// using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;

double
nuclear_nuclear_repulsion_energy(const std::vector<libint2::Atom> &atoms);

std::array<double, 3> compute_electronic_energy_expectation_value(const matrix2d &dens_mat,
                                                   const matrix2d &T, const matrix2d &Vne,
                                                   const matrix2d &G);
