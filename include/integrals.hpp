#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <libint2.hpp>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;

matrix2d overlap(const libint2::BasisSet& obs);
matrix2d kinetic(const libint2::BasisSet& obs);
matrix2d electron_nuclear_attraction(const libint2::BasisSet& obs, const std::vector<libint2::Atom>& atoms);

tensor4d electron_electron_repulsion(const libint2::BasisSet& obs);
