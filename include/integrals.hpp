#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <libint2.hpp>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;

matrix2d overlap(libint2::BasisSet obs);
matrix2d kinetic(libint2::BasisSet obs);
matrix2d electron_nuclear_attraction(libint2::BasisSet obs, std::vector<libint2::Atom> atoms);

tensor4d electron_electron_repulsion(libint2::BasisSet obs);
