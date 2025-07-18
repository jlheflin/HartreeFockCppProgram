#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <libint2.hpp>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;

double
scf_cycle(const std::tuple<matrix2d, matrix2d, matrix2d, tensor4d> &molecular_terms,
          const std::tuple<double, uint> &scf_parameters, const libint2::BasisSet &obs, const std::vector<libint2::Atom> &atoms, int charge = 0,
bool diis_enabled = true);
