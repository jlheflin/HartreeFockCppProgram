#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11//Tensor>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;

matrix2d 
compute_density_matrix(const matrix2d  &mos, int n_occ);

matrix2d compute_G(const matrix2d &dens_mat, const tensor4d &Vee);
