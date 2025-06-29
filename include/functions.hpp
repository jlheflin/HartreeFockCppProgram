#pragma once
#include <unsupported/Eigen/CXX11//Tensor>
#include <libint2.hpp>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;

int test();


matrix2d overlap(libint2::BasisSet obs);
matrix2d kinetic(libint2::BasisSet obs);
matrix2d electron_nuclear_attraction(libint2::BasisSet obs, std::vector<libint2::Atom> atoms);

tensor4d electron_electron_repulsion(libint2::BasisSet obs);
double
nuclear_nuclear_repulsion_energy(std::vector<libint2::Atom> atoms);
matrix2d 
compute_density_matrix(const matrix2d  &mos, int n_occ);

matrix2d compute_G(matrix2d  dens_mat, tensor4d Vee);

double compute_electronic_energy_expectation_value(matrix2d dens_mat,
                                                   matrix2d T, matrix2d Vne,
                                                   matrix2d G);
double
scf_cycle(std::tuple<matrix2d, matrix2d, matrix2d, tensor4d> molecular_terms,
          std::tuple<double, int> scf_parameters, libint2::BasisSet obs, std::vector<libint2::Atom> atoms, int charge = 0);

