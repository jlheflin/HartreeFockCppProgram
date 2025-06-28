#pragma once
#include <classes.hpp>
#include <unsupported/Eigen/CXX11//Tensor>
#include <nlohmann/json.hpp>
#include <libint2.hpp>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;
using json = nlohmann::json;

int test();


matrix2d overlap(libint2::BasisSet obs);
matrix2d kinetic(libint2::BasisSet obs);
matrix2d electron_nuclear_attraction(libint2::BasisSet obs, std::vector<libint2::Atom> atoms);

tensor4d electron_electron_repulsion(molecule mol);
double
nuclear_nuclear_repulsion_energy(const std::vector<coord_type> &coord_list,
                                 const std::vector<int> &Z_list);
matrix2d 
compute_density_matrix(const matrix2d  &mos, int n_occ);

matrix2d compute_G(matrix2d  dens_mat, tensor4d Vee);

double compute_electronic_energy_expectation_value(matrix2d dens_mat,
                                                   matrix2d T, matrix2d Vne,
                                                   matrix2d G);
double
scf_cycle(std::tuple<matrix2d, matrix2d, matrix2d, tensor4d> molecular_terms,
          std::tuple<double, int> scf_parameters, molecule mol);

std::vector<atomic_orbital> ao_basis_from_file(json& basis_data, std::vector<coord_type>& coords);
