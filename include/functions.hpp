#pragma once
#include <classes.hpp>
#include <iomanip>

using tensor4d = std::vector<std::vector<std::vector<std::vector<double>>>>;
using matrix2d = std::vector<std::vector<double>>;

inline void print_2d_matrix(std::vector<std::vector<double>> mat) {
  for (auto row : mat) {
    for (auto num : row) {
      std::cout << std::setw(10) << std::setprecision(6) << std::fixed << num
                << " ";
    }
    std::cout << "\n";
  }
}

inline void print_tensor4d(tensor4d tensor) {
  for (auto set1 : tensor) {
    for (auto set2 : set1) {
      for (auto set3 : set2) {
        for (auto num : set3) {
          std::cout << std::setw(10) << std::setprecision(6) << std::fixed
                    << num << " ";
        }
      }
      std::cout << "\n";
    }
  }
}

std::vector<std::vector<double>> overlap(molecule mol);
std::vector<std::vector<double>> kinetic(molecule mol);
std::vector<std::vector<double>>
electron_nuclear_attraction(molecule mol, std::vector<int> Z_list);

tensor4d electron_electron_repulsion(molecule mol);
double
nuclear_nuclear_repulsion_energy(const std::vector<coord_type> &coord_list,
                                 const std::vector<int> &Z_list);
std::vector<std::vector<double>>
compute_density_matrix(const std::vector<std::vector<double>> &mos, int n_occ);

std::vector<std::vector<double>>
compute_G(std::vector<std::vector<double>> dens_mat, tensor4d Vee);

double compute_electronic_energy_expectation_value(matrix2d dens_mat,
                                                   matrix2d T, matrix2d Vne,
                                                   matrix2d G);
double
scf_cycle(std::tuple<matrix2d, matrix2d, matrix2d, tensor4d> molecular_terms,
          std::tuple<double, int> scf_parameters, molecule mol);
