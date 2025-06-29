#include <energy.hpp>

double
nuclear_nuclear_repulsion_energy(std::vector<libint2::Atom> atoms) {

  int natoms = atoms.size();

  double E_NN = 0.;

  for (int i = 0; i < natoms; i++) {
    for (int j = i + 1; j < natoms; j++) {
      auto Rijx = atoms[i].x - atoms[j].x;
      auto Rijy = atoms[i].y - atoms[j].y;
      auto Rijz = atoms[i].z - atoms[j].z;

      auto Rijx2 = Rijx * Rijx;
      auto Rijy2 = Rijy * Rijy;
      auto Rijz2 = Rijz * Rijz;
      auto Rij = std::sqrt(Rijx2 + Rijy2 + Rijz2);
      E_NN += (atoms[i].atomic_number * atoms[j].atomic_number) / Rij;
    }
  }
  return E_NN;
}

std::array<double, 3> compute_electronic_energy_expectation_value(matrix2d dens_mat,
                                                   matrix2d T, matrix2d Vne,
                                                   matrix2d G) {

  auto Hcore = T + Vne;
  auto electronic_energy = 0.0;
  double E_one = 0.0;
  double E_two = 0.0;
  auto nbasis_functions = dens_mat.rows();

  for (int i = 0; i < nbasis_functions; i++) {
    for (int j = 0; j < nbasis_functions; j++) {
      E_one += dens_mat(i,j) * Hcore(i,j);
      E_two += dens_mat(i,j) * 0.5 * G(i,j);
      electronic_energy += dens_mat(i,j) * (Hcore(i,j) + 0.5 * G(i,j));
    }
  }
  return {electronic_energy, E_one, E_two};
}
