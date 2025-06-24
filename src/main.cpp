#include <classes.hpp>
#include <Eigen/Eigen>
#include <functions.hpp>
#include <iomanip>

int main() {

  auto six_31g_h_basis_alpha =
      std::vector{0.1873113696E+02, 0.2825394365E+01, 0.6401216923E+00, 0.1612777588E+00};
  auto six_31g_h_basis_coeff =
      std::vector{0.3349460434E-01, 0.2347269535E+00, 0.8137573261E+00, 1.0000000};

  atomic_orbital h1_1s, h2_1s, h1_2s, h2_2s;

  for (int i = 0; i < 3; i++) {
    auto pg_1 =
        primitive_gaussian(six_31g_h_basis_alpha[i], six_31g_h_basis_coeff[i]);
    auto pg_2 =
        primitive_gaussian(six_31g_h_basis_alpha[i], six_31g_h_basis_coeff[i]);
    h1_1s.push_back(pg_1);
    h2_1s.push_back(pg_2);
  }
  for (int i = 3; i < 4; i++) {
    auto pg_1 =
        primitive_gaussian(six_31g_h_basis_alpha[i], six_31g_h_basis_coeff[i]);
    auto pg_2 =
        primitive_gaussian(six_31g_h_basis_alpha[i], six_31g_h_basis_coeff[i]);
    h1_2s.push_back(pg_1);
    h2_2s.push_back(pg_2);
  }
  h1_1s.coords = {0., 0., 0.};
  h2_1s.coords = {0., 0., 1.};
  h1_2s.coords = {0., 0., 0.};
  h2_2s.coords = {0., 0., 1.};

  molecule mol = {h1_1s, h2_1s, h1_2s, h2_2s};
  mol.Z_list = {1, 1};
  mol.coord_list = {h1_1s.coords, h2_1s.coords};

  auto S = overlap(mol);
  auto T = kinetic(mol);
  auto V_ne = electron_nuclear_attraction(mol, mol.Z_list);
  auto V_ee = electron_electron_repulsion(mol);
  auto E_NN = nuclear_nuclear_repulsion_energy(mol.coord_list, mol.Z_list);
  auto molecular_terms = std::make_tuple(S, T, V_ne, V_ee);
  auto scf_parameters = std::make_tuple(1e-5, 20);
  auto electronic_energy = scf_cycle(molecular_terms, scf_parameters, mol);
  auto total_energy = electronic_energy + E_NN;

  // std::cout << "Overlap Matrix:\n";
  // std::cout << S << std::endl;
  // std::cout << "Kinetic Matrix:\n";
  // std::cout << T << std::endl;
  // std::cout << "V_ne Matrix:\n";
  // std::cout << V_ne << std::endl;
  // std::cout << "V_ee Matrix:\n";
  // std::cout << V_ee << std::endl;
  // std::cout << "E_NN Value: " << E_NN << std::endl;
  // std::cout << "Electronic energy: " << electronic_energy << std::endl;
  std::cout << std::setprecision(16) << "Total energy: " << total_energy
            << std::endl;

  return 0;
}
