#include <classes.hpp>
#include <functions.hpp>

int main() {

  auto sto_3g_h_basis_alpha =
      std::vector{0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00};
  auto sto_3g_h_basis_coeff =
      std::vector{0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00};

  atomic_orbital h1_s, h2_s;

  for (int i = 0; i < 3; i++) {
    auto pg_1 =
        primitive_gaussian(sto_3g_h_basis_alpha[i], sto_3g_h_basis_coeff[i]);
    auto pg_2 =
        primitive_gaussian(sto_3g_h_basis_alpha[i], sto_3g_h_basis_coeff[i]);
    h1_s.push_back(pg_1);
    h2_s.push_back(pg_2);
  }
  h1_s.coords = {0., 0., 0.};
  h2_s.coords = {0., 0., 1.};

  molecule mol = {h1_s, h2_s};
  mol.Z_list = {1, 1};
  mol.coord_list = {h1_s.coords, h2_s.coords};

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
  // print_2d_matrix(S);
  // std::cout << "Kinetic Matrix:\n";
  // print_2d_matrix(T);
  // std::cout << "V_ne Matrix:\n";
  // print_2d_matrix(V_ne);
  // std::cout << "V_ee Matrix:\n";
  // print_tensor4d(V_ee);
  // std::cout << "E_NN Value: " << E_NN << std::endl;
  // std::cout << "Electronic energy: " << electronic_energy << std::endl;
  std::cout << std::setprecision(16) << "Total energy: " << total_energy
            << std::endl;

  return 0;
}
