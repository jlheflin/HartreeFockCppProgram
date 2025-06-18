#include <classes.hpp>
#include <functions.hpp>
#include <iostream>

int main() {

  auto sto_3g_h_basis_alpha =
      std::vector{0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00};
  auto sto_3g_h_basis_coeff =
      std::vector{0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00};

  atom h1("H", 1, 0, 0, 0);
  atom h2("H", 1, 0, 0, 1);

  molecule mol({h1, h2});
  mol.push_back(h1);
  std::cout << mol[0] << std::endl;
  std::cout << mol << std::endl;

  atomic_orbital h1_s, h2_s;

  for (int i = 0; i < 3; i++) {
    auto pg_1 =
        primitive_gaussian(sto_3g_h_basis_alpha[i], sto_3g_h_basis_coeff[i],
                           std::vector{std::make_tuple(h1.x, h1.y, h1.z)});
    auto pg_2 =
        primitive_gaussian(sto_3g_h_basis_alpha[i], sto_3g_h_basis_coeff[i],
                           std::vector{std::make_tuple(h2.x, h2.y, h2.z)});
    h1_s.push_back(pg_1);
    h2_s.push_back(pg_2);
  }

  molecular_system h2_sys;
  orbital_set h2_aobs;

  h2_aobs.push_back(h1_s);
  h2_aobs.push_back(h2_s);

  h2_sys.push_back(std::make_tuple(mol, h2_aobs));

  std::cout << "Running overlap" << std::endl;
  auto overlap_matrix = overlap(h2_aobs);
  std::cout << "Overlap done" << std::endl;

  for (int i = 0; i < overlap_matrix.size(); i++) {
    for (int j = 0; j < overlap_matrix.size(); j++) {
      std::cout << overlap_matrix[i][j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << overlap_matrix[0][0] << std::endl;
  std::cout << h1_s << std::endl;
  std::cout << h2_s << std::endl;
}
