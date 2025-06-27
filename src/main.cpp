#include <classes.hpp>
#include <Eigen/Eigen>
#include <functions.hpp>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

int main() {

  json basis;

  std::ifstream file("6-31g.1.json");

  if (!file.is_open()){
      throw std::runtime_error("Failed to open basis file!");
  }

  file >> basis;

  std::vector<coord_type> coord_list;
  coord_list.push_back({0., 0., 0.});
  coord_list.push_back({0., 0., 1.});

  auto h_ao_basis =  ao_basis_from_file(basis, coord_list);

  molecule mol;
  for (const atomic_orbital& ao : h_ao_basis) {
      mol.push_back(ao);
  }
  mol.Z_list = {1, 1};

  for (const coord_type& coord : coord_list) {
      mol.coord_list.push_back(coord);
  }

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
  std::cout << std::setprecision(17) << "Total energy: " << total_energy
            << std::endl;

  return 0;
}
