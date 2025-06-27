#include <Eigen/Eigen>
#include <classes.hpp>
#include <fstream>
#include <functions.hpp>
#include <iomanip>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main() {

  std::array<std::string, 5> files = {"sto-3g.1.json", "sto-6g.1.json",
                                      "3-21g.1.json", "6-31g.1.json",
                                      "6-311g.0.json"};

  std::vector<coord_type> coord_list;
  coord_list.push_back({0., 0., 0.});
  coord_list.push_back({0., 0., 1.});

  std::map<std::string, std::vector<atomic_orbital>> basis_sets;
  std::map<std::string, molecule> molecules;

  for(auto filename : files) {
    json json_data;
    std::ifstream file(filename);
    file >> json_data;
    file.close();

    basis_sets[filename] = ao_basis_from_file(json_data, coord_list);

    for (atomic_orbital ao : basis_sets[filename]) {
      molecules[filename].push_back(ao);
    }
    for (coord_type coord : coord_list) {
      molecules[filename].coord_list.push_back(coord);
    }

    molecules[filename].Z_list = {1, 1};


  }

  for (std::string filename : files) {
    auto mol = molecules[filename];
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
  std::cout << std::setprecision(17) << "Basis: " << filename << "\nTotal energy: " << total_energy << "\n"
            << std::endl;
  }
  return 0;
}
