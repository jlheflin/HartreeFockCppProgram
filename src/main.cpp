#include <Eigen/Eigen>
#include <fstream>
#include <functions.hpp>
#include <iomanip>


int main() {
  libint2::initialize();

  std::vector<std::string> basis_sets;
  basis_sets.push_back("sto-3g");
  basis_sets.push_back("sto-6g");
  basis_sets.push_back("3-21g");
  basis_sets.push_back("6-31g");
  basis_sets.push_back("6-311gss");
  // basis_sets.push_back("cc-pvqz");
  // basis_sets.push_back("cc-pv6z");


  std::ifstream xyz("h2.xyz");
  std::vector<libint2::Atom> atoms = libint2::read_dotxyz(xyz);

  for (std::string basis : basis_sets) {
    libint2::BasisSet obs(basis, atoms);
    // obs.set_pure(false);
    auto S = overlap(obs); 
    // std::cout << "Overlap Matrix:\n";
    // std::cout << S << "\n" << std::endl;
    auto T = kinetic(obs);
    // std::cout << "Kinetic Matrix:\n";
    // std::cout << T << "\n" << std::endl;
    auto V_ne = electron_nuclear_attraction(obs, atoms);
    // std::cout << "V_ne Matrix:\n";
    // std::cout << V_ne << "\n" << std::endl;
    auto V_ee = electron_electron_repulsion(obs);
    // std::cout << "V_ee Matrix:\n";
    // std::cout << V_ee << "\n" << std::endl;
    auto E_NN = nuclear_nuclear_repulsion_energy(atoms);
    // std::cout << "E_NN Value: " << E_NN << std::endl;
    auto molecular_terms = std::make_tuple(S, T, V_ne, V_ee);
    auto scf_parameters = std::make_tuple(1e-8, 50);
    auto electronic_energy = scf_cycle(molecular_terms, scf_parameters, obs, atoms);
    // std::cout << "Electronic energy: " << electronic_energy << std::endl;
    auto total_energy = electronic_energy + E_NN;

    std::cout << std::setprecision(17) << "Basis: " << basis << "\nTotal energy: " << total_energy << "\n"
              << std::endl;
  }
  libint2::finalize();
  return 0;
}
