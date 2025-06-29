#include <Eigen/Eigen>
#include <fstream>
#include <functions.hpp>
#include <spdlog/spdlog.h>
#include <cxxopts.hpp>

std::string to_string_mat(const Eigen::MatrixXd& mat) {
  std::stringstream ss;
  ss << mat;
  return ss.str();
}
std::string to_string_ten(const Eigen::Tensor<double, 4>& ten) {
  std::stringstream ss;
  ss << ten;
  return ss.str();
}

int main(int argc, char* argv[]) {
  try {
    cxxopts::Options options("program", "Hartree Fock Toy Code");
    options.add_options()
      ("l,log-level", "Set logging level", cxxopts::value<std::string>()->default_value("info"))
      ("h,help", "Print usage");
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      return 0;
    }

    std::string log_level = result["log-level"].as<std::string>();
    if (log_level == "debug") {
      spdlog::set_level(spdlog::level::debug);
    } else {
      spdlog::set_level(spdlog::level::info);
    }
  } catch (const cxxopts::exceptions::exception& e) {
      spdlog::error("Error parsing options: {}", e.what());
      return 1;
  }
  
  spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%5l%$] %v");

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
    spdlog::info("Basis: {}", basis);
    libint2::BasisSet obs(basis, atoms);
    auto S = overlap(obs);
    spdlog::trace("\nOverlap Matrix:\n{}", to_string_mat(S));
    auto T = kinetic(obs);
    spdlog::trace("\nKinetic Matrix:\n{}", to_string_mat(T));
    auto V_ne = electron_nuclear_attraction(obs, atoms);
    spdlog::trace("\nV_ne Matrix:\n{}", to_string_mat(V_ne));
    auto V_ee = electron_electron_repulsion(obs);
    spdlog::trace("\nV_ee Matrix:\n{}", to_string_ten(V_ee));
    auto E_NN = nuclear_nuclear_repulsion_energy(atoms);
    spdlog::debug("E_nn Value: {}", E_NN);
    auto molecular_terms = std::make_tuple(S, T, V_ne, V_ee);
    auto scf_parameters = std::make_tuple(1e-8, 50);
    auto electronic_energy = scf_cycle(molecular_terms, scf_parameters, obs, atoms);
    spdlog::debug("Electronic Energy: {}", electronic_energy);
    auto total_energy = electronic_energy + E_NN;

    spdlog::info("Total energy: {:.16f}", total_energy);
  }
  libint2::finalize();
  return 0;
}
