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
  std::string  xyz_file, basis;
  uint max_iterations = 20;
  double tolerance = 1e-6;
  bool diis_enabled = true;
  int charge = 0;
  
  try {
    std::map<std::string, spdlog::level::level_enum> log_levels = {
      {"debug", spdlog::level::debug},
      {"info", spdlog::level::info},
      {"warn", spdlog::level::warn},
      {"error", spdlog::level::err},
      {"trace", spdlog::level::trace}
    };
    cxxopts::Options options("program", "Hartree Fock Toy Code");
    options.add_options()
      ("l,log-level", "Set logging level [debug, info]", cxxopts::value<std::string>()->default_value("info"))
      ("f,xyz-file", "XYZ File to use", cxxopts::value<std::string>()->default_value("h2.xyz"))
      ("b,basis", "basis set to use", cxxopts::value<std::string>()->default_value("sto-3g"))
      ("c,charge", "charge of system", cxxopts::value<int>()->default_value("0"))
      ("i,max-iterations", "Max number of SCF iterations", cxxopts::value<unsigned int>()->default_value("20"))
      ("t,tolerance", "Tolerance for dE during SCF", cxxopts::value<double>()->default_value("1e-6"))
      ("d,diis", "DIIS enabled: [true, false]", cxxopts::value<std::string>()->default_value("true"))
      ("h,help", "Print usage");
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      return 0;
    }

    std::string log_level = result["log-level"].as<std::string>();
    xyz_file = result["xyz-file"].as<std::string>();
    basis = result["basis"].as<std::string>();
    max_iterations = result["max-iterations"].as<unsigned int>();
    tolerance = result["tolerance"].as<double>();
    if (result["diis"].as<std::string>() == "false") {
      diis_enabled = false;
    }
    charge = result["charge"].as<int>();

    auto it = log_levels.find(log_level);
    if (it != log_levels.end()) {
      spdlog::set_level(it->second);
    } else {
      spdlog::error("Unknown logging level: {}", log_level);
      return 1;
    }
  } catch (const cxxopts::exceptions::exception& e) {
      spdlog::error("Error parsing options: {}", e.what());
      return 1;
  }
  
  spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%5l%$] %v");

  std::ifstream xyz(xyz_file);
  spdlog::info("Using XYZ File: {}", xyz_file);
  if (!xyz) {
    spdlog::error("XYZ File not found: {}", xyz_file);
    return 1;
  }

  libint2::initialize();
  std::vector<libint2::Atom> atoms = libint2::read_dotxyz(xyz);

  spdlog::info("Basis: {}", basis);
  libint2::BasisSet obs(basis, atoms);
  if (obs.size() == 0) {
    spdlog::error("Invalid Basis Set: {}", basis);
    return 1;
  }
  auto S = overlap(obs);
  spdlog::trace("\nOverlap Matrix:\n{}", to_string_mat(S));
  auto T = kinetic(obs);
  spdlog::trace("\nKinetic Matrix:\n{}", to_string_mat(T));
  auto V_ne = electron_nuclear_attraction(obs, atoms);
  spdlog::trace("\nV_ne Matrix:\n{}", to_string_mat(V_ne));
  auto V_ee = electron_electron_repulsion(obs);
  spdlog::trace("\nV_ee Matrix:\n{}", to_string_ten(V_ee));
  auto E_NN = nuclear_nuclear_repulsion_energy(atoms);
  auto molecular_terms = std::make_tuple(S, T, V_ne, V_ee);
  auto scf_parameters = std::make_tuple(tolerance, max_iterations);
  auto electronic_energy = scf_cycle(molecular_terms, scf_parameters, obs, atoms, charge, diis_enabled);
  spdlog::debug("Electronic Energy: {}", electronic_energy);
  spdlog::debug("E_nn Value: {}", E_NN);
  auto total_energy = electronic_energy + E_NN;

  spdlog::info("Total energy: {:.16f}", total_energy);

  libint2::finalize();
  return 0;
}
