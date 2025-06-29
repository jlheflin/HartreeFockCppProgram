#include <scf.hpp>
#include <density.hpp>
#include <energy.hpp>
#include <spdlog/spdlog.h>

struct DIISManager {
  std::vector<matrix2d> fock_history;
  std::vector<matrix2d> error_history;
  int max_diis = 6;

  void add(const matrix2d& fock, const matrix2d& error) {
    fock_history.push_back(fock);
    error_history.push_back(error);

    if (fock_history.size() > max_diis) {
      fock_history.erase(fock_history.begin());
      error_history.erase(error_history.begin());
    }
  }

  matrix2d extrapolate() {
    int n = error_history.size();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
    rhs(n) = -1.0;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        B(i, j) = (error_history[i].array() * error_history[j].array()).sum();
      }
      B(i, n) = B(n, i) = -1.0;
    }

    Eigen::VectorXd coeffs = B.fullPivLu().solve(rhs);

    matrix2d fock_new = matrix2d::Zero(fock_history[0].rows(), fock_history[0].cols());
    for (int i = 0; i < n; i++) {
      fock_new += coeffs(i) * fock_history[i];
    }
    return fock_new;
  }

  bool ready() const { return fock_history.size() >= 2; }
};

double
scf_cycle(std::tuple<matrix2d, matrix2d, matrix2d, tensor4d> molecular_terms,
          std::tuple<double, uint> scf_parameters, libint2::BasisSet obs, std::vector<libint2::Atom> atoms, int charge, bool diis_enabled) {
  auto [S, T, Vne, Vee] = molecular_terms;
  auto [tolerance, max_iter] = scf_parameters;
  auto electronic_energy = 0.0;
  bool converged = false;
  DIISManager diis;

  int nbasis_functions = obs.nbf();
  matrix2d dens_mat(nbasis_functions, nbasis_functions);
  dens_mat.setZero();

  int electrons = 0;
  for (const auto& atom : atoms) {
    electrons += atom.atomic_number;
  }
  electrons -= charge;
  int n_occ = electrons / 2;

  Eigen::SelfAdjointEigenSolver<matrix2d> es(S);
  if (es.info() != Eigen::Success) {
    throw std::runtime_error("Man this stuff sucks at S");
  }

  Eigen::VectorXd evals = es.eigenvalues();
  Eigen::MatrixXd evecs = es.eigenvectors();

  for (int i = 0; i < evals.size(); i++) {
      evals[i] = 1.0 / std::sqrt(evals[i]);
  }

  matrix2d S_inv_sqrt = evecs * evals.asDiagonal() * evecs.transpose();

  for (int scf_step = 0; scf_step < max_iter; scf_step++) {
    auto electronic_energy_old = electronic_energy;
    auto G = compute_G(dens_mat, Vee);
    auto F = (T + Vne + G).eval();

    matrix2d err = F * dens_mat * S - S * dens_mat * F;
    if (scf_step > 0) {
      diis.add(F, err);
    }

    if (diis_enabled && diis.ready()) {
      F = diis.extrapolate();
    }

    matrix2d F_unitS = S_inv_sqrt * F * S_inv_sqrt;
    matrix2d F_unitS_eigen = F_unitS;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F_unitS_eigen);

    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Eigenvalue decomp failed");
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    auto mos = S_inv_sqrt * eigenvectors;

    dens_mat = compute_density_matrix(mos, n_occ);

    auto [E_energy, E_one, E_two] =
        compute_electronic_energy_expectation_value(dens_mat, T, Vne, G);
    electronic_energy = E_energy;

    double de = electronic_energy - electronic_energy_old;
    bool using_diis = (diis_enabled && diis.ready());

    spdlog::debug("SCF Step: {}\tDIIS Enabled: {}\tElectronic Energy: {}\tdE: {}", scf_step, using_diis, electronic_energy, de);

    if (std::fabs(de) < tolerance) {
      converged = true;
      spdlog::debug("One-electron Energy: {}", E_one);
      spdlog::debug("Two-electron Energy: {}", E_two);
      break;
    }
  }
  if (!converged) {
    spdlog::warn("SCF did not converge, values are not accurate!");
  }
  return electronic_energy;
}
