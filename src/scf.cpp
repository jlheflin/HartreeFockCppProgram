#include <scf.hpp>
#include <density.hpp>
#include <energy.hpp>

double
scf_cycle(std::tuple<matrix2d, matrix2d, matrix2d, tensor4d> molecular_terms,
          std::tuple<double, int> scf_parameters, libint2::BasisSet obs, std::vector<libint2::Atom> atoms, int charge) {
  auto [S, T, Vne, Vee] = molecular_terms;
  auto [tolerance, max_iter] = scf_parameters;
  auto electronic_energy = 0.0;

  int nbasis_functions = obs.nbf();
  matrix2d dens_mat(nbasis_functions, nbasis_functions);
  dens_mat.setZero();

  int electrons = 0;
  for (const auto& atom : atoms) {
    electrons += atom.atomic_number;
  }
  electrons -= charge;
  int n_occ = electrons / 2;

  for (int scf_step = 0; scf_step < max_iter; scf_step++) {
    auto electronic_energy_old = electronic_energy;
    auto G = compute_G(dens_mat, Vee);
    auto F = T + Vne + G;

    Eigen::SelfAdjointEigenSolver<matrix2d> es(S);
    if (es.info() != Eigen::Success) {
      throw std::runtime_error("Man this stuff sucks at S");
    }

    Eigen::VectorXd evals = es.eigenvalues();
    Eigen::MatrixXd evecs = es.eigenvectors();


    double cond_number = es.eigenvalues().maxCoeff() / es.eigenvalues().minCoeff();

    double eps = 1e-12;

    for (int i = 0; i < evals.size(); i++) {
        evals[i] = 1.0 / std::sqrt(evals[i]);
    }
    matrix2d S_inv_sqrt = evecs * evals.asDiagonal() * evecs.transpose();

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

    electronic_energy =
        compute_electronic_energy_expectation_value(dens_mat, T, Vne, G);

    if (std::fabs(electronic_energy - electronic_energy_old) < tolerance) {
      break;
    }
  }
  return electronic_energy;
}
