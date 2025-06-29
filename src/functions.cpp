#include "functions.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <classes.hpp>
#include <unsupported/Eigen/MatrixFunctions>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;
using json = nlohmann::json;

std::vector<atomic_orbital> ao_basis_from_file(json& basis_data, std::vector<coord_type>& coords) {
  std::vector<atomic_orbital> aos;
  json shells = basis_data["elements"]["1"]["electron_shells"];

  for(auto& coord : coords) {
    for(auto& entry : shells) {
      auto num_pg = entry["coefficients"][0].size();
      atomic_orbital ao;
      for(int i = 0; i < num_pg; i++) {
        double coeff = std::stod(entry["coefficients"][0][i].get<std::string>());
        double alpha = std::stod(entry["exponents"][i].get<std::string>());
        primitive_gaussian pg(alpha, coeff);
        ao.push_back(pg);
      }
      ao.coords = coord;
      aos.push_back(ao);
    }
  }
  return aos;
}

double dot(const coord_type &a, const coord_type &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}


matrix2d overlap(libint2::BasisSet obs) {
  auto nbasis = obs.nbf();
  matrix2d S;
  S.setZero(nbasis, nbasis);

  libint2::Engine s_engine(libint2::Operator::overlap, obs.max_nprim(), obs.max_l());
  // auto shell2bf = obs.shell2bf();
  const auto& buf_vec = s_engine.results();
  const auto shell2bf = obs.shell2bf();
  for (int i = 0; i < obs.size(); i++) {
    int bf1 = shell2bf[i];
    int n1 = obs[i].size();
    for (int j = 0; j < obs.size(); j++) {
      int bf2 = shell2bf[j];
      int n2 = obs[j].size();
      s_engine.compute(obs[i], obs[j]);
      auto ints_shellset = buf_vec[0];
      if (ints_shellset == nullptr)
        continue;


      for (int k = 0; k < n1; k++) {
        for (int l = 0; l < n2; l++) {
          S(bf1 + k, bf2 + l) += ints_shellset[k * n2 * l];
          // std::cout << "S(" << bf1+k << "," << bf2+l << "): " << ints_shellset[k * n2 * l] << std::endl;
        }
      }
    }
  }
  return S;
}

matrix2d kinetic(libint2::BasisSet obs) {
  int nbasis = obs.nbf();

  matrix2d  T(nbasis, nbasis);
  T.setZero();

  libint2::Engine t_engine(libint2::Operator::kinetic, obs.max_nprim(), obs.max_l());
  // auto shell2bf = obs.shell2bf();
  const auto& buf_vec = t_engine.results();
  const auto shell2bf = obs.shell2bf();
  for (int i = 0; i < obs.size(); i++) {
    int bf1 = shell2bf[i];
    int n1 = obs[i].size();
    for (int j = 0; j < obs.size(); j++) {
      int bf2 = shell2bf[j];
      int n2 = obs[j].size();
      t_engine.compute(obs[i], obs[j]);
      auto ints_shellset = buf_vec[0];
      if (ints_shellset == nullptr)
        continue;

      for (int k = 0; k < n1; k++) {
        for (int l = 0; l < n2; l++) {
          T(bf1 + k, bf2 + k) += ints_shellset[k * n2 * l];
        }
      }
    }
  }
  return T;
}

matrix2d electron_nuclear_attraction(libint2::BasisSet obs, std::vector<libint2::Atom> atoms) {
  int nbasis = obs.nbf();
  matrix2d V_ne(nbasis,nbasis);
  V_ne.setZero();

  std::vector<std::pair<libint2::scalar_type, std::array<libint2::scalar_type, 3>>> q;
  for (const auto& atom : atoms) {
    q.push_back({static_cast<libint2::scalar_type>(atom.atomic_number), {{atom.x, atom.y, atom.z}}});
  }

  libint2::Engine n_engine(libint2::Operator::nuclear, obs.max_nprim(), obs.max_l());
  n_engine.set_params(q);

  // auto shell2bf = obs.shell2bf();
  const auto& buf_vec = n_engine.results();
  const auto shell2bf = obs.shell2bf();
  for (int i = 0; i < obs.size(); i++) {
    int bf1 = shell2bf[i];
    int n1 = obs[i].size();
    for (int j = 0; j < obs.size(); j++) {
      int bf2 = shell2bf[j];
      int n2 = obs[j].size();
      n_engine.compute(obs[i], obs[j]);
      auto ints_shellset = buf_vec[0];
      if (ints_shellset == nullptr)
        continue;

      for (int k = 0; k < n1; k++) {
        for (int l = 0; l < n2; l++) {
          V_ne(bf1 + k, bf2 + l) += ints_shellset[k * n2 * l];
        }
      }
    }
  }
  return V_ne;
}

tensor4d electron_electron_repulsion(libint2::BasisSet obs) {
  int nbasis = obs.nbf();

  tensor4d V_ee(nbasis, nbasis, nbasis, nbasis);
  V_ee.setZero();


  libint2::Engine ee_engine(libint2::Operator::coulomb, obs.max_nprim(), obs.max_l());
  auto shell2bf = obs.shell2bf();
  const auto& buf_vec = ee_engine.results();

  for (int i = 0; i < obs.size(); i++) {
    const auto bf1 = shell2bf[i];
    const auto n1 = obs[i].size();
    for (int j = 0; j < obs.size(); j++) {
      const auto bf2 = shell2bf[j];
      const auto n2 = obs[j].size();
      for (int k = 0; k < obs.size(); k++) {
        const auto bf3 = shell2bf[k];
        const auto n3 = obs[k].size();
        for (int l = 0; l < obs.size(); l++) {
          const auto bf4 = shell2bf[l];
          const auto n4 = obs[l].size();

          ee_engine.compute(obs[i], obs[j], obs[k], obs[l]);
          auto ints_shellset = buf_vec[0];
          if (ints_shellset == nullptr)
            continue;

          int nprimitives_i = obs.size();
          int nprimitives_j = obs.size();
          int nprimitives_k = obs.size();
          int nprimitives_l = obs.size();

          for (int ii = 0; ii < n1; ii++) {
            for (int jj = 0; jj < n2; jj++) {
              for (int kk = 0; kk < n3; kk++) {
                for (int ll = 0; ll < n4; ll++) {

                  const size_t p = bf1 + ii;
                  const size_t q = bf2 + jj;
                  const size_t r = bf3 + kk;
                  const size_t s = bf4 + ll;

                  V_ee(p,q, r, s) += ints_shellset[ii * n2 * n3 * n4 + jj * n3 * n4 + kk * n4 + ll];
                }
              }
            }
          }
        }
      }
    }
  }
  return V_ee;
}

double
nuclear_nuclear_repulsion_energy(std::vector<libint2::Atom> atoms) {

  int natoms = atoms.size();

  double E_NN = 0.;

  for (int i = 0; i < natoms; i++) {
    for (int j = i + 1; j < natoms; j++) {
      auto Rijx = atoms[i].x - atoms[j].x;
      auto Rijy = atoms[i].y - atoms[j].y;
      auto Rijz = atoms[i].z - atoms[j].z;

      auto Rijx2 = Rijx * Rijx;
      auto Rijy2 = Rijy * Rijy;
      auto Rijz2 = Rijz * Rijz;
      auto Rij = std::sqrt(Rijx2 + Rijy2 + Rijz2);
      E_NN += (atoms[i].atomic_number * atoms[j].atomic_number) / Rij;
    }
  }
  return E_NN;
}

matrix2d compute_density_matrix(const matrix2d &mos, int n_occ) {
  int nbasis_functions = mos.rows();

  matrix2d density_matrix(nbasis_functions, nbasis_functions);
  density_matrix.setZero();

  double occupation = 2.0;

  for (int i = 0; i < nbasis_functions; i++) {
    for (int j = 0; j < nbasis_functions; j++) {
      for (int o = 0; o < n_occ; o++) {
        auto C = mos(i,o);
        auto C_dagger = mos(j,o);
        density_matrix(i,j) += occupation * C * C_dagger;
      }
    }
  }
  return density_matrix;
}

matrix2d compute_G(matrix2d dens_mat, tensor4d Vee) {
  int nbasis_functions = dens_mat.rows();
  matrix2d G(nbasis_functions, nbasis_functions);
  G.setZero();

  for (int i = 0; i < nbasis_functions; i++) {
    for (int j = 0; j < nbasis_functions; j++) {
      for (int k = 0; k < nbasis_functions; k++) {
        for (int l = 0; l < nbasis_functions; l++) {
          auto density = dens_mat(k,l);
          auto J = Vee(i,j,k,l);
          auto K = Vee(i,l,k,j);
          G(i,j) += density * (J - 0.5 * K);
        }
      }
    }
  }
  return G;
}

double compute_electronic_energy_expectation_value(matrix2d dens_mat,
                                                   matrix2d T, matrix2d Vne,
                                                   matrix2d G) {

  auto Hcore = T + Vne;
  auto electronic_energy = 0.0;
  auto nbasis_functions = dens_mat.rows();

  for (int i = 0; i < nbasis_functions; i++) {
    for (int j = 0; j < nbasis_functions; j++) {
      electronic_energy += dens_mat(i,j) * (Hcore(i,j) + 0.5 * G(i,j));
    }
  }
  return electronic_energy;
}

double
scf_cycle(std::tuple<matrix2d, matrix2d, matrix2d, tensor4d> molecular_terms,
          std::tuple<double, int> scf_parameters, libint2::BasisSet obs) {
  auto [S, T, Vne, Vee] = molecular_terms;
  auto [tolerance, max_iter] = scf_parameters;
  auto electronic_energy = 0.0;

  int nbasis_functions = obs.nbf();

  matrix2d dens_mat(nbasis_functions, nbasis_functions);
  dens_mat.setZero();

  for (int scf_step = 0; scf_step < max_iter; scf_step++) {
    auto electronic_energy_old = electronic_energy;
    auto G = compute_G(dens_mat, Vee);
    auto F = T + Vne + G;

    // auto S_eigen = S;
    // auto S_eigen_inv = S_eigen.inverse();
    // auto S_eigen_inv_sqrt = S_eigen_inv.sqrt();
    // auto S_inv_sqrt = S_eigen_inv_sqrt;

    Eigen::SelfAdjointEigenSolver<matrix2d> es(S);
    if (es.info() != Eigen::Success) {
      throw std::runtime_error("Man this stuff sucks at S");
    }

    Eigen::VectorXd evals = es.eigenvalues();
    Eigen::MatrixXd evecs = es.eigenvectors();

    // std::cout << "Eigenvalues of S:\n" << es.eigenvalues().transpose() << std::endl;

    double cond_number = es.eigenvalues().maxCoeff() / es.eigenvalues().minCoeff();
    // std::cout << "Condition number of S = " << cond_number << std::endl;

    double eps = 1e-12;

    for (int i = 0; i < evals.size(); i++) {
      if (evals[i] < eps) {
        // std::cerr << "Warning: overlap eigenvalue " << i << " = " << evals[i] << " < " << " (basis nearly linearly dependent)\n";
        evals[i] = 0.0;
      } else {
        evals[i] = 1.0 / std::sqrt(evals[i]);
      }
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

    dens_mat = compute_density_matrix(mos, 1);

    electronic_energy =
        compute_electronic_energy_expectation_value(dens_mat, T, Vne, G);

    // std::cout << "SCF Step:\t" << scf_step << "\tElectronic energy:\t" << electronic_energy << "\tdE:\t" << std::fabs(electronic_energy - electronic_energy_old) << std::endl;

    if (std::fabs(electronic_energy - electronic_energy_old) < tolerance) {
      break;
    }
  }
  return electronic_energy;
}
