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
  auto nbasis = obs.size();
  matrix2d S;
  S.setZero(nbasis, nbasis);

  libint2::Engine s_engine(libint2::Operator::overlap, obs.max_nprim(), obs.max_l());
  // auto shell2bf = obs.shell2bf();
  const auto& buf_vec = s_engine.results();
  for (int i = 0; i < obs.size(); i++) {
    for (int j = 0; j < obs.size(); j++) {
      s_engine.compute(obs[i], obs[j]);
      auto ints_shellset = buf_vec[0];
      if (ints_shellset == nullptr)
        continue;

      auto n1 = obs[i].size();
      auto n2 = obs[j].size();

      for (int k = 0; k < n1; k++) {
        for (int l = 0; l < n2; l++) {
          S(i, j) += ints_shellset[k * n2 * l];
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
  for (int i = 0; i < obs.size(); i++) {
    for (int j = 0; j < obs.size(); j++) {
      t_engine.compute(obs[i], obs[j]);
      auto ints_shellset = buf_vec[0];
      if (ints_shellset == nullptr)
        continue;

      auto n1 = obs[i].size();
      auto n2 = obs[j].size();

      for (int k = 0; k < n1; k++) {
        for (int l = 0; l < n2; l++) {
          T(i, j) += ints_shellset[k * n2 * l];
        }
      }
    }
  }
  return T;
}

double boys(double x, int n) {
  if (x == 0) {
    return 1.0 / (2 * n + 1);
  } else {
    return boost::math::gamma_p(n + 0.5, x) * boost::math::tgamma(n + 0.5) *
           (1.0 / (2 * pow(x, (n + 0.5))));
  }
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
  for (int i = 0; i < obs.size(); i++) {
    for (int j = 0; j < obs.size(); j++) {
      n_engine.compute(obs[i], obs[j]);
      auto ints_shellset = buf_vec[0];
      if (ints_shellset == nullptr)
        continue;

      auto n1 = obs[i].size();
      auto n2 = obs[j].size();

      for (int k = 0; k < n1; k++) {
        for (int l = 0; l < n2; l++) {
          V_ne(i, j) += ints_shellset[k * n2 * l];
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

                  const size_t i = bf1 + ii;
                  const size_t j = bf2 + jj;
                  const size_t k = bf3 + kk;
                  const size_t l = bf4 + ll;

                  V_ee(i,j, k, l) = ints_shellset[((ii * n2 * jj) * n3 * kk) * n4 + ll];
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
nuclear_nuclear_repulsion_energy(const std::vector<coord_type> &coord_list,
                                 const std::vector<int> &Z_list) {
  assert(coord_list.size() == Z_list.size());

  int natoms = Z_list.size();

  double E_NN = 0.;

  for (int i = 0; i < natoms; i++) {
    auto Zi = Z_list[i];
    for (int j = 0; j < natoms; j++) {
      if (j > i) {
        auto Zj = Z_list[j];

        auto Rijx = coord_list[i].x - coord_list[j].x;
        auto Rijy = coord_list[i].y - coord_list[j].y;
        auto Rijz = coord_list[i].z - coord_list[j].z;

        auto Rijx2 = Rijx * Rijx;
        auto Rijy2 = Rijy * Rijy;
        auto Rijz2 = Rijz * Rijz;
        auto Rij = std::sqrt(Rijx2 + Rijy2 + Rijz2);
        E_NN += (Zi * Zj) / Rij;
      }
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
          std::tuple<double, int> scf_parameters, molecule mol) {
  auto [S, T, Vne, Vee] = molecular_terms;
  auto [tolerance, max_iter] = scf_parameters;
  auto electronic_energy = 0.0;

  int nbasis_functions = mol.size();

  matrix2d dens_mat(nbasis_functions, nbasis_functions);
  dens_mat.setZero();

  for (int scf_step = 0; scf_step < max_iter; scf_step++) {
    auto electronic_energy_old = electronic_energy;
    auto G = compute_G(dens_mat, Vee);
    auto F = T + Vne + G;

    auto S_eigen = S;
    auto S_eigen_inv = S_eigen.inverse();
    // std::cout << "S_eigen_inv: \n" << S_eigen_inv << std::endl;
    auto S_eigen_inv_sqrt = S_eigen_inv.sqrt();
    auto S_inv_sqrt = S_eigen_inv_sqrt;

    auto F_unitS = S_inv_sqrt * F * S_inv_sqrt;
    auto F_unitS_eigen = F_unitS;

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

    if (std::fabs(electronic_energy - electronic_energy_old) < tolerance) {
      break;
    }
  }
  return electronic_energy;
}
