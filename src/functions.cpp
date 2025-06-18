#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include <classes.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <format>

using tensor4d = std::vector<std::vector<std::vector<std::vector<double>>>>;
using matrix2d = std::vector<std::vector<double>>;

double dot(const coord_type &a, const coord_type &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

matrix2d add_matrices(const matrix2d &A, const matrix2d &B) {
  assert(A.size() == B.size());
  assert(A[0].size() == B[0].size());

  int rows = A.size();
  int cols = A[0].size();

  matrix2d result(rows, std::vector<double>(cols, 0.0));

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      result[i][j] = A[i][j] + B[i][j];
    }
  }
  return result;
}

matrix2d mult_matrices(const matrix2d &A, const matrix2d &B) {
  int rows = A.size();
  int cols = A[0].size();
  int inner = B.size();

  assert(A[0].size() == inner);

  matrix2d result(rows, std::vector<double>(cols, 0.0));

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      for (int k = 0; k < inner; k++) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return result;
}

Eigen::MatrixXd to_eigen(const std::vector<std::vector<double>> &mat) {
  int rows = mat.size();
  int cols = mat[0].size();
  Eigen::MatrixXd eigen_mat(rows, cols);
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j)
      eigen_mat(i, j) = mat[i][j];
  return eigen_mat;
}

std::vector<std::vector<double>> from_eigen(const Eigen::MatrixXd &mat) {
  std::vector<std::vector<double>> result(mat.rows(),
                                          std::vector<double>(mat.cols()));
  for (int i = 0; i < mat.rows(); ++i)
    for (int j = 0; j < mat.cols(); ++j)
      result[i][j] = mat(i, j);
  return result;
}

std::vector<std::vector<double>> overlap(molecule mol) {
  int nbasis = mol.size();

  std::vector<std::vector<double>> S(nbasis, std::vector<double>(nbasis, 0.0));

  for (int i = 0; i < nbasis; i++) {
    for (int j = 0; j < nbasis; j++) {
      int nprimitives_i = mol[i].size();
      int nprimitives_j = mol[j].size();

      for (int k = 0; k < nprimitives_i; k++) {
        for (int l = 0; l < nprimitives_j; l++) {
          auto N = mol[i][k].A * mol[j][l].A;
          auto p = mol[i][k].alpha + mol[j][l].alpha;
          auto q = (mol[i][k].alpha * mol[j][l].alpha) / p;
          auto Q = mol[i].coords - mol[j].coords;
          auto Q2 = dot(Q, Q);
          S[i][j] += N * mol[i][k].coeff * mol[j][l].coeff * std::exp(-q * Q2) *
                     pow((M_PI / p), (3. / 2.));
        }
      }
    }
  }
  return S;
}

std::vector<std::vector<double>> kinetic(molecule mol) {
  int nbasis = mol.size();

  std::vector<std::vector<double>> T(nbasis, std::vector<double>(nbasis, 0.0));

  for (int i = 0; i < nbasis; i++) {
    for (int j = 0; j < nbasis; j++) {
      int nprimitives_i = mol[i].size();
      int nprimitives_j = mol[j].size();

      for (int k = 0; k < nprimitives_i; k++) {
        for (int l = 0; l < nprimitives_j; l++) {
          auto N = mol[i][k].A * mol[j][l].A;
          auto cacb = mol[i][k].coeff * mol[j][l].coeff;

          auto p = mol[i][k].alpha + mol[j][l].alpha;
          auto P =
              mol[i][k].alpha * mol[i].coords + mol[j][l].alpha * mol[j].coords;
          auto Pp = P / p;
          auto PG = Pp - mol[j].coords;
          auto PGx2 = PG.x * PG.x;
          auto PGy2 = PG.y * PG.y;
          auto PGz2 = PG.z * PG.z;

          auto q = (mol[i][k].alpha * mol[j][l].alpha) / p;
          auto Q = mol[i].coords - mol[j].coords;
          auto Q2 = dot(Q, Q);

          auto s = std::exp(-q * Q2) * pow((M_PI / p), 1.5) * N * cacb;

          T[i][j] += 3.0 * mol[j][l].alpha * s;
          T[i][j] -=
              2.0 * mol[j][l].alpha * mol[j][l].alpha * s * (PGx2 + 0.5 / p);
          T[i][j] -=
              2.0 * mol[j][l].alpha * mol[j][l].alpha * s * (PGy2 + 0.5 / p);
          T[i][j] -=
              2.0 * mol[j][l].alpha * mol[j][l].alpha * s * (PGz2 + 0.5 / p);
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

std::vector<std::vector<double>>
electron_nuclear_attraction(molecule mol, std::vector<int> Z_list) {
  int natoms = Z_list.size();
  int nbasis = mol.size();

  std::vector<coord_type> coordinates;
  for (int i = 0; i < nbasis; i++) {
    coordinates.push_back(mol[i].coords);
  }

  std::sort(coordinates.begin(), coordinates.end(),
            [](const coord_type &a, const coord_type &b) {
              if (a.x != b.x)
                return a.x < b.x;
              if (a.y != b.y)
                return a.y < b.y;
              return a.z < b.z;
            });

  coordinates.erase(std::unique(coordinates.begin(), coordinates.end(),
                                [](const coord_type &a, const coord_type &b) {
                                  return a.x == b.x && a.y == b.y && a.z == b.z;
                                }),
                    coordinates.end());

  std::vector<std::vector<double>> V_ne(nbasis,
                                        std::vector<double>(nbasis, 0.0));

  for (int atom = 0; atom < natoms; atom++) {
    for (int i = 0; i < nbasis; i++) {
      for (int j = 0; j < nbasis; j++) {
        int nprimitives_i = mol[i].size();
        int nprimitives_j = mol[j].size();

        for (int k = 0; k < nprimitives_i; k++) {
          for (int l = 0; l < nprimitives_j; l++) {

            auto N = mol[i][k].A * mol[j][l].A;
            auto cacb = mol[i][k].coeff * mol[j][l].coeff;

            auto p = mol[i][k].alpha + mol[j][l].alpha;
            auto P = mol[i][k].alpha * mol[i].coords +
                     mol[j][l].alpha * mol[j].coords;
            auto Pp = P / p;

            auto PG = P / p - coordinates[atom];
            auto PG2 = dot(PG, PG);

            auto q = mol[i][k].alpha * mol[j][l].alpha / p;
            auto Q = mol[i].coords - mol[j].coords;
            auto Q2 = dot(Q, Q);

            V_ne[i][j] += N * cacb * -Z_list[atom] * (2.0 * M_PI / p) *
                          std::exp(-q * Q2) * boys(p * PG2, 0);
          }
        }
      }
    }
  }
  return V_ne;
}

tensor4d electron_electron_repulsion(molecule mol) {
  int nbasis = mol.size();

  tensor4d V_ee(nbasis,
                std::vector<std::vector<std::vector<double>>>(
                    nbasis, std::vector<std::vector<double>>(
                                nbasis, std::vector<double>(nbasis, 0.0))));

  for (int i = 0; i < nbasis; i++) {
    for (int j = 0; j < nbasis; j++) {
      for (int k = 0; k < nbasis; k++) {
        for (int l = 0; l < nbasis; l++) {
          int nprimitives_i = mol[i].size();
          int nprimitives_j = mol[j].size();
          int nprimitives_k = mol[k].size();
          int nprimitives_l = mol[l].size();

          for (int ii = 0; ii < nprimitives_i; ii++) {
            for (int jj = 0; jj < nprimitives_j; jj++) {
              for (int kk = 0; kk < nprimitives_k; kk++) {
                for (int ll = 0; ll < nprimitives_l; ll++) {

                  auto N =
                      mol[i][ii].A * mol[j][jj].A * mol[k][kk].A * mol[l][ll].A;
                  auto cicjckcl = mol[i][ii].coeff * mol[j][jj].coeff *
                                  mol[k][kk].coeff * mol[l][ll].coeff;

                  auto pij = mol[i][ii].alpha + mol[j][jj].alpha;
                  auto pkl = mol[k][kk].alpha + mol[l][ll].alpha;

                  auto Pij = mol[i][ii].alpha * mol[i].coords +
                             mol[j][jj].alpha * mol[j].coords;
                  auto Pkl = mol[k][kk].alpha * mol[k].coords +
                             mol[l][ll].alpha * mol[l].coords;

                  auto Ppij = Pij / pij;
                  auto Ppkl = Pkl / pkl;

                  auto PpijPpkl = Ppij - Ppkl;
                  auto PpijPpkl2 = dot(PpijPpkl, PpijPpkl);
                  auto denom = 1.0 / pij + 1.0 / pkl;

                  auto qij = mol[i][ii].alpha * mol[j][jj].alpha / pij;
                  auto qkl = mol[k][kk].alpha * mol[l][ll].alpha / pkl;

                  auto Qij = mol[i].coords - mol[j].coords;
                  auto Qkl = mol[k].coords - mol[l].coords;

                  auto Q2ij = dot(Qij, Qij);
                  auto Q2kl = dot(Qkl, Qkl);

                  auto term1 = 2.0 * M_PI * M_PI / (pij * pkl);
                  auto term2 = std::sqrt(M_PI / (pij + pkl));
                  auto term3 = std::exp(-qij * Q2ij);
                  auto term4 = std::exp(-qkl * Q2kl);

                  V_ee[i][j][k][l] += N * cicjckcl * term1 * term2 * term3 *
                                      term4 * boys(PpijPpkl2 / denom, 0);
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

std::vector<std::vector<double>>
compute_density_matrix(const std::vector<std::vector<double>> &mos, int n_occ) {
  int nbasis_functions = mos.size();

  std::vector<std::vector<double>> density_matrix(
      nbasis_functions, std::vector<double>(nbasis_functions, 0.0));

  double occupation = 2.0;

  for (int i = 0; i < nbasis_functions; i++) {
    for (int j = 0; j < nbasis_functions; j++) {
      for (int o = 0; o < n_occ; o++) {
        auto C = mos[i][o];
        auto C_dagger = mos[j][o];
        density_matrix[i][j] += occupation * C * C_dagger;
      }
    }
  }
  return density_matrix;
}

std::vector<std::vector<double>>
compute_G(std::vector<std::vector<double>> dens_mat, tensor4d Vee) {
  int nbasis_functions = dens_mat.size();
  std::vector<std::vector<double>> G(
      nbasis_functions, std::vector<double>(nbasis_functions, 0.0));

  for (int i = 0; i < nbasis_functions; i++) {
    for (int j = 0; j < nbasis_functions; j++) {
      for (int k = 0; k < nbasis_functions; k++) {
        for (int l = 0; l < nbasis_functions; l++) {
          auto density = dens_mat[k][l];
          auto J = Vee[i][j][k][l];
          auto K = Vee[i][l][k][j];
          G[i][j] += density * (J - 0.5 * K);
        }
      }
    }
  }
  return G;
}

double compute_electronic_energy_expectation_value(matrix2d dens_mat,
                                                   matrix2d T, matrix2d Vne,
                                                   matrix2d G) {

  auto Hcore = add_matrices(T, Vne);
  auto electronic_energy = 0.0;
  auto nbasis_functions = dens_mat.size();

  for (int i = 0; i < nbasis_functions; i++) {
    for (int j = 0; j < nbasis_functions; j++) {
      electronic_energy += dens_mat[i][j] * (Hcore[i][j] + 0.5 * G[i][j]);
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

  matrix2d dens_mat(nbasis_functions,
                    std::vector<double>(nbasis_functions, 0.0));

  for (int scf_step = 0; scf_step < max_iter; scf_step++) {
    auto electronic_energy_old = electronic_energy;

    auto G = compute_G(dens_mat, Vee);

    auto F = add_matrices(add_matrices(T, Vne), G);

    auto S_eigen = to_eigen(S);
    auto S_eigen_inv = S_eigen.inverse();
    auto S_eigen_inv_sqrt = S_eigen_inv.sqrt();
    auto S_inv_sqrt = from_eigen(S_eigen_inv_sqrt);

    auto F_unitS = mult_matrices(S_inv_sqrt, mult_matrices(F, S_inv_sqrt));
    auto F_unitS_eigen = to_eigen(F_unitS);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F_unitS_eigen);

    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Eigenvalue decomp failed");
    }

    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    auto mos = mult_matrices(S_inv_sqrt, from_eigen(eigenvectors));
    auto density_matrix = compute_density_matrix(mos, 1);
    electronic_energy =
        compute_electronic_energy_expectation_value(density_matrix, T, Vne, G);

    std::cout << "Iteration: " << scf_step
              << " , Electronic energy: " << electronic_energy << std::endl;

    if (std::fabs(electronic_energy - electronic_energy_old) < tolerance) {
      break;
    }
  }
  return electronic_energy;
}
