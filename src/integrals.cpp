#include <integrals.hpp>

using tensor4d = Eigen::Tensor<double, 4>;
using matrix2d = Eigen::MatrixXd;


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
          S(bf1 + k, bf2 + l) += ints_shellset[k * n2 + l];
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
          T(bf1 + k, bf2 + l) += ints_shellset[k * n2 + l];
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
          V_ne(bf1 + k, bf2 + l) += ints_shellset[k * n2 + l];
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
