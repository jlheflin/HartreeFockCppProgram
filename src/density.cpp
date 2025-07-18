#include <density.hpp>

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

matrix2d compute_G(const matrix2d &dens_mat, const tensor4d &Vee) {
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
