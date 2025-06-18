#include <classes.hpp>

std::vector<std::vector<double>> overlap(orbital_set aobs) {
  int nbasis = aobs.size();

  std::vector<std::vector<double>> S(nbasis, std::vector<double>(nbasis, 0.0));

  for (int i = 0; i < nbasis; i++) {
    for (int j = 0; j < nbasis; j++) {
      int nprimitives_i = aobs.size();
      int nprimitives_j = aobs.size();

      std::cout << nprimitives_i << std::endl;
    }
  }

  return S;
  // for (int i = 0; i < nbasis; i++) {
  //   for (int j = 0; j < nbasis; j++) {

  //   }
  // }
}
