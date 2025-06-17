#include <classes.hpp>

std::vector<std::vector<double>> overlap(orbital_set aobs) {
  int nbasis = aobs.size();

  for (auto gaussian_set : aobs)

    std::vector<std::vector<double>> S(nbasis,
                                       std::vector<double>(nbasis, 0.0));

  return S;
  // for (int i = 0; i < nbasis; i++) {
  //   for (int j = 0; j < nbasis; j++) {

  //   }
  // }
}
