# Hartree Fock Python Program Notes

## Types/Containers

- primitive_gaussian = (alpha, coefficient, coordinates, l1, l2, l3)
  - class
  - The l1, l2, l3 variables are not used in this case
  - in the code the main variables are H1_pg1a, H1_pg1b, H1_pg1c

- zlist = \[Z1, Z2\]
  - Just a list of the atomic numbers of the atoms involved in the system

- number_occupied_orbitals = \#
  - Number of occupied orbitals

- atom_coods = \[\[ x, y, z \], \[ x, y, z \]\]
  - An array of the atom coordinates of the system

- atomic_orbital = \[primitive_gaussian, ... \]
  - A list of primitive gaussians that make up the atomic orbital of the atom
  - In the python code it is H1_1s and H2_1s, with the list including H1_pg1a, H1_pg1b, and H1_pg1c

- molecule = \[atomic_orbital1, atomic_orbital2, ...\]
  - Includes the atomic orbitals of the system
  - Breakdown:
    - molecule contains
      - atomic orbitals which are lists of
        - primitive gaussians, which have
          - alphas, coeffs, coords, and l1, l2, l3

## Containers

f
