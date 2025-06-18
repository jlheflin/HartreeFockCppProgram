# Hartree Fock Python Program Notes

## Types/Containers

- primitive_gaussian = (alpha, coefficient, coordinates, l1, l2, l3)
  - class
  - The l1, l2, l3 variables are not used in this case
  - in the code the main variables are H1_pg1a, H1_pg1b, H1_pg1c
  - Instance A = ( 2.0 * alpha / math.pi ) ^ 0.75

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

## Functions

### overlap(molecule)
- Input: molecule
- Output: an N x N matrix, where N is the number of atomic orbtials
  - N is determined by nbasis, which is allocated via len(molecule)

```python
for i in range(nbasis):
  for j in range(nbasis):
    # Accesses the atomic orbital, returning the number of primitive gaussians associated with orbital
    nprimitives_i = len(molecule[i])
    nprimitives_j = len(molecule[j])

    for k in range(nprimitives_i):
      for l in range(nprimitives_j):
      # Loops over all of the primitives
        N = molecule[i][k].A * molecule[j][l].A
        p = molecule[i][k].alpha + molecule[j][l].alpha
        q = molecule[i][k].alpha * molecule[j][l].alpha / p
        Q = molecule[i][k].coordinates - molecule[j][l].coordinates
        Q2 = np.dot(Q,Q)

        S[i,j] += N * molecule[i][k].coeff * molecule[j][l].coeff * math.exp(-q*Q2) * (math.pi/p) ^ (3/2)
```

### kinetic(molecule)
- Input: molecule
- Output: an N x N matrix

```python
```
