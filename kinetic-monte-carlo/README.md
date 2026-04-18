# Kinetic Monte Carlo — Fe Vacancy Hopping in Maghemite

Fortran implementation of a kinetic Monte Carlo (KMC) simulation modeling
cation vacancy hopping in a maghemite (γ-Fe₂O₃) superblock. Developed to
investigate cation disorder and vacancy migration pathways relevant to
surface reactivity and magnetic properties.

## Physical Background

Maghemite contains two distinct Fe sublattice sites:
- **Tetrahedral (Td)** sites — labeled `Mn` in the input XYZ file
- **Octahedral (Oh)** sites — labeled `Fe`
- **V1 vacancies** — tetrahedral vacancies, labeled `B`
- **V2 vacancies** — octahedral vacancies, labeled `C`

Vacancy hopping probability is governed by the Arrhenius equation:

    P_hop = exp(-Ea / kB·T)

where Ea is the migration activation energy, kB is the Boltzmann constant
(8.617e-5 eV/K), and T is temperature (default: 300 K).

## Files

| File | Description |
|------|-------------|
| `hopping_v1.f90` | Original KMC — single activation energy, basic random flip |
| `hopping_v2_multipath.f90` | Improved — 10 migration pathways, weighted path selection, statistical output, dynamic memory allocation |

## Key Improvements: v1 → v2

- **Multi-pathway kinetics**: 10 physically distinct migration pathways
  (4 for V1, 6 for V2) with DFT-derived activation energies (0.76–0.83 eV)
- **Weighted path selection**: cumulative probability weighting replaces
  uniform random selection
- **Dynamic allocation**: allocatable arrays replace static fixed-size arrays
- **Statistical reporting**: per-pathway migration counts, relative fractions,
  and accumulated probabilities
- **Robust I/O**: iostat error handling and implicit none for type safety

## Input

- `100_100superblock.xyz` — XYZ coordinate file (~680,000 atoms)
- Atom labels: `Mn` (Td-Fe), `Fe` (Oh-Fe), `B` (V1 vacancy), `C` (V2 vacancy)

## Output

| File | Contents |
|------|----------|
| `outfile1.data` | Full atom list from input |
| `v1neigh.xyz` | V1 vacancies with Oh-Fe neighbors |
| `v2neigh.xyz` | V2 vacancies with Oh-Fe neighbors |
| `switched.xyz` | Final structure after vacancy hopping |
| stdout | Pathway statistics table |

## Compilation

```bash
# gfortran
gfortran -O2 -o hopping hopping_v2_multipath.f90

# Intel Fortran
ifort -O2 -o hopping hopping_v2_multipath.f90
```

## Usage

```bash
./hopping
```

Ensure `100_100superblock.xyz` is in the working directory before running.

## Parameters to Adjust

```fortran
temperature = 300.0   ! Temperature in Kelvin
numattempts = int(nv1v2 * 10)  ! Number of MC steps

real(kindr), dimension(10) :: activation_energies = (/ &
    0.79, 0.79, 0.76, 0.83, &        ! V1 pathways (I, IV, II, III)
    0.76, 0.76, 0.76, 0.79, 0.79, 0.79 /)  ! V2 pathways (V-X)
```

Update activation energies from your DFT-NEB calculations.

## Related Publication

- Salawu O.A., et al. "Physisorption and Ortho-Para Conversion of H₂
  on γ-Fe₂O₃ (001)." J. Phys. Chem. C, 2025, 129, 12679.
