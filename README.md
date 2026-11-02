# IntraHB

**Intramolecular hydrogen-bond energies and cooperativity by the Molecular Tailoring Approach (MTA).**

Developed by **Deepak Patkar**, based on the method of Deshmukh, Bartolotti & Gadre,
*J. Phys. Chem. A* 2008, 112, 312–321.

IntraHB finds the intramolecular hydrogen bonds in a molecule, cuts out the donor and acceptor groups,
writes the Gaussian inputs for every capped fragment, and turns the finished calculations into
hydrogen-bond energies and cooperativity contributions.

- Automatic H-bond detection (distance, angle and ring-size criteria), or manual definition
- Automatic donor / acceptor groups (‑OH, ‑CHO, ‑COOH, ring systems, …), overridable per bond
- Detection of H-bond networks and cooperativity energies
- One Gaussian job per distinct structure, with no duplicates
- Pure Python 3.9+, standard library only

---

## Method

For a hydrogen bond D‑H···A in a molecule **M**, the donor group **Fd** and the acceptor group **Fa**
are removed. Each broken bond is capped with a hydrogen atom placed along the original bond, and
nothing is re-optimised:

$$E_{HB} = E(M\text{-}F_d) + E(M\text{-}F_a) - E(M\text{-}F_d\text{-}F_a) - E(M)$$

If the bond is part of a network, the same expression is evaluated in **M′**, the molecule with the
groups of the neighbouring H-bonds removed. The difference is the cooperativity contribution:

$$E_{coop} = E_{HB}(M) - E_{HB}(M')$$

> Deshmukh, M. M.; Bartolotti, L. J.; Gadre, S. R. *J. Phys. Chem. A* **2008**, 112, 312–321.
> [doi:10.1021/jp076316b](https://doi.org/10.1021/jp076316b)

### Fragment naming

| Name        | Meaning                                  |
|-------------|------------------------------------------|
| `M`         | parent molecule                          |
| `F1, F2 …`  | groups cut out of M (donor / acceptor)   |
| `M-F1`      | M with F1 removed and capped             |
| `M-F1-F3`   | M with F1 and F3 removed and capped      |

Gaussian inputs are named `<molecule>_<structure>.gjf`, for example `rahb_M-F1-F2.gjf`.

---

## Installation

```bash
git clone https://github.com/<user>/IntraHB.git
cd IntraHB
python intrahb.py --help
```

No packages to install. Running jobs needs Gaussian (`g16`) on the `PATH`.

---

## Quick start

```bash
python intrahb.py INPUT.cfg --mode detect   # list H-bonds and fragments
python intrahb.py INPUT.cfg --mode gen      # write the Gaussian inputs
python intrahb.py INPUT.cfg --mode run      # write inputs and run g16
python intrahb.py INPUT.cfg --mode cal      # read logs, compute energies
```

| Mode     | What it does                                                        |
|----------|---------------------------------------------------------------------|
| `detect` | prints H-bonds, fragments and the job list; writes nothing          |
| `gen`    | writes `<molecule>/<molecule>_*.gjf` and a summary                  |
| `run`    | as `gen`, then runs every job that has no finished log (`njob` in parallel) |
| `cal`    | reads the `.log` files and writes `<molecule>/<molecule>_summary.txt` |

### Using existing calculations

If the fragment calculations were already done elsewhere, under any file names, point IntraHB to that
folder:

```bash
python intrahb.py INPUT.cfg --mode cal --ref path/to/old_calculation
```

Each required structure is matched to a reference input by geometry. Matching ignores atom order,
rotation and reflection, and the energy is read from the matching `.log` file.

---

## Input

Supported structure formats: `.xyz` and Gaussian inputs (`.gjf`, `.com`, `.inp`).

### `INPUT.cfg`

```ini
[SYSTEM]
xyz_file = [rahb.xyz]           # one or more files, wildcards allowed
nproc    = 8
mem      = 12GB
outdir   = .                    # one sub-folder per molecule is created here

[THEORY]
method       = MP2
basis        = aug-cc-pVTZ
keywords     = scf=tight
charge       = auto             # auto: from a Gaussian input, else 0
multiplicity = auto             # auto: from a Gaussian input, else from electron count

[HB]
hb_donor        = N,O,F
hb_acceptor     = N,O,F,S,Cl
hb_distance_min = 1.2           # H···A (Å)
hb_distance_max = 2.5
hb_angle_min    = 100           # D-H···A (degrees)
hb_angle_max    = 180
hb_min_bonds    = 3             # min. covalent bonds between D and A

[FRAGMENT]
cap_distance       = 1.0        # cap H distance (Å)
network_link_bonds = 2          # HBs whose D/A atoms are this close (in bonds) form a network

[RUN]
mode        = gen               # detect | gen | run | cal
njob        = 1
run_command = g16 {input}       # executed inside the molecule folder
```

Section names are only for readability. Text after ` #` is a comment.

### Manual H-bonds

Any `hbN` line switches off automatic detection. Atom numbers are 1-based:

```ini
hb1 = 3 4 10                    # D H A
hb2 = 11 12 5  agroup=1,3-5     # optional dgroup= / agroup= override the automatic groups
```

### How groups are chosen

Walking along the covalent path between donor and acceptor, a group is cut at the first bond that
does not lie in a ring:

- **Donor:** D‑H plus whatever hangs on D (‑OH → OH, ring N‑H → the ring).
- **Acceptor:** a terminal atom without H (C=O, N=O) brings its carrier atom along (‑CHO, ‑COOH).
  Any other acceptor is cut at the first such bond (OH, OMe, ring).

If you need a different cut, give it explicitly with `dgroup=` / `agroup=`.

---

## Example: 2-RAHB (C₄H₄O₄), MP2/aug-cc-pVTZ

```
IntraHB 1.0  -  intramolecular H-bond energies by MTA        developed by Deepak Patkar
based on: Deshmukh, Bartolotti, Gadre, J. Phys. Chem. A 2008, 112, 312-321
...

FRAGMENTS  (cut out of M, open valences capped with H at 1.00 A)
  F1   O3 H4                     HB1 donor
  F2   C8 H9 O10                 HB1 acceptor
  F3   O11 H12                   HB2 donor
  F4   C1 O3 H4 O5               HB2 acceptor

HYDROGEN BONDS
  HB    D-H...A             r(D-H)  r(H..A)  r(D..A)    angle      E_HB  E_HB(M')   E_coop   network
  HB1   O3-H4...O10          0.995    1.734    2.652   151.55     9.420     6.560    2.860   HB2
  HB2   O11-H12...O5         1.012    1.631    2.556   149.68    17.607    15.297    2.310   HB1

HB1  O3-H4...O10   donor F1, acceptor F2,   M' = M-F3
              in M                            in M' (network broken)
  E(M)        M                  -455.02705822    M-F3                 -379.85905501
  E(M-Fd)     M-F1               -379.84856394    M-F1-F3              -304.68477772
  E(M-Fa)     M-F2               -341.85377151    M-F2-F3              -266.69427289
  E(M-Fd-Fa)  M-F1-F2            -266.69028905    M-F1-F2-F3           -191.53044962
  E_HB                                   9.420                                 6.560    E_coop = 2.860
```

All energies are in kcal/mol. The summary file also lists every job with its status and energy,
followed by the coordinates of all capped structures.

---

## Energies read from Gaussian logs

The highest level available in **every** log is used, in this order:
`CCSD(T)` → `MP4(SDTQ)` → `MP3` → `MP2` → `SCF` (HF / DFT).

## Notes and limitations

- Capped geometries are not re-optimised, as the method requires.
- Bifurcated H-bonds (one H shared by two acceptors) are reported in the summary notes, because
  their energies are not independent.
- A note is printed if a capped fragment's electron count does not fit the chosen multiplicity.

## Author

**Deepak Patkar** wrote IntraHB as an implementation of the Molecular Tailoring Approach for
intramolecular hydrogen bonds, following the reference paper given under [Method](#method).
