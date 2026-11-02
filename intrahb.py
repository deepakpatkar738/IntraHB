#!/usr/bin/env python3
"""
IntraHB - intramolecular hydrogen-bond energies by the Molecular Tailoring Approach (MTA).

Author: Deepak Patkar
Based on: Deshmukh, Bartolotti & Gadre, J. Phys. Chem. A 2008, 112, 312-321.

For a hydrogen bond D-H...A inside a molecule M, the donor group Fd and the acceptor
group Fa are cut out, open valences are capped with H, nothing is re-optimised:

    E_HB   = E(M-Fd) + E(M-Fa) - E(M-Fd-Fa) - E(M)

If the bond is part of a network, the same expression is evaluated in M', the molecule
with the groups of the neighbouring bonds removed, and the difference is cooperativity:

    E_coop = E_HB(M) - E_HB(M')

Usage
    python intrahb.py [INPUT.cfg]                     mode from the cfg
    python intrahb.py INPUT.cfg --mode detect         detect | gen | run | cal
    python intrahb.py INPUT.cfg --mode cal --ref DIR  energies from an existing calculation

Only the Python standard library is required.
"""
from __future__ import annotations

import argparse
import configparser
import math
import re
import subprocess
import sys
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime
from itertools import combinations
from pathlib import Path

__author__ = "Deepak Patkar"
__version__ = "1.0"

HARTREE_TO_KCAL = 627.5094740631
BOND_TOLERANCE = 1.2  # bonded if distance <= 1.2 x (sum of covalent radii)
MODES = ("detect", "gen", "run", "cal")

ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
            "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I").split()
ATOMIC_NUMBER = {sym: z for z, sym in enumerate(ELEMENTS, 1)}
COVALENT_RADIUS = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
                   "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Se": 1.20, "Br": 1.20, "I": 1.39}
MAX_NEIGHBOURS = {"N": 3, "O": 2, "S": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1}  # beyond: no lone pair

_NUMBER = r"(-?\d+\.\d+(?:[DEde][+-]?\d+)?)"
ENERGY_PATTERNS = [  # highest level first; the best level present in every log is used
    ("CCSD(T)", re.compile(r"CCSD\(T\)\s*=\s*" + _NUMBER)),
    ("MP4", re.compile(r"MP4\(SDTQ\)\s*=\s*" + _NUMBER)),
    ("MP3", re.compile(r"EUMP3\s*=\s*" + _NUMBER)),
    ("MP2", re.compile(r"EUMP2\s*=\s*" + _NUMBER)),
    ("SCF", re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*" + _NUMBER)),
]


class InputError(Exception):
    """A problem with the user's input that stops the analysis of one molecule."""


# ----------------------------------------------------------------------------- settings
@dataclass
class Settings:
    structures: list[Path]
    outdir: Path
    mode: str = "gen"
    nproc: int = 8
    mem: str = "8GB"
    method: str = "MP2"
    basis: str = "aug-cc-pVTZ"
    keywords: str = ""
    charge: int | None = None
    mult: int | None = None
    command: str = "g16 {input}"
    njob: int = 1
    donors: set[str] = field(default_factory=lambda: {"N", "O", "F"})
    acceptors: set[str] = field(default_factory=lambda: {"N", "O", "F", "S", "Cl"})
    hb_distance: tuple[float, float] = (1.2, 2.5)
    hb_angle: tuple[float, float] = (100.0, 180.0)
    min_bonds: int = 3
    link_bonds: int = 2
    cap_distance: float = 1.0
    manual: dict[int, str] = field(default_factory=dict)

    @property
    def route(self) -> str:
        return f"#p {self.method}/{self.basis} {self.keywords}".strip()


def load_settings(path: Path) -> Settings:
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",), interpolation=None, strict=False)
    parser.read_string("[_]\n" + path.read_text(encoding="utf-8-sig", errors="replace"))
    raw = {key: value.strip() for name in parser.sections() for key, value in parser[name].items()}
    get = raw.get

    def optional_int(key):
        value = get(key, "auto")
        return None if value.lower() == "auto" else int(value)

    def elements(key, default):
        return {e.capitalize() for e in re.split(r"[,\s]+", get(key, default)) if e}

    base = path.parent
    structures = []
    for token in re.findall(r"[^\s,\[\]]+", get("xyz_file", "")):
        structures += sorted(base.glob(token)) if any(c in token for c in "*?") else [base / token]

    mode = get("mode", "gen").lower()
    if mode not in MODES:
        raise InputError(f"mode must be one of {', '.join(MODES)}")
    return Settings(
        structures=structures,
        outdir=base / get("outdir", "."),
        mode=mode,
        nproc=int(get("nproc", 8)),
        mem=get("mem", "8GB"),
        method=get("method", "MP2"),
        basis=get("basis", "aug-cc-pVTZ"),
        keywords=get("keywords", ""),
        charge=optional_int("charge"),
        mult=optional_int("multiplicity"),
        command=get("run_command", "g16 {input}"),
        njob=max(1, int(get("njob", 1))),
        donors=elements("hb_donor", "N,O,F"),
        acceptors=elements("hb_acceptor", "N,O,F,S,Cl"),
        hb_distance=(float(get("hb_distance_min", 1.2)), float(get("hb_distance_max", 2.5))),
        hb_angle=(float(get("hb_angle_min", 100)), float(get("hb_angle_max", 180))),
        min_bonds=int(get("hb_min_bonds", 3)),
        link_bonds=int(get("network_link_bonds", 2)),
        cap_distance=float(get("cap_distance", 1.0)),
        manual={int(k[2:]): v for k, v in raw.items() if re.fullmatch(r"hb\d+", k)},
    )


# ----------------------------------------------------------------------------- molecule
class Molecule:
    """Atoms, coordinates and the covalent bond graph."""

    def __init__(self, name, symbols, coords, charge=None, mult=None):
        self.name, self.symbols, self.coords = name, symbols, coords
        self.charge, self.mult = charge, mult
        self.n = len(symbols)
        self.neighbours = self._connect()

    @classmethod
    def read(cls, path: Path) -> Molecule:
        """Read an .xyz file or a Gaussian input (.gjf / .com / .inp)."""
        text = path.read_text(encoding="utf-8", errors="replace")
        charge = mult = None
        if path.suffix.lower() == ".xyz":
            lines = text.splitlines()
            atom_lines = lines[2:2 + int(lines[0].split()[0])]
        else:  # blocks: link0 + route | title | charge-multiplicity + atoms
            blocks = re.split(r"\n\s*\n", text.strip())
            if len(blocks) < 3:
                raise InputError(f"{path.name}: not a Gaussian input")
            head, *atom_lines = blocks[2].splitlines()
            charge, mult = (int(x) for x in head.split()[:2])
        symbols, coords = [], []
        for line in atom_lines:
            parts = line.split()
            symbol = re.match(r"[A-Z][a-z]?", parts[0].capitalize())
            if not symbol or symbol.group() not in ATOMIC_NUMBER:
                raise InputError(f"{path.name}: unknown element in line '{line.strip()}'")
            symbols.append(symbol.group())
            coords.append(tuple(float(x) for x in parts[-3:]))  # also skips a Gaussian freeze flag
        return cls(path.stem, symbols, coords, charge, mult)

    def _connect(self):
        bonded = [set() for _ in range(self.n)]
        for i, j in combinations(range(self.n), 2):
            ri = COVALENT_RADIUS.get(self.symbols[i], 1.5)
            rj = COVALENT_RADIUS.get(self.symbols[j], 1.5)
            if self.distance(i, j) <= BOND_TOLERANCE * (ri + rj):
                bonded[i].add(j)
                bonded[j].add(i)
        for h in range(self.n):  # a hydrogen keeps only its closest partner
            if self.symbols[h] == "H" and len(bonded[h]) > 1:
                keep = min(bonded[h], key=lambda j: self.distance(h, j))
                for j in bonded[h] - {keep}:
                    bonded[j].discard(h)
                bonded[h] = {keep}
        return bonded

    def label(self, i):
        return f"{self.symbols[i]}{i + 1}"

    def distance(self, i, j):
        return math.dist(self.coords[i], self.coords[j])

    def angle(self, i, j, k):
        u = [a - b for a, b in zip(self.coords[i], self.coords[j])]
        v = [a - b for a, b in zip(self.coords[k], self.coords[j])]
        cos = sum(a * b for a, b in zip(u, v)) / (math.hypot(*u) * math.hypot(*v))
        return math.degrees(math.acos(max(-1.0, min(1.0, cos))))

    def path(self, start, goal):
        """Shortest covalent path start -> goal (list of atoms), or None."""
        previous, queue = {start: None}, deque([start])
        while queue:
            atom = queue.popleft()
            if atom == goal:
                route = []
                while atom is not None:
                    route.append(atom)
                    atom = previous[atom]
                return route[::-1]
            for nxt in sorted(self.neighbours[atom] - previous.keys()):
                previous[nxt] = atom
                queue.append(nxt)
        return None

    def side(self, start, cut):
        """Atoms reachable from `start` without crossing the bond `cut`."""
        cut, seen, queue = set(cut), {start}, deque([start])
        while queue:
            atom = queue.popleft()
            for nxt in self.neighbours[atom] - seen:
                if {atom, nxt} != cut:
                    seen.add(nxt)
                    queue.append(nxt)
        return frozenset(seen)


# ----------------------------------------------------------------------------- hydrogen bonds
@dataclass
class HBond:
    index: int
    donor: int
    hydrogen: int
    acceptor: int
    path: list[int]                       # covalent path donor -> acceptor
    donor_group: frozenset = frozenset()
    acceptor_group: frozenset = frozenset()
    fd: int = 0                           # fragment number of the donor group
    fa: int = 0                           # fragment number of the acceptor group
    partners: list[int] = field(default_factory=list)

    def name(self, mol):
        return f"{mol.label(self.donor)}-{mol.label(self.hydrogen)}...{mol.label(self.acceptor)}"

    def geometry(self, mol):
        d, h, a = self.donor, self.hydrogen, self.acceptor
        return mol.distance(d, h), mol.distance(h, a), mol.distance(d, a), mol.angle(d, h, a)


def has_lone_pair(mol, atom):
    symbol, nbrs = mol.symbols[atom], mol.neighbours[atom]
    nitro = symbol == "N" and sum(mol.symbols[j] == "O" for j in nbrs) >= 2
    return len(nbrs) <= MAX_NEIGHBOURS.get(symbol, 8) and not nitro


def find_hbonds(mol, s):
    (dmin, dmax), (amin, amax) = s.hb_distance, s.hb_angle
    found = []
    for h in range(mol.n):
        if mol.symbols[h] != "H" or len(mol.neighbours[h]) != 1:
            continue
        (d,) = mol.neighbours[h]
        if mol.symbols[d] not in s.donors:
            continue
        for a in range(mol.n):
            if a in (d, h) or mol.symbols[a] not in s.acceptors or not has_lone_pair(mol, a):
                continue
            if not (dmin <= mol.distance(h, a) <= dmax and amin <= mol.angle(d, h, a) <= amax):
                continue
            path = mol.path(d, a)
            if path and len(path) - 1 >= s.min_bonds:
                found.append(HBond(len(found) + 1, d, h, a, path))
    return found


def parse_atoms(spec, n):
    """'1,3,5-7' (1-based) -> frozenset of 0-based indices."""
    atoms = set()
    for token in filter(None, re.split(r"[,\s]+", spec)):
        lo, _, hi = token.partition("-")
        atoms.update(range(int(lo) - 1, int(hi or lo)))
    if not all(0 <= i < n for i in atoms):
        raise InputError(f"atom number out of range in '{spec}'")
    return frozenset(atoms)


def manual_hbonds(mol, s):
    """hbN = D H A [dgroup=..] [agroup=..]   (1-based atom numbers)."""
    hbs = []
    for n, spec in sorted(s.manual.items()):
        tokens = spec.split()
        try:
            d, h, a = (int(t) - 1 for t in tokens[:3])
        except ValueError:
            raise InputError(f"hb{n}: expected 'D H A' atom numbers, got '{spec}'") from None
        if not all(0 <= i < mol.n for i in (d, h, a)):
            raise InputError(f"hb{n}: atom number out of range 1..{mol.n}")
        if mol.symbols[h] != "H" or d not in mol.neighbours[h]:
            raise InputError(f"hb{n}: atom {h + 1} is not a hydrogen bonded to atom {d + 1}")
        path = mol.path(d, a)
        if path is None:
            raise InputError(f"hb{n}: donor and acceptor are not covalently connected")
        hb = HBond(n, d, h, a, path)
        for token in tokens[3:]:
            key, _, value = token.partition("=")
            if key.lower() not in ("dgroup", "agroup"):
                raise InputError(f"hb{n}: unknown option '{token}'")
            setattr(hb, "donor_group" if key.lower() == "dgroup" else "acceptor_group",
                    parse_atoms(value, mol.n))
        hbs.append(hb)
    return hbs


# ----------------------------------------------------------------------------- fragments
def cut_group(mol, walk, skip=0):
    """Group on walk[0]'s side of the first non-ring bond along `walk` (after `skip` of them)."""
    for u, v in zip(walk, walk[1:]):
        side = mol.side(walk[0], (u, v))
        if walk[-1] in side:      # bond lies in a ring
            continue
        if skip == 0:
            return side
        skip -= 1
    return None


def donor_group(mol, hb):
    """D-H plus whatever hangs on D:  -OH -> OH, ring N-H -> the ring."""
    return cut_group(mol, hb.path)


def acceptor_group(mol, hb):
    """A terminal acceptor without H (C=O, N=O) takes its carrier atom along (-CHO, -COOH);
    otherwise the group is cut at the first non-ring bond (OH, OMe, ring)."""
    walk, a = hb.path[::-1], hb.acceptor
    terminal = len(mol.neighbours[a]) == 1 and mol.symbols[next(iter(mol.neighbours[a]))] != "H"
    return (terminal and cut_group(mol, walk, skip=1)) or cut_group(mol, walk)


def assign_fragments(mol, hbs):
    """Set donor / acceptor groups and number the distinct groups F1, F2, ..."""
    numbers = {}
    for hb in hbs:
        tag = f"HB{hb.index} ({hb.name(mol)})"
        hb.donor_group = hb.donor_group or donor_group(mol, hb)
        hb.acceptor_group = hb.acceptor_group or acceptor_group(mol, hb)
        if not hb.donor_group or not hb.acceptor_group:
            raise InputError(f"{tag}: no bond separates donor and acceptor; give dgroup=/agroup=")
        if not {hb.donor, hb.hydrogen} <= hb.donor_group or hb.acceptor not in hb.acceptor_group:
            raise InputError(f"{tag}: donor group needs D and H, acceptor group needs A")
        if hb.donor_group & hb.acceptor_group:
            raise InputError(f"{tag}: donor and acceptor groups overlap")
        hb.fd = numbers.setdefault(hb.donor_group, len(numbers) + 1)
        hb.fa = numbers.setdefault(hb.acceptor_group, len(numbers) + 1)
    return {k: atoms for atoms, k in numbers.items()}


def link_partners(mol, hbs, max_bonds):
    """Two HBs form a network if one's D/A is within `max_bonds` bonds of the other's."""
    for p, q in combinations(hbs, 2):
        routes = [mol.path(x, y) for x in (p.donor, p.acceptor) for y in (q.donor, q.acceptor)]
        gap = min((len(r) - 1 for r in routes if r), default=math.inf)
        if gap <= max_bonds:
            p.partners.append(q.index)
            q.partners.append(p.index)


@dataclass
class Job:
    """One capped structure M-Fi-Fj... = one Gaussian calculation."""
    name: str
    labels: tuple[int, ...]
    removed: frozenset
    atoms: list[tuple[str, tuple]]
    caps: list[str]
    gjf: Path | None = None
    log: Path | None = None
    energy: float | None = None
    status: str = "pending"


class Tailor:
    """Builds capped structures; structures with the same removed atoms are one job."""

    def __init__(self, mol, groups, cap_distance):
        self.mol, self.groups, self.cap = mol, groups, cap_distance
        self.jobs: dict[frozenset, Job] = {}

    def job(self, labels) -> Job:
        labels = tuple(sorted(set(labels)))
        removed = frozenset().union(*(self.groups[k] for k in labels))
        if removed not in self.jobs:
            name = "-".join(["M"] + [f"F{k}" for k in labels])
            self.jobs[removed] = Job(name, labels, removed, *self._capped(removed, name))
        return self.jobs[removed]

    def _capped(self, removed, name):
        """Kept atoms in parent order; each removed atom is replaced by an H cap on every
        kept neighbour, placed along the broken bond."""
        mol, atoms, caps = self.mol, [], []
        for i in range(mol.n):
            if i not in removed:
                atoms.append((mol.symbols[i], mol.coords[i]))
                continue
            for k in sorted(mol.neighbours[i] - removed):
                if mol.symbols[k] == "H":
                    raise InputError(f"{name}: {mol.label(k)} would lose its only bond")
                bond = [b - a for a, b in zip(mol.coords[k], mol.coords[i])]
                scale = self.cap / math.hypot(*bond)
                atoms.append(("H", tuple(a + scale * b for a, b in zip(mol.coords[k], bond))))
                caps.append(f"H{len(atoms)} on {mol.label(k)}")
        return atoms, caps


@dataclass
class Scheme:
    """The four structures giving E_HB in M, and the same four in M' (network broken)."""
    hb: HBond
    full: tuple[Job, ...]
    reduced: tuple[Job, ...] | None
    broken: tuple[int, ...]

    @staticmethod
    def energy(quartet):
        if quartet is None or any(j.energy is None for j in quartet):
            return None
        m, d, a, da = (j.energy for j in quartet)
        return (d + a - da - m) * HARTREE_TO_KCAL

    @property
    def e_full(self):
        return self.energy(self.full)

    @property
    def e_reduced(self):
        return self.energy(self.reduced)

    @property
    def e_coop(self):
        full, reduced = self.e_full, self.e_reduced
        return None if full is None or reduced is None else full - reduced


def build_scheme(hb, hbs, groups, tailor):
    def quartet(base):
        return tuple(tailor.job(base + extra) for extra in ([], [hb.fd], [hb.fa], [hb.fd, hb.fa]))

    protected = set(hb.path) | {hb.hydrogen}  # never cut through the bond being measured
    broken = sorted({k for p in hbs if p.index in hb.partners
                     for k in (p.fd, p.fa) if not groups[k] & protected})
    return Scheme(hb, quartet([]), quartet(broken) if broken else None, tuple(broken))


# ----------------------------------------------------------------------------- Gaussian
def write_gjf(job, s, charge, mult, title):
    lines = [f"%nprocshared={s.nproc}", f"%mem={s.mem}", s.route, "", title, "", f"{charge} {mult}"]
    lines += [f"{sym:<2} {x:16.9f} {y:16.9f} {z:16.9f}" for sym, (x, y, z) in job.atoms]
    job.gjf.write_text("\n".join(lines) + "\n\n", encoding="ascii")


def read_log(path):
    """-> (normal termination?, {level: energy in Hartree})"""
    if path is None or not path.is_file():
        return False, {}
    text = path.read_text(encoding="utf-8", errors="replace")
    energies = {}
    for level, pattern in ENERGY_PATTERNS:
        hits = pattern.findall(text)
        if hits:
            energies[level] = float(hits[-1].upper().replace("D", "E"))
    return "Normal termination" in text, energies


def load_energies(jobs):
    """Fill job.energy with the highest level found in every finished log."""
    finished = {}
    for job in jobs:
        normal, energies = read_log(job.log)
        if job.log is None or not job.log.is_file():
            job.status = "no log"
        elif not normal:
            job.status = "failed"
        else:
            job.status, finished[job.name] = "ok", (job, energies)
    level = next((lv for lv, _ in ENERGY_PATTERNS
                  if finished and all(lv in e for _, e in finished.values())), None)
    for job, energies in finished.values():
        job.energy = energies.get(level)
    return level


def run_jobs(jobs, s, workdir):
    pending = [job for job in jobs if not read_log(job.log)[0]]
    print(f"  running {len(pending)} Gaussian job(s), {s.njob} at a time")

    def launch(job):
        command = s.command.format(input=job.gjf.name)
        return job, subprocess.run(command, shell=True, cwd=workdir).returncode

    with ThreadPoolExecutor(s.njob) as pool:
        for job, code in pool.map(launch, pending):
            print(f"    {job.gjf.name}: {'done' if code == 0 else f'FAILED (exit {code})'}")


def shape(atoms):
    """Fingerprint independent of atom order, rotation and reflection."""
    return (sorted(sym for sym, _ in atoms),
            sorted(math.dist(p, q) for (_, p), (_, q) in combinations(atoms, 2)))


def map_reference(jobs, ref_dir):
    """Point every job to the log of a geometrically identical input in `ref_dir`."""
    references = []
    for inp in sorted(ref_dir.iterdir()):
        log = inp.with_suffix(".log")
        if inp.suffix.lower() in (".gjf", ".com", ".inp") and read_log(log)[0]:
            ref = Molecule.read(inp)
            references.append((shape(list(zip(ref.symbols, ref.coords))), log))
    for job in jobs:
        syms, dists = shape(job.atoms)
        job.log = next((log for (rs, rd), log in references
                        if rs == syms and all(abs(a - b) < 1e-3 for a, b in zip(rd, dists))), None)


# ----------------------------------------------------------------------------- report
@dataclass
class Analysis:
    mol: Molecule
    source: Path
    settings: Settings
    charge: int
    mult: int
    hbs: list[HBond]
    groups: dict[int, frozenset]
    schemes: list[Scheme]
    jobs: list[Job]
    level: str | None = None
    notes: list[str] = field(default_factory=list)


def fmt(value, digits=3):
    return "-" if value is None else f"{value:.{digits}f}"


def report(an: Analysis, mode: str, details=True) -> str:
    mol, s = an.mol, an.settings
    atoms_text = lambda atoms: " ".join(mol.label(i) for i in sorted(atoms))
    roles = {k: [] for k in an.groups}
    for hb in an.hbs:
        roles[hb.fd].append(f"HB{hb.index} donor")
        roles[hb.fa].append(f"HB{hb.index} acceptor")

    out = [
        f"IntraHB {__version__}  -  intramolecular H-bond energies by MTA        developed by {__author__}",
        "based on: Deshmukh, Bartolotti, Gadre, J. Phys. Chem. A 2008, 112, 312-321",
        "=" * 100,
        f"Molecule   {mol.name}   ({an.source.name}, {mol.n} atoms, charge {an.charge}, multiplicity {an.mult})",
        f"Theory     {s.route}",
        f"Energies   {an.level or 'not available yet'} (kcal/mol)      mode {mode}      {datetime.now():%Y-%m-%d %H:%M}",
        "",
        f"FRAGMENTS  (cut out of M, open valences capped with H at {s.cap_distance:.2f} A)",
    ]
    out += [f"  F{k:<4}{atoms_text(atoms):<26}{', '.join(roles[k])}" for k, atoms in an.groups.items()]

    out += ["", "HYDROGEN BONDS",
            f"  {'HB':<6}{'D-H...A':<18}{'r(D-H)':>8}{'r(H..A)':>9}{'r(D..A)':>9}{'angle':>9}"
            f"{'E_HB':>10}{"E_HB(M')":>10}{'E_coop':>9}   network"]
    for sc in an.schemes:
        hb = sc.hb
        dh, ha, da, ang = hb.geometry(mol)
        partners = ", ".join(f"HB{p}" for p in hb.partners) or "-"
        out.append(f"  {'HB%d' % hb.index:<6}{hb.name(mol):<18}{dh:8.3f}{ha:9.3f}{da:9.3f}{ang:9.2f}"
                   f"{fmt(sc.e_full):>10}{fmt(sc.e_reduced):>10}{fmt(sc.e_coop):>9}   {partners}")

    for sc in an.schemes:
        hb = sc.hb
        out += ["", f"HB{hb.index}  {hb.name(mol)}   donor F{hb.fd}, acceptor F{hb.fa}" +
                (f",   M' = M-{'-'.join(f'F{k}' for k in sc.broken)}" if sc.reduced else "")]
        if sc.reduced:
            out.append(f"  {'':<12}{'in M':<32}{"in M' (network broken)"}")
        for i, term in enumerate(("E(M)", "E(M-Fd)", "E(M-Fa)", "E(M-Fd-Fa)")):
            row = f"  {term:<12}{sc.full[i].name:<16}{fmt(sc.full[i].energy, 8):>16}"
            if sc.reduced:
                row += f"    {sc.reduced[i].name:<18}{fmt(sc.reduced[i].energy, 8):>16}"
            out.append(row)
        row = f"  {'E_HB':<12}{'':<16}{fmt(sc.e_full):>16}"
        if sc.reduced:
            row += f"    {'':<18}{fmt(sc.e_reduced):>16}    E_coop = {fmt(sc.e_coop)}"
        elif hb.partners:
            row += "    (network cannot be broken without cutting this HB)"
        else:
            row += "    (isolated HB, no cooperativity)"
        out.append(row)

    if an.notes:
        out += ["", "NOTES"] + [f"  ! {note}" for note in an.notes]
    if not details:
        return "\n".join(out) + "\n"

    out += ["", "JOBS", f"  {'name':<20}{'atoms':>6}{'caps':>6}  {'status':<9}{'energy (Ha)':>17}   file"]
    for job in an.jobs:
        out.append(f"  {job.name:<20}{len(job.atoms):>6}{len(job.caps):>6}  {job.status:<9}"
                   f"{fmt(job.energy, 8):>17}   {job.log.name if job.log else '-'}")

    out += ["", "COORDINATES  (Angstrom)"]
    for job in an.jobs:
        removed = f"removed {atoms_text(job.removed)};  caps {', '.join(job.caps)}" if job.removed else "parent"
        out += ["", f"# {job.name}   {removed}", str(len(job.atoms)), f"{mol.name} {job.name}"]
        out += [f"{sym:<2} {x:16.9f} {y:16.9f} {z:16.9f}" for sym, (x, y, z) in job.atoms]
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------------- driver
def analyse(path: Path, s: Settings) -> Analysis:
    mol = Molecule.read(path)
    charge = s.charge if s.charge is not None else (mol.charge or 0)
    electrons = lambda symbols: sum(ATOMIC_NUMBER[x] for x in symbols) - charge
    mult = s.mult or mol.mult or (1 if electrons(mol.symbols) % 2 == 0 else 2)

    hbs = manual_hbonds(mol, s) if s.manual else find_hbonds(mol, s)
    if not hbs:
        raise InputError(f"{path.name}: no intramolecular hydrogen bond found with the current thresholds")
    groups = assign_fragments(mol, hbs)
    link_partners(mol, hbs, s.link_bonds)

    tailor = Tailor(mol, groups, s.cap_distance)
    schemes = [build_scheme(hb, hbs, groups, tailor) for hb in hbs]
    jobs = sorted(tailor.jobs.values(), key=lambda j: (len(j.labels), j.labels))

    an = Analysis(mol, path, s, charge, mult, hbs, groups, schemes, jobs)
    for job in jobs:
        if electrons(sym for sym, _ in job.atoms) % 2 != (mult - 1) % 2:
            an.notes.append(f"{job.name}: electron count does not fit multiplicity {mult}")
    for h in {hb.hydrogen for hb in hbs}:
        shared = [f"HB{hb.index}" for hb in hbs if hb.hydrogen == h]
        if len(shared) > 1:
            an.notes.append(f"{mol.label(h)} is shared by {', '.join(shared)} (bifurcated HB)")
    return an


def process(path: Path, s: Settings, mode: str, ref: Path | None):
    an = analyse(path, s)
    print(f"\n=== {path.name}: {len(an.hbs)} HB(s), {len(an.groups)} fragment(s), {len(an.jobs)} job(s)")
    if mode == "detect":
        print(report(an, mode, details=False))
        return

    workdir = s.outdir / an.mol.name
    workdir.mkdir(parents=True, exist_ok=True)
    for job in an.jobs:
        job.gjf = workdir / f"{an.mol.name}_{job.name}.gjf"
        job.log = job.gjf.with_suffix(".log")
        if mode in ("gen", "run"):
            removed = " ".join(an.mol.label(i) for i in sorted(job.removed)) or "none"
            write_gjf(job, s, an.charge, an.mult, f"{an.mol.name} {job.name}  removed: {removed}")
    if mode == "run":
        run_jobs(an.jobs, s, workdir)
    if ref:
        map_reference(an.jobs, ref)
    if mode in ("run", "cal"):
        an.level = load_energies(an.jobs)

    summary = workdir / f"{an.mol.name}_summary.txt"
    summary.write_text(report(an, mode), encoding="utf-8")
    print(report(an, mode, details=False))
    print(f"  written: {summary}" + (f"  (+ {len(an.jobs)} .gjf inputs)" if mode in ("gen", "run") else ""))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Intramolecular H-bond energies and cooperativity by MTA")
    ap.add_argument("config", nargs="?", default="INPUT.cfg", type=Path, help="config file (default INPUT.cfg)")
    ap.add_argument("--mode", choices=MODES, help="override the mode set in the config")
    ap.add_argument("--ref", type=Path, help="folder with finished Gaussian inputs + logs to take energies from")
    args = ap.parse_args(argv)

    try:
        if not args.config.is_file():
            raise InputError(f"config file '{args.config}' not found")
        s = load_settings(args.config)
        if not s.structures:
            raise InputError("no structure given (xyz_file in the config)")
        if args.ref and not args.ref.is_dir():
            raise InputError(f"--ref: '{args.ref}' is not a folder")
    except InputError as e:
        print(f"ERROR: {e}")
        return 1

    status = 0
    for path in s.structures:
        try:
            if not path.is_file():
                raise InputError(f"structure file '{path}' not found")
            process(path, s, args.mode or s.mode, args.ref)
        except InputError as e:
            print(f"ERROR: {e}")
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
