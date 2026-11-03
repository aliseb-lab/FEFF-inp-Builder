#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FEFF .inp Builder
=================
Turns crystal and molecule structure files into FEFF6 feff.inp files for Artemis.

Reads   .cif / .mmcif       crystal structures (symmetry is expanded for you)
        .res / .ins         SHELX files
        .xyz                Mercury / Avogadro / ORCA style xyz
        .pdb / .ent         PDB (first model, altloc A)
        .mol / .sdf         MDL molfiles (V2000 and V3000)
        .mol2               Tripos mol2
        .gjf / .com         Gaussian input (Cartesian)
        POSCAR / CONTCAR / .vasp
        .inp                an existing feff.inp (to re-center or re-cut it)

Run     Double-click this file (Windows: Python 3.8+ from python.org).
        Or from a terminal:
          python feff_inp_builder.pyw                        open the window
          python feff_inp_builder.pyw a.cif b.xyz            open with files loaded
          python feff_inp_builder.pyw --batch a.cif -a Co    no window, just write .inp files
          python feff_inp_builder.pyw --batch -h             all batch options

        (On macOS/Linux use python3; on Windows, py works if python isn't found.)

Only the Python standard library is used. Optional: `pip install tkinterdnd2`
turns on drag-and-drop of files onto the window.

Author:  Ali Seblini
License: MIT (see LICENSE)
"""

import math
import os
import re
import sys

APP_NAME = "FEFF .inp Builder"
VERSION = "1.0"

# --------------------------------------------------------------------------
# Periodic table data
# --------------------------------------------------------------------------

SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te "
    "I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir "
    "Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No "
    "Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
Z_OF = {s.upper(): i for i, s in enumerate(SYMBOLS) if i}

# Covalent radii (Angstrom), Cordero et al. 2008, indexed by Z (0 = unused).
COV = [0.0,
       0.31, 0.28,
       1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
       1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
       2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22,
       1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
       2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
       1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
       2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92,
       1.92, 1.89, 1.90, 1.87, 1.87,
       1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48,
       1.40, 1.50, 1.50,
       2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69]

NONMETALS = {1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36,
             51, 52, 53, 54, 85, 86}

BOND_TOL = 0.45      # bonded if d <= r_cov(A) + r_cov(B) + BOND_TOL (Mercury-like)
MIN_BOND = 0.40      # anything closer than this is an overlap, not a bond
EDGES = {"K": 1, "L1": 2, "L2": 3, "L3": 4}
FEFF6_MAX_ATOMS = 500   # FEFF 6L stops with "Too many atoms in ATOMS list" above this
KEEP_PARTS = {"", ".", "?", "0", "1", "-1", "a", "A"}


def zof(el):
    return Z_OF.get(el.upper(), 0)


def cov(el):
    z = zof(el)
    return COV[z] if 0 < z < len(COV) else 1.50


def _proper(sym):
    z = Z_OF.get(sym.upper())
    return SYMBOLS[z] if z else None


def element_from_label(lab):
    """Element from an atom label like 'Co1', 'C18A', 'N5', 'Cl2', 'H1A'."""
    m = re.match(r"\s*([A-Za-z])([A-Za-z]?)", str(lab))
    if not m:
        return None
    a, b = m.group(1), m.group(2)
    if b and b.islower() and _proper(a + b):
        return _proper(a + b)
    if _proper(a):
        return _proper(a)
    if b and _proper(a + b):
        return _proper(a + b)
    return None


def element_from_symbol(tok):
    """Element from a symbol-ish token: 'Co', 'CO', 'Co2+', '27', 'C.ar'."""
    tok = str(tok).strip()
    if not tok:
        return None
    if tok.isdigit():
        z = int(tok)
        return SYMBOLS[z] if 0 < z < len(SYMBOLS) else None
    m = re.match(r"[A-Za-z]+", tok)
    if not m:
        return None
    alpha = m.group(0)
    if alpha.upper() in ("D", "T"):
        return "H"
    if alpha == tok and _proper(alpha):       # a bare symbol, any case
        return _proper(alpha)
    return element_from_label(tok)


_NUM = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")


def num(s):
    """Float from a CIF-style number, dropping uncertainties: '1.234(5)' -> 1.234."""
    if s is None:
        return None
    m = _NUM.match(str(s))
    return float(m.group(1)) if m else None


def is_int(s):
    return re.match(r"^[-+]?\d+$", s) is not None


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# --------------------------------------------------------------------------
# Vector helpers
# --------------------------------------------------------------------------

def dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def vnorm(u):
    return math.sqrt(dot(u, u))


def lattice_from_params(a, b, c, al, be, ga):
    al, be, ga = math.radians(al), math.radians(be), math.radians(ga)
    ca, cb, cg, sg = math.cos(al), math.cos(be), math.cos(ga), math.sin(ga)
    cx = c * cb
    cy = c * (ca - cb * cg) / sg
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return ((a, 0.0, 0.0), (b * cg, b * sg, 0.0), (cx, cy, cz))


def f2c(f, L):
    return (f[0] * L[0][0] + f[1] * L[1][0] + f[2] * L[2][0],
            f[0] * L[0][1] + f[1] * L[1][1] + f[2] * L[2][1],
            f[0] * L[0][2] + f[1] * L[1][2] + f[2] * L[2][2])


def c2f(c, L):
    """Cartesian -> fractional for lattice rows L."""
    a, b, cc = L
    V = dot(a, cross(b, cc))
    ra, rb, rc = cross(b, cc), cross(cc, a), cross(a, b)
    return (dot(c, ra) / V, dot(c, rb) / V, dot(c, rc) / V)


def recip_lengths(L):
    a, b, c = L
    V = abs(dot(a, cross(b, c)))
    return (vnorm(cross(b, c)) / V, vnorm(cross(c, a)) / V, vnorm(cross(a, b)) / V)


def _cell_key(p, size):
    return (int(math.floor(p[0] / size)), int(math.floor(p[1] / size)), int(math.floor(p[2] / size)))


def _grid(points, size):
    g = {}
    for i, p in enumerate(points):
        g.setdefault(_cell_key(p, size), []).append(i)
    return g


def _near(g, p, size):
    kx, ky, kz = _cell_key(p, size)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for j in g.get((kx + dx, ky + dy, kz + dz), ()):
                    yield j


# --------------------------------------------------------------------------
# Symmetry operators
# --------------------------------------------------------------------------

def _frac_val(s):
    s = s.strip()
    if not s:
        return 1.0
    if "/" in s:
        n, d = s.split("/", 1)
        return float(n or 1) / float(d)
    return float(s)


def parse_symop(text):
    """'-x+1/2, y+1/2, -z+1/2' -> (rotation rows, translation)."""
    s = text.replace(" ", "").replace("'", "").replace('"', "").lower()
    parts = s.split(",")
    if len(parts) != 3:
        raise ValueError("bad symmetry operator: %r" % text)
    R, T = [], []
    for p in parts:
        row, t = [0.0, 0.0, 0.0], 0.0
        for term in re.findall(r"[+-]?[^+-]+", p):
            sign = -1.0 if term.startswith("-") else 1.0
            term = term.lstrip("+-")
            var = None
            for v in "xyz":
                if v in term:
                    var = v
                    term = term.replace(v, "").replace("*", "")
            if var:
                row["xyz".index(var)] += sign * _frac_val(term)
            else:
                t += sign * _frac_val(term)
        R.append(tuple(row))
        T.append(t)
    return (tuple(R), tuple(T))


IDENTITY = (((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)), (0.0, 0.0, 0.0))


def apply_op(op, f):
    R, T = op
    return tuple(R[k][0] * f[0] + R[k][1] * f[1] + R[k][2] * f[2] + T[k] for k in range(3))


def expand_sites(sites, ops, L, tol=0.1):
    """Asymmetric unit -> full unit cell. Returns [el, label, frac, part, site_index]."""
    out = []
    tol2 = tol * tol
    for si, (el, lab, f, part) in enumerate(sites):
        pts = []
        for op in ops:
            g = apply_op(op, f)
            g = tuple(v - math.floor(v) for v in g)
            dup = False
            for p in pts:
                d = [g[k] - p[k] for k in range(3)]
                d = [v - round(v) for v in d]
                c = f2c(d, L)
                if dot(c, c) < tol2:
                    dup = True
                    break
            if not dup:
                pts.append(g)
        for p in pts:
            out.append((el, lab, p, part, si))
    return out


# --------------------------------------------------------------------------
# Structure container
# --------------------------------------------------------------------------

class Structure(object):
    """Either molecular (Cartesian atoms) or periodic (lattice + asym sites)."""

    def __init__(self, name, path=""):
        self.name = name
        self.path = path
        self.atoms = []        # molecular: (el, label, (x,y,z), part)
        self.lattice = None    # periodic: 3 row vectors (A)
        self.sites = []        # periodic: (el, label, (fa,fb,fc), part)
        self.ucell = []        # periodic: (el, label, frac, part, site_index)
        self.notes = []
        self.default_abs = None
        self._bonds = None

    @property
    def periodic(self):
        return self.lattice is not None

    def entries(self):
        return self.sites if self.periodic else self.atoms

    def elements(self):
        return sorted({e[0] for e in self.entries()}, key=lambda e: (-zof(e), e))

    def default_element(self):
        ents = self.entries()
        if self.default_abs is not None and self.default_abs < len(ents):
            return ents[self.default_abs][0]
        els = self.elements()
        metals = [e for e in els if zof(e) not in NONMETALS]
        return metals[0] if metals else (els[0] if els else None)

    def indices_of(self, el):
        return [i for i, e in enumerate(self.entries()) if e[0] == el]

    def default_index(self, el):
        idx = self.indices_of(el)
        if self.default_abs in idx:
            return self.default_abs
        return idx[0] if idx else None

    def label_of(self, i):
        return self.entries()[i][1]

    def describe(self):
        n = len(self.entries())
        if self.periodic:
            return "%d sites, %d in cell" % (n, len(self.ucell))
        return "%d atoms" % n


def _finish_periodic(st, sites, ops, cellpars=None, L=None):
    if L is None:
        L = lattice_from_params(*cellpars)
    st.lattice = L
    st.sites = sites
    st.ucell = expand_sites(sites, ops, L)


def _autolabel(atoms):
    """Give unlabeled atoms Mercury-style labels: C1, C2, N1 ..."""
    count = {}
    out = []
    for el, lab, p, part in atoms:
        if not lab:
            count[el] = count.get(el, 0) + 1
            lab = "%s%d" % (el, count[el])
        out.append((el, lab, p, part))
    return out


# --------------------------------------------------------------------------
# CIF / mmCIF
# --------------------------------------------------------------------------

def _cif_tokens(text):
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < n and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield ("\n".join(buf), True)
            continue
        j, m = 0, len(line)
        while j < m:
            ch = line[j]
            if ch.isspace():
                j += 1
                continue
            if ch == "#":
                break
            if ch in "'\"":
                k = j + 1
                while k < m and not (line[k] == ch and (k + 1 == m or line[k + 1].isspace())):
                    k += 1
                yield (line[j + 1:k], True)
                j = k + 1
                continue
            k = j
            while k < m and not line[k].isspace():
                k += 1
            yield (line[j:k], False)
            j = k
        i += 1


def _cif_norm(tag):
    return tag.lower().replace(".", "_")


def parse_cif_blocks(text):
    toks = list(_cif_tokens(text))
    blocks, cur = [], None
    i, n = 0, len(toks)

    def new_block(name):
        b = {"name": name, "tags": {}, "loops": []}
        blocks.append(b)
        return b

    while i < n:
        t, quoted = toks[i]
        tl = t.lower()
        if not quoted and tl.startswith("data_"):
            cur = new_block(t[5:])
            i += 1
            continue
        if cur is None:
            cur = new_block("")
        if not quoted and tl == "loop_":
            i += 1
            names = []
            while i < n and not toks[i][1] and toks[i][0].startswith("_"):
                names.append(_cif_norm(toks[i][0]))
                i += 1
            vals = []
            while i < n:
                t2, q2 = toks[i]
                t2l = t2.lower()
                if not q2 and (t2.startswith("_") or t2l == "loop_" or
                               t2l.startswith(("data_", "save_", "global_", "stop_"))):
                    break
                vals.append(t2)
                i += 1
            if names:
                loop = {nm: [] for nm in names}
                nc = len(names)
                for r in range(len(vals) // nc):
                    for c, nm in enumerate(names):
                        loop[nm].append(vals[r * nc + c])
                cur["loops"].append(loop)
            continue
        if not quoted and t.startswith("_"):
            val = toks[i + 1][0] if i + 1 < n else ""
            cur["tags"][_cif_norm(t)] = val
            i += 2
            continue
        i += 1
    return blocks


def _cif_get(block, *names):
    for nm in names:
        if nm in block["tags"]:
            return block["tags"][nm]
    for loop in block["loops"]:
        for nm in names:
            if nm in loop and loop[nm]:
                return loop[nm][0]
    return None


def _cif_list(block, *names):
    for loop in block["loops"]:
        for nm in names:
            if nm in loop:
                return loop[nm]
    for nm in names:
        if nm in block["tags"]:
            return [block["tags"][nm]]
    return None


def _cif_loop_with(block, key):
    for loop in block["loops"]:
        if key in loop:
            return loop
    return None


def _cif_block_structure(block, name):
    st = Structure(name)
    fl = _cif_loop_with(block, "_atom_site_fract_x")
    cl = _cif_loop_with(block, "_atom_site_cartn_x")
    cell = [num(_cif_get(block, t)) for t in (
        "_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")]
    skipped = 0
    if fl is not None and None in cell and cl is None:
        missing = [t for t, v in zip(("a", "b", "c", "alpha", "beta", "gamma"), cell) if v is None]
        raise ValueError("fractional coordinates but no cell parameters (missing %s)" % ", ".join(missing))
    if fl is not None and None not in cell:
        for k in ("_atom_site_fract_y", "_atom_site_fract_z"):
            if k not in fl:
                raise ValueError("the atom_site loop has no %s column" % k)
        nrow = len(fl["_atom_site_fract_x"])
        labels = fl.get("_atom_site_label") or [""] * nrow
        types = fl.get("_atom_site_type_symbol")
        parts = fl.get("_atom_site_disorder_group") or ["."] * nrow
        occ = fl.get("_atom_site_occupancy")
        sites, partial = [], 0
        for r in range(nrow):
            x = num(fl["_atom_site_fract_x"][r])
            y = num(fl["_atom_site_fract_y"][r])
            z = num(fl["_atom_site_fract_z"][r])
            el = element_from_symbol(types[r]) if types and types[r] not in ("?", ".") else None
            if el is None:
                el = element_from_label(labels[r])
            if el is None or None in (x, y, z):
                skipped += 1
                continue
            if (occ and num(occ[r]) is not None and num(occ[r]) < 0.999 and
                    parts[r] in (".", "?", "")):
                partial += 1
            sites.append((el, labels[r] or el, (x, y, z), parts[r]))
        ops_txt = _cif_list(block, "_symmetry_equiv_pos_as_xyz",
                            "_space_group_symop_operation_xyz",
                            "_space_group_symop_operation_xyz_")
        ops = []
        for o in ops_txt or []:
            try:
                ops.append(parse_symop(o))
            except Exception:
                st.notes.append("Ignored unreadable symmetry operator %r." % o)
        if not ops:
            sg = _cif_get(block, "_symmetry_space_group_name_h-m",
                          "_space_group_name_h-m_alt") or ""
            ops = [IDENTITY]
            if sg.replace(" ", "").upper() not in ("", "P1", "?"):
                st.notes.append("No symmetry operators in the CIF (space group %s); "
                                "using the listed atoms only." % sg)
        if partial:
            st.notes.append("%d site(s) have partial occupancy but no disorder group; "
                            "check them in the Atoms list." % partial)
        _finish_periodic(st, sites, ops, cellpars=cell)
    elif cl is not None:
        nrow = len(cl["_atom_site_cartn_x"])
        labels = (cl.get("_atom_site_label") or cl.get("_atom_site_auth_atom_id") or
                  cl.get("_atom_site_label_atom_id") or [""] * nrow)
        types = cl.get("_atom_site_type_symbol")
        alts = cl.get("_atom_site_label_alt_id")
        models = cl.get("_atom_site_pdbx_pdb_model_num")
        first_model = models[0] if models else None
        for r in range(nrow):
            if models and models[r] != first_model:
                continue
            if alts and alts[r] not in (".", "?", "A", "1"):
                continue
            x, y, z = (num(cl["_atom_site_cartn_%s" % k][r]) for k in "xyz")
            el = element_from_symbol(types[r]) if types and types[r] not in ("?", ".") else None
            if el is None:
                el = element_from_label(labels[r])
            if el is None or None in (x, y, z):
                skipped += 1
                continue
            st.atoms.append((el, labels[r], (x, y, z), "."))
    else:
        return None
    if skipped:
        st.notes.append("Skipped %d atom(s) with no recognizable element or coordinates." % skipped)
    return st


def parse_cif(text, stem):
    blocks = parse_cif_blocks(text)
    out = []
    for b in blocks:
        st = _cif_block_structure(b, stem)
        if st is not None and st.entries():
            st.name = b["name"] or stem
            out.append(st)
    if len(out) == 1:
        out[0].name = stem
    else:
        for st in out:
            st.name = "%s_%s" % (stem, st.name)
    return out


# --------------------------------------------------------------------------
# SHELX .res / .ins
# --------------------------------------------------------------------------

SHELX_CMDS = set((
    "TITL CELL ZERR LATT SYMM SFAC DISP UNIT LAUE REM MORE TIME END HKLF OMIT SHEL "
    "BASF TWIN EXTI SWAT HOPE MERG SPEC RESI MOVE ANIS AFIX HFIX FRAG FEND EXYZ EADP "
    "EQIV CONN PART BIND FREE DFIX DANG BUMP SAME SADI CHIV FLAT DELU SIMU DEFS ISOR "
    "NCSY SUMP L.S. CGLS BLOC DAMP STIR WGHT FVAR BOND CONF MPLA RTAB HTAB LIST ACTA "
    "SIZE TEMP WPDB FMAP GRID PLAN MOLE RIGU XNPD NEUT ABIN ANSC ANSR LONE PRIG WIGL "
    "BEDE SADM").split())

LATT_CENTERING = {
    1: [],
    2: [(0.5, 0.5, 0.5)],
    3: [(2 / 3., 1 / 3., 1 / 3.), (1 / 3., 2 / 3., 2 / 3.)],
    4: [(0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    5: [(0.0, 0.5, 0.5)],
    6: [(0.5, 0.0, 0.5)],
    7: [(0.5, 0.5, 0.0)],
}


def _shelx_coord(v):
    while v > 5:
        v -= 10.0
    while v < -5:
        v += 10.0
    return v


def parse_shelx(text, stem):
    lines, buf = [], ""
    for raw in text.splitlines():
        line = (buf + " " + raw.strip()) if buf else raw
        buf = ""
        if line.rstrip().endswith("="):
            buf = line.rstrip()[:-1]
            continue
        lines.append(line)
    cell, latt, symm, sfac, sites, part = None, 1, [], [], [], "0"
    st = Structure(stem)
    for line in lines:
        toks = line.split()
        if not toks:
            continue
        cmd = toks[0].upper()
        if cmd in ("END", "HKLF"):
            break
        if cmd == "CELL" and len(toks) >= 8:
            cell = [float(v) for v in toks[2:8]]
            continue
        if cmd == "LATT" and len(toks) >= 2:
            latt = int(toks[1])
            continue
        if cmd == "SYMM":
            symm.append(parse_symop(" ".join(toks[1:])))
            continue
        if cmd == "SFAC":
            if len(toks) > 2 and is_float(toks[2]):
                sfac.append(element_from_symbol(toks[1]))
            else:
                sfac.extend(element_from_symbol(t) for t in toks[1:])
            continue
        if cmd == "PART":
            part = toks[1] if len(toks) > 1 else "0"
            continue
        if cmd in SHELX_CMDS:
            continue
        if (len(toks) >= 5 and is_int(toks[1]) and is_float(toks[2]) and
                is_float(toks[3]) and is_float(toks[4])):
            if re.match(r"^Q\d+$", toks[0], re.I):
                continue
            k = int(toks[1])
            el = sfac[k - 1] if 0 < k <= len(sfac) else element_from_label(toks[0])
            if el is None:
                continue
            f = tuple(_shelx_coord(float(v)) for v in toks[2:5])
            sites.append((el, toks[0], f, part))
    if cell is None or not sites:
        return []
    ops = [IDENTITY] + symm
    if latt > 0:
        ops = ops + [(tuple(tuple(-v for v in row) for row in R), tuple(-t for t in T)) for R, T in ops]
    cen = LATT_CENTERING.get(abs(latt), [])
    ops = ops + [(R, tuple(T[k] + c[k] for k in range(3))) for c in cen for R, T in ops]
    _finish_periodic(st, sites, ops, cellpars=cell)
    return [st]


# --------------------------------------------------------------------------
# VASP POSCAR / CONTCAR
# --------------------------------------------------------------------------

def parse_poscar(text, stem):
    lines = [l for l in text.splitlines()]
    scale = float(lines[1].split()[0])
    L = [tuple(float(x) for x in lines[k].split()[:3]) for k in (2, 3, 4)]
    if scale < 0:
        V = abs(dot(L[0], cross(L[1], L[2])))
        scale = (abs(scale) / V) ** (1.0 / 3.0)
    L = tuple(tuple(v * scale for v in row) for row in L)
    toks = lines[5].split()
    if all(re.match(r"^[A-Za-z]", t) for t in toks):
        species = [element_from_symbol(t) for t in toks]
        counts = [int(t) for t in lines[6].split()]
        k = 7
    else:
        counts = [int(t) for t in toks]
        species = [element_from_symbol(t) for t in lines[0].split()[:len(counts)]]
        k = 6
    if lines[k].strip()[:1] in "sS":
        k += 1
    cart = lines[k].strip()[:1] in "cCkK"
    k += 1
    sites = []
    for el, cnt in zip(species, counts):
        for _ in range(cnt):
            v = [float(x) for x in lines[k].split()[:3]]
            k += 1
            if cart:
                v = c2f([x * scale for x in v], L)
            sites.append((el, "", tuple(v), "."))
    st = Structure(stem)
    _finish_periodic(st, _autolabel(sites), [IDENTITY], L=L)
    return [st]


# --------------------------------------------------------------------------
# Molecular formats
# --------------------------------------------------------------------------

def _xyz_atom(toks):
    if len(toks) < 4 or not all(is_float(t) for t in toks[1:4]):
        return None
    tok = toks[0]
    el = element_from_symbol(tok)
    if el is None:
        return None
    label = "" if (tok.isalpha() or tok.isdigit()) else tok   # keep real labels like C18
    return (el, label, tuple(float(t) for t in toks[1:4]), ".")


def parse_xyz(text, stem):
    lines = text.splitlines()
    atoms = []
    if lines and lines[0].split() and is_int(lines[0].split()[0]):
        n = int(lines[0].split()[0])
        body = lines[2:2 + n]
    else:
        body = lines
    for line in body:
        a = _xyz_atom(line.split())
        if a:
            atoms.append(a)
    st = Structure(stem)
    st.atoms = _autolabel(atoms)
    return [st] if st.atoms else []


def parse_loose(text, stem):
    """Any text with 'El x y z' (or Gaussian 'El 0 x y z') lines."""
    atoms = []
    for line in text.splitlines():
        toks = line.replace(",", " ").split()
        if len(toks) < 4:
            continue
        tok = toks[0].split("(")[0]
        alpha = re.sub(r"[^A-Za-z]", "", tok)
        if not (tok.isdigit() or (alpha and _proper(alpha) and tok[:len(alpha)] == alpha)):
            continue
        el = element_from_symbol(tok)
        if el is None:
            continue
        nums = toks[1:]

        def coords(ts):      # real coordinates have decimal points
            return all(is_float(t) and re.search(r"[.eE]", t) for t in ts)
        if len(nums) >= 4 and nums[0] in ("0", "-1") and coords(nums[1:4]):
            xyz = nums[1:4]
        elif coords(nums[:3]):
            xyz = nums[:3]
        else:
            continue
        atoms.append((el, "", tuple(float(t) for t in xyz), "."))
    st = Structure(stem)
    st.atoms = _autolabel(atoms)
    return [st] if st.atoms else []


def parse_gaussian(text, stem):
    """Gaussian input: Link 0 / route / blank / title / blank / charge-mult / atoms / blank."""
    lines = text.splitlines()
    i = 0
    while i < len(lines) and (lines[i].strip().startswith("%") or not lines[i].strip()):
        i += 1
    while i < len(lines) and lines[i].strip():          # route section (may span lines)
        i += 1
    i += 1
    while i < len(lines) and lines[i].strip():          # title
        i += 1
    i += 1
    if i >= len(lines) or len(lines[i].split()) < 2:
        raise ValueError("no charge/multiplicity line")
    i += 1
    atoms = []
    for line in lines[i:]:
        if not line.strip():
            break
        toks = line.replace(",", " ").split()
        el = element_from_symbol(toks[0].split("(")[0].split("-")[0])
        if el is None:
            continue
        nums = toks[1:]
        if len(nums) == 4 and nums[0] in ("0", "-1"):
            nums = nums[1:]                              # frozen-atom flag
        if len(nums) != 3 or not all(is_float(t) for t in nums):
            raise ValueError("this looks like Z-matrix input; save Cartesian coordinates instead")
        atoms.append((el, "", tuple(float(t) for t in nums), "."))
    st = Structure(stem)
    st.atoms = _autolabel(atoms)
    return [st] if st.atoms else []


def _pdb_element(line):
    el = line[76:78].strip() if len(line) >= 78 else ""
    if el:
        e = element_from_symbol(el)
        if e:
            return e
    s = re.sub(r"[^A-Za-z]", "", line[12:14])
    if len(s) == 2 and _proper(s):
        return _proper(s)
    return element_from_label(line[12:16].strip().lstrip("0123456789"))


def parse_pdb(text, stem):
    atoms = []
    for line in text.splitlines():
        rec = line[:6].strip().upper()
        if rec == "ENDMDL" and atoms:
            break
        if rec not in ("ATOM", "HETATM") or len(line) < 54:
            continue
        if line[16:17] not in (" ", "A", "1", ""):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        el = _pdb_element(line)
        if el is None:
            continue
        atoms.append((el, line[12:16].strip(), xyz, "."))
    st = Structure(stem)
    st.atoms = _autolabel(atoms)
    return [st] if st.atoms else []


def _molfile(lines):
    if len(lines) < 4:
        return None, []
    name = lines[0].strip()
    atoms = []
    if "V3000" in lines[3]:
        inside = False
        for line in lines[4:]:
            u = line.upper()
            if "BEGIN ATOM" in u:
                inside = True
                continue
            if "END ATOM" in u:
                break
            if inside and line.startswith("M  V30"):
                t = line.split()
                el = element_from_symbol(t[3])
                if el:
                    atoms.append((el, "", (float(t[4]), float(t[5]), float(t[6])), "."))
        return name, atoms
    n = int(lines[3][0:3])
    for line in lines[4:4 + n]:
        try:
            xyz = (float(line[0:10]), float(line[10:20]), float(line[20:30]))
            sym = line[31:34].strip()
        except ValueError:
            t = line.split()
            xyz = (float(t[0]), float(t[1]), float(t[2]))
            sym = t[3]
        el = element_from_symbol(sym)
        if el:
            atoms.append((el, "", xyz, "."))
    return name, atoms


def parse_mol(text, stem):
    text = text.replace("\r\n", "\n")
    out = []
    for rec in text.split("$$$$"):
        if rec.startswith("\n"):
            rec = rec[1:]
        if not rec.strip():
            continue
        try:
            name, atoms = _molfile(rec.split("\n"))
        except Exception:
            continue
        if atoms:
            st = Structure(stem)
            st.atoms = _autolabel(atoms)
            st._recname = name
            out.append(st)
    if len(out) > 1:
        for i, st in enumerate(out):
            st.name = "%s_%d" % (stem, i + 1)
    return out


def parse_mol2(text, stem):
    out, cur, section = [], None, ""
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip().upper()[9:]
            if section == "MOLECULE":
                cur = Structure(stem)
                out.append(cur)
            continue
        if section != "ATOM" or not line.strip():
            continue
        t = line.split()
        if len(t) < 6 or not all(is_float(v) for v in t[2:5]):
            continue
        el = element_from_symbol(t[5].split(".")[0])
        if el is None:
            continue
        if cur is None:
            cur = Structure(stem)
            out.append(cur)
        cur.atoms.append((el, t[1], (float(t[2]), float(t[3]), float(t[4])), "."))
    out = [s for s in out if s.atoms]
    for s in out:
        s.atoms = _autolabel(s.atoms)
    if len(out) > 1:
        for i, st in enumerate(out):
            st.name = "%s_%d" % (stem, i + 1)
    return out


def parse_feff_inp(text, stem):
    pots, atoms, section, abs_i = {}, [], "", None
    for raw in text.splitlines():
        line = raw.split("*")[0].split("!")[0].strip()
        if not line:
            continue
        t = line.split()
        key = t[0].upper()
        if key in ("POTENTIALS", "ATOMS"):
            section = key
            continue
        if key == "END":
            break
        if re.match(r"^[A-Z]", key):
            section = ""
            continue
        if section == "POTENTIALS" and len(t) >= 2 and is_int(t[0]) and is_int(t[1]):
            pots[int(t[0])] = SYMBOLS[int(t[1])]
        elif section == "ATOMS" and len(t) >= 4 and is_int(t[3]):
            ip = int(t[3])
            el = pots.get(ip)
            if el is None:
                continue
            atoms.append((el, t[4] if len(t) > 4 else "", (float(t[0]), float(t[1]), float(t[2])), "."))
            if ip == 0 and abs_i is None:
                abs_i = len(atoms) - 1          # the old absorber is the default again
    st = Structure(stem)
    st.atoms = _autolabel(atoms)
    st.default_abs = abs_i
    return [st] if st.atoms else []


# --------------------------------------------------------------------------
# Loading / format detection
# --------------------------------------------------------------------------

FILE_TYPES = (".cif", ".mmcif", ".mcif", ".res", ".ins", ".xyz", ".pdb", ".ent",
              ".mol", ".sdf", ".sd", ".mol2", ".gjf", ".com", ".gau", ".vasp",
              ".poscar", ".inp", ".txt")


def _sniff(text, stem, ext):
    base = stem.upper()
    if ext in (".cif", ".mmcif", ".mcif"):
        return "cif"
    if ext in (".res", ".ins"):
        return "shelx"
    if ext == ".xyz":
        return "xyz"
    if ext in (".pdb", ".ent"):
        return "pdb"
    if ext in (".mol", ".sdf", ".sd"):
        return "mol"
    if ext == ".mol2":
        return "mol2"
    if ext in (".gjf", ".com", ".gau"):
        return "gaussian"
    if ext in (".vasp", ".poscar") or "POSCAR" in base or "CONTCAR" in base:
        return "poscar"
    head = text[:20000]
    hl = head.splitlines()
    if (len(hl) > 7 and len(hl[1].split()) == 1 and is_float(hl[1].split()[0]) and
            all(len(hl[k].split()) >= 3 and all(is_float(t) for t in hl[k].split()[:3]) for k in (2, 3, 4))):
        return "poscar"
    if ext == ".inp" or ext == "":
        if re.search(r"^\s*POTENTIALS", head, re.M | re.I) and re.search(r"^\s*ATOMS", head, re.M | re.I):
            return "feff"
    if re.search(r"^data_", head, re.M):
        return "cif"
    if "@<TRIPOS>" in head:
        return "mol2"
    if re.search(r"^(ATOM  |HETATM)", head, re.M):
        return "pdb"
    if re.search(r"V[23]000", head):
        return "mol"
    if re.search(r"^CELL\s", head, re.M) and re.search(r"^SFAC\s", head, re.M):
        return "shelx"
    first = head.lstrip().split("\n", 1)[0].split()
    if first and is_int(first[0]):
        return "xyz"
    return "loose"


KIND_NAMES = {"cif": "CIF", "shelx": "SHELX", "xyz": "XYZ", "pdb": "PDB", "mol": "MDL molfile",
              "mol2": "Mol2", "poscar": "POSCAR", "feff": "feff.inp", "gaussian": "Gaussian input",
              "loose": "plain coordinates"}

PARSERS = {
    "cif": parse_cif, "shelx": parse_shelx, "xyz": parse_xyz, "pdb": parse_pdb,
    "mol": parse_mol, "mol2": parse_mol2, "poscar": parse_poscar,
    "feff": parse_feff_inp, "gaussian": parse_gaussian, "loose": parse_loose,
}


def load_structures(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    base = os.path.basename(path)
    stem, ext = os.path.splitext(base)
    kind = _sniff(text, stem, ext.lower())
    if kind == "feff" and stem.lower().endswith("_feff"):
        stem = stem[:-5]
    try:
        structs = PARSERS[kind](text, stem)
    except Exception as exc:
        raise ValueError("%s: could not read it as %s: %s" % (base, KIND_NAMES.get(kind, kind), exc))
    if not structs and kind in ("pdb", "mol", "mol2"):
        structs = parse_loose(text, stem)     # cell-based formats never fall back: coordinates would be wrong
    if not structs:
        raise ValueError("%s: no atoms found" % base)
    for st in structs:
        st.path = os.path.abspath(path)
        st.kind = kind
    return structs


# --------------------------------------------------------------------------
# Cluster building
# --------------------------------------------------------------------------

class Options(object):
    def __init__(self, edge="K", s02=1.0, rmax=5.0, radius=10.0, hydrogens=False,
                 molecule=True, drop_minor=True, distances=False):
        self.edge = edge
        self.s02 = s02
        self.rmax = rmax
        self.radius = radius
        self.hydrogens = hydrogens
        self.molecule = molecule
        self.drop_minor = drop_minor
        self.distances = distances


class CAtom(object):
    __slots__ = ("el", "label", "x", "y", "z", "d", "is_abs", "part", "ipot",
                 "tag", "key", "excluded", "n")

    def __init__(self, el, label, x, y, z, is_abs, part):
        self.el, self.label, self.x, self.y, self.z = el, label, x, y, z
        self.d = math.sqrt(x * x + y * y + z * z)
        self.is_abs, self.part = is_abs, part
        self.ipot, self.tag, self.excluded, self.n = 0, "", False, 0
        self.key = (label, round(x, 2), round(y, 2), round(z, 2))


def _bfs_points(points, els, start):
    if not points:
        return []
    cut_max = 2 * max(cov(e) for e in set(els)) + BOND_TOL
    g = _grid(points, cut_max)
    seen, stack = {start}, [start]
    min2 = MIN_BOND * MIN_BOND
    while stack:
        i = stack.pop()
        p = points[i]
        for j in _near(g, p, cut_max):
            if j in seen:
                continue
            q = points[j]
            d2 = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
            c = cov(els[i]) + cov(els[j]) + BOND_TOL
            if min2 < d2 <= c * c:
                seen.add(j)
                stack.append(j)
    return sorted(seen)


def _cluster_molecular(st, idx, opt):
    atoms = st.atoms
    p0 = atoms[idx][2]
    if opt.molecule:
        keep = _bfs_points([a[2] for a in atoms], [a[0] for a in atoms], idx)
    else:
        keep = range(len(atoms))
    out = []
    for i in keep:
        el, lab, p, part = atoms[i]
        out.append(CAtom(el, lab, p[0] - p0[0], p[1] - p0[1], p[2] - p0[2], i == idx, part))
    notes = []
    if opt.molecule and len(out) == 1:
        notes.append("Nothing is bonded to the absorber; switch to 'Everything within radius'.")
    return out, notes


def _ucell_bonds(st):
    if st._bonds is not None:
        return st._bonds
    L, uc = st.lattice, st.ucell
    carts = [f2c(u[2], L) for u in uc]
    els = [u[0] for u in uc]
    cut_max = 2 * max(cov(e) for e in set(els)) + BOND_TOL
    n = [max(1, int(math.ceil(cut_max * r))) for r in recip_lengths(L)]
    pts, owner = [], []
    for ta in range(-n[0], n[0] + 1):
        for tb in range(-n[1], n[1] + 1):
            for tc in range(-n[2], n[2] + 1):
                T = f2c((ta, tb, tc), L)
                for j, c in enumerate(carts):
                    pts.append((c[0] + T[0], c[1] + T[1], c[2] + T[2]))
                    owner.append((j, (ta, tb, tc)))
    g = _grid(pts, cut_max)
    min2 = MIN_BOND * MIN_BOND
    bonds = [[] for _ in uc]
    for i, c in enumerate(carts):
        for k in _near(g, c, cut_max):
            j, t = owner[k]
            if j == i and t == (0, 0, 0):
                continue
            q = pts[k]
            d2 = (c[0] - q[0]) ** 2 + (c[1] - q[1]) ** 2 + (c[2] - q[2]) ** 2
            cc = cov(els[i]) + cov(els[j]) + BOND_TOL
            if min2 < d2 <= cc * cc:
                bonds[i].append((j, t))
    st._bonds = bonds
    return bonds


def _cluster_periodic(st, site, opt):
    L, uc = st.lattice, st.ucell
    notes = []
    f0 = st.sites[site][2]
    origin = f2c(f0, L)
    i0, t0, best = None, (0, 0, 0), 1e9
    for i, u in enumerate(uc):
        if u[4] != site:
            continue
        d = [f0[k] - u[2][k] for k in range(3)]
        r = [int(round(v)) for v in d]
        dc = f2c([d[k] - r[k] for k in range(3)], L)
        dd = dot(dc, dc)
        if dd < best:
            i0, t0, best = i, tuple(r), dd
    out = []
    if opt.molecule:
        bonds = _ucell_bonds(st)
        limit = opt.radius + 6.0 if opt.radius > 0 else 30.0
        seen, stack, cut = {(i0, t0)}, [(i0, t0)], False
        cell_of = {i0: t0}      # a discrete molecule holds each unit-cell atom once
        network = False
        while stack:
            i, t = stack.pop()
            for j, dt in bonds[i]:
                tt = (t[0] + dt[0], t[1] + dt[1], t[2] + dt[2])
                node = (j, tt)
                if node in seen:
                    continue
                if cell_of.setdefault(j, tt) != tt:
                    network = True
                f = uc[j][2]
                c = f2c((f[0] + tt[0], f[1] + tt[1], f[2] + tt[2]), L)
                dx, dy, dz = c[0] - origin[0], c[1] - origin[1], c[2] - origin[2]
                if dx * dx + dy * dy + dz * dz > limit * limit:
                    cut = True
                    continue
                seen.add(node)
                stack.append(node)
        if network:
            notes.append("The bonds repeat through the crystal (an extended solid or polymer), so "
                         "there is no single molecule; 'Everything within radius' is usually the "
                         "better choice here.")
        if cut and opt.radius <= 0:
            notes.append("The bonded fragment is polymeric or very large; it was cut at %.0f A." % limit)
        if len(seen) == 1:
            notes.append("Nothing is bonded to the absorber; switch to 'Everything within radius'.")
        for j, t in seen:
            el, lab, f, part, s = uc[j]
            c = f2c((f[0] + t[0], f[1] + t[1], f[2] + t[2]), L)
            out.append(CAtom(el, lab, c[0] - origin[0], c[1] - origin[1], c[2] - origin[2],
                             (j, t) == (i0, t0), part))
    else:
        R = opt.radius
        if R <= 0:
            R = 8.0
            notes.append("A crystal has no edge, so radius 0 was treated as 8 A.")
        n = [int(math.ceil(R * r)) + 1 for r in recip_lengths(L)]
        rel = []
        for u in uc:
            c = f2c(u[2], L)
            rel.append((c[0] - origin[0], c[1] - origin[1], c[2] - origin[2]))
        R2 = R * R
        for ta in range(-n[0], n[0] + 1):
            for tb in range(-n[1], n[1] + 1):
                for tc in range(-n[2], n[2] + 1):
                    T = f2c((ta, tb, tc), L)
                    for j, c in enumerate(rel):
                        x, y, z = c[0] + T[0], c[1] + T[1], c[2] + T[2]
                        if x * x + y * y + z * z <= R2:
                            u = uc[j]
                            out.append(CAtom(u[0], u[1], x, y, z,
                                             j == i0 and (ta, tb, tc) == t0, u[3]))
    return out, notes


def _natural(label):
    """Sort key so H2 comes before H10."""
    return [(0, int(t), "") if t.isdigit() else (1, 0, t) for t in re.split(r"(\d+)", label or "") if t]


def _fix_representatives(used, opt, notes):
    """FEFF6 builds each unique potential around the FIRST atom listed with that ipot.
    If that atom has no bonded neighbour (e.g. a metal at the cluster edge whose ligands
    were cut off) FEFF6 crashes in the muffin-tin step. Put a bonded atom first instead,
    or drop the element if every copy is unbonded and beyond RMAX."""
    if len(used) < 2:
        return used
    pts = [(a.x, a.y, a.z) for a in used]
    size = 2 * max(cov(a.el) for a in used) + BOND_TOL
    g = _grid(pts, size)
    min2 = MIN_BOND * MIN_BOND
    memo = {}

    def bonded(i):
        if i not in memo:
            a, p, ok = used[i], pts[i], False
            for j in _near(g, p, size):
                if j == i:
                    continue
                q = pts[j]
                d2 = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
                c = cov(a.el) + cov(used[j].el) + BOND_TOL
                if min2 < d2 <= c * c:
                    ok = True
                    break
            memo[i] = ok
        return memo[i]

    order = list(range(len(used)))
    drop = set()
    for el in sorted({a.el for a in used if not a.is_abs}):
        slots = [k for k, i in enumerate(order) if used[i].el == el and not used[i].is_abs]
        if bonded(order[slots[0]]):
            continue
        good = next((k for k in slots[1:] if bonded(order[k])), None)
        if good is not None:
            order.insert(slots[0], order.pop(good))
            notes.append("Listed a bonded %s first so FEFF6 can build its potential "
                         "(the nearest %s sits at the cluster edge)." % (el, el))
        elif all(used[order[k]].d > opt.rmax for k in slots):
            drop.update(order[k] for k in slots)
            notes.append("Dropped %d %s atom(s): all are beyond RMAX with no bonded neighbours, "
                         "which crashes FEFF6." % (len(slots), el))
        else:
            notes.append("No %s atom has a bonded neighbour in the cluster; FEFF6 may fail "
                         "building its potential." % el)
    return [used[i] for i in order if i not in drop]


def _fmt(v):
    s = ("%.4f" % v).rstrip("0")
    return s + "0" if s.endswith(".") else s


def _tag(label, el):
    t = re.sub(r"[^A-Za-z0-9_]", "", label or "")[:6]
    return t or el


class Result(object):
    pass


def build(st, idx, opt, excluded=()):
    """Make the FEFF cluster and feff.inp text for absorber entry `idx`."""
    if st.periodic:
        raw, notes = _cluster_periodic(st, idx, opt)
    else:
        raw, notes = _cluster_molecular(st, idx, opt)
    notes = list(st.notes) + notes
    abs_atoms = [a for a in raw if a.is_abs]
    if not abs_atoms:
        el, lab = st.entries()[idx][0], st.entries()[idx][1]
        raw.append(CAtom(el, lab, 0.0, 0.0, 0.0, True, "."))
    elif len(abs_atoms) > 1:
        for a in abs_atoms[1:]:
            a.is_abs = False

    if opt.drop_minor:
        n0 = len(raw)
        raw = [a for a in raw if a.is_abs or str(a.part) in KEEP_PARTS]
        if len(raw) < n0:
            notes.append("Dropped %d atom(s) from minor disorder parts." % (n0 - len(raw)))
    if not opt.hydrogens:
        raw = [a for a in raw if a.is_abs or a.el != "H"]
    if opt.radius > 0:
        raw = [a for a in raw if a.is_abs or a.d <= opt.radius + 1e-6]
    raw.sort(key=lambda a: (not a.is_abs, round(a.d, 4), -zof(a.el), _natural(a.label)))

    excluded = set(excluded)
    for a in raw:
        a.excluded = (not a.is_abs) and a.key in excluded
    used = [a for a in raw if not a.excluded]
    used = _fix_representatives(used, opt, notes)
    absorber = used[0]

    pot_els = sorted({a.el for a in used if not a.is_abs}, key=lambda e: (-zof(e), e))
    ipot = {e: i + 1 for i, e in enumerate(pot_els)}
    potentials = [(0, zof(absorber.el), absorber.el)] + [(ipot[e], zof(e), e) for e in pot_els]
    for k, a in enumerate(used):
        a.ipot = 0 if a.is_abs else ipot[a.el]
        a.tag = absorber.el if a.is_abs else _tag(a.label, a.el)
        a.n = k + 1
    for a in raw:
        if a.excluded:
            a.ipot, a.tag, a.n = ipot.get(a.el, 0), _tag(a.label, a.el), 0

    if len(pot_els) > 7:
        notes.append("%d unique potentials: FEFF6 only allows ipot 1-7, so drop or merge "
                     "some elements." % len(pot_els))
    if len(used) > FEFF6_MAX_ATOMS:
        cut_r = used[FEFF6_MAX_ATOMS - 1].d
        notes.append("%d atoms: FEFF6 refuses more than %d. Lower the radius below %.1f A "
                     "(RMAX %s only needs about %s A)." % (len(used), FEFF6_MAX_ATOMS, cut_r,
                                                           _fmt(opt.rmax), _fmt(opt.rmax + 2)))
    if opt.radius > 0 and opt.radius < opt.rmax:
        notes.append("Cluster radius is smaller than RMAX, so paths beyond %s A are missing." % _fmt(opt.radius))
    pts = [(a.x, a.y, a.z) for a in used]
    g = _grid(pts, 1.0)
    clash = set()
    for i, p in enumerate(pts):
        for j in _near(g, p, 1.0):
            if j > i:
                q = pts[j]
                if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 < 0.36:
                    clash.add("%s/%s" % (used[i].tag, used[j].tag))
    if clash:
        notes.append("Atoms closer than 0.6 A (disorder?): %s. FEFF will fail on these; "
                     "exclude one of each pair." % ", ".join(sorted(clash)[:6]))

    # nearest-shell summary
    others = [a for a in used if not a.is_abs]
    shell = ""
    if others:
        d1 = others[0].d
        first = [a for a in others if a.d <= d1 + 0.35]
        groups = {}
        for a in first:
            groups.setdefault(a.el, []).append(a.d)
        bits = []
        for e in sorted(groups, key=lambda e: (-len(groups[e]), e)):
            ds = groups[e]
            if max(ds) - min(ds) < 0.0005:
                bits.append("%d %s @ %.3f A" % (len(ds), e, ds[0]))
            else:
                bits.append("%d %s @ %.3f-%.3f A" % (len(ds), e, min(ds), max(ds)))
        shell = ", ".join(bits)

    abs_label = st.label_of(idx)
    mode = "bonded molecule" if opt.molecule else "all atoms"
    radius_txt = ("r <= %s A" % _fmt(opt.radius)) if opt.radius > 0 else "no radius cut"
    title1 = ("%s centered at %s EXAFS feff.inp" % (st.name, abs_label))[:76]
    title2 = ("from %s: %s, %s, H %s" % (os.path.basename(st.path) or st.name, mode, radius_txt,
                                         "on" if opt.hydrogens else "off"))[:76]

    L = []
    L.append(" TITLE " + title1)
    L.append(" TITLE " + title2)
    L.append("")
    L.append(" HOLE   %d  %s *2nd number is S0^2" % (EDGES[opt.edge], _fmt(opt.s02)))
    L.append(" *         mphase,mpath,mfeff,mchi")
    L.append(" CONTROL   1      1     1     1")
    L.append(" PRINT     1      0     0     0")
    L.append("")
    L.append(" RMAX   %s" % _fmt(opt.rmax))
    L.append("")
    L.append(" POTENTIALS")
    L.append("  * ipot   Z      tag")
    for ip, z, e in potentials:
        L.append("     %-5d %-6d %s" % (ip, z, e))
    L.append("")
    L.append(" ATOMS                   *%d atoms" % len(used))
    L.append(" *   x          y          z     ipot tag" + ("        distance" if opt.distances else ""))
    for a in used:
        x, y, z = (round(v, 6) + 0.0 for v in (a.x, a.y, a.z))     # + 0.0 turns -0.0 into 0.0
        line = " %11.6f%11.6f%11.6f  %d  %-6s" % (x, y, z, a.ipot, a.tag)
        if opt.distances:
            line += "  %9.5f" % a.d
        L.append(line.rstrip())
    L.append("")
    L.append("END")

    res = Result()
    res.text = "\n".join(L) + "\n"
    res.atoms = raw
    res.used = used
    res.potentials = potentials
    res.notes = notes
    res.shell = shell
    res.abs_label = abs_label
    res.filename = "%s_%s_feff.inp" % (_safe(st.name), _safe(abs_label))
    return res


def _safe(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_") or "structure"


def structure_files_in(folder):
    """Every file in a folder that looks like a structure (not .inp/.txt)."""
    out = []
    for f in sorted(os.listdir(folder)):
        fp = os.path.join(folder, f)
        ext = os.path.splitext(f)[1].lower()
        if os.path.isfile(fp) and ((ext in FILE_TYPES and ext not in (".txt", ".inp")) or
                                   f.upper().startswith(("POSCAR", "CONTCAR"))):
            out.append(fp)
    return out


def unique_path(st, res, folder, taken):
    """Output path that doesn't collide with one already written in this run."""
    path = os.path.join(folder, res.filename)
    if path in taken and taken[path] != st.path:      # x.cif and x.res both loaded
        ext = os.path.splitext(st.path)[1].lstrip(".") or "2"
        path = os.path.join(folder, res.filename.replace("_feff.inp", "_%s_feff.inp" % ext))
    k = 2
    base = path
    while path in taken:                              # same label twice in one file
        path = base.replace("_feff.inp", "_%d_feff.inp" % k)
        k += 1
    res.filename = os.path.basename(path)
    taken[path] = st.path
    return path


def write_inp(res, folder):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, res.filename)
    with open(path, "w") as fh:
        fh.write(res.text)
    return path


# --------------------------------------------------------------------------
# Batch (command line) mode
# --------------------------------------------------------------------------

def run_batch(argv):
    import argparse
    ap = argparse.ArgumentParser(
        prog="feff_inp_builder --batch",
        description="Write FEFF6 feff.inp files for Artemis from structure files.")
    ap.add_argument("files", nargs="+", help="structure files (.cif .xyz .pdb .mol .sdf .mol2 .res .ins ...)")
    ap.add_argument("-a", "--absorber", help="absorber element (default: heaviest metal)")
    ap.add_argument("--site", default=None,
                    help="absorber site label, e.g. Co1 (default: first one); 'all' = one file per site")
    ap.add_argument("--edge", default="K", choices=sorted(EDGES), help="edge (default K)")
    ap.add_argument("--s02", type=float, default=1.0, help="S0^2 on the HOLE card (default 1.0)")
    ap.add_argument("--rmax", type=float, default=5.0, help="RMAX in A (default 5.0)")
    ap.add_argument("-r", "--radius", type=float, default=10.0,
                    help="keep atoms within this many A of the absorber; 0 = keep all (default 10)")
    ap.add_argument("--hydrogens", action="store_true", help="keep H atoms (default: drop them)")
    ap.add_argument("--crystal", action="store_true",
                    help="keep every atom inside the radius, not just the bonded molecule")
    ap.add_argument("--keep-disorder", action="store_true", help="keep minor disorder parts")
    ap.add_argument("--distances", action="store_true", help="add a distance column to ATOMS")
    ap.add_argument("-o", "--outdir", help="output folder (default: next to each input file)")
    args = ap.parse_args([a for a in argv if a != "--batch"])
    if args.rmax < 0.5:
        ap.error("--rmax must be at least 0.5")
    if args.s02 < 0:
        ap.error("--s02 can't be negative")
    if args.radius < 0:
        ap.error("--radius can't be negative (0 means no cut)")

    opt = Options(edge=args.edge, s02=args.s02, rmax=args.rmax, radius=args.radius,
                  hydrogens=args.hydrogens, molecule=not args.crystal,
                  drop_minor=not args.keep_disorder, distances=args.distances)
    status = 0
    import glob
    files, taken = [], {}
    for f in args.files:
        if os.path.isdir(f):          # a folder: every structure file in it
            files.extend(structure_files_in(f))
        elif any(ch in f for ch in "*?["):   # Windows cmd doesn't expand *.cif, so do it here
            hits = sorted(h for h in glob.glob(f) if os.path.isfile(h))
            if not hits:
                print("!! no files match %s" % f)
                status = 1
            files.extend(hits)
        else:
            files.append(f)
    for f in files:
        try:
            structs = load_structures(f)
        except Exception as exc:
            print("!! %s" % exc)
            status = 1
            continue
        for st in structs:
            el = args.absorber or st.default_element()
            el = element_from_symbol(el) or el
            idx = st.indices_of(el)
            if not idx:
                print("!! %s: no %s atoms" % (st.name, el))
                status = 1
                continue
            if args.site and args.site.lower() == "all":
                chosen = idx
            elif args.site:
                chosen = [i for i in idx if st.label_of(i).lower() == args.site.lower()]
                if not chosen:
                    print("!! %s: no site labelled %s" % (st.name, args.site))
                    status = 1
                    continue
            else:
                chosen = [st.default_index(el)]
            for i in chosen:
                res = build(st, i, opt)
                path = unique_path(st, res, args.outdir or os.path.dirname(st.path), taken)
                write_inp(res, os.path.dirname(path))
                print("%s -> %s  (%d atoms; %s)" % (st.name, path, len(res.used), res.shell))
                for n in res.notes:
                    print("   note: %s" % n)
    return status


# --------------------------------------------------------------------------
# Window
# --------------------------------------------------------------------------

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    import tkinter.font as tkfont
except ImportError:          # pragma: no cover
    tk = None


class App(object):
    def __init__(self, root, dnd=False):
        self.root = root
        self.structs = []
        self.state = {}
        self.result = None
        self._pending = None
        self._site_map = []
        root.title("%s %s" % (APP_NAME, VERSION))
        self.k = max(1.0, root.winfo_fpixels("1i") / 96.0)     # HiDPI scale factor
        sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
        root.geometry("%dx%d" % (min(int(1200 * self.k), sw - 60), min(int(760 * self.k), sh - 90)))
        root.minsize(960, 600)
        root.report_callback_exception = self._on_error

        self.v_el = tk.StringVar()
        self.v_site = tk.StringVar()
        self.v_edge = tk.StringVar(value="K")
        self.v_s02 = tk.StringVar(value="1.0")
        self.v_rmax = tk.StringVar(value="5.0")
        self.v_radius = tk.StringVar(value="10.0")
        self.v_mode = tk.StringVar(value="molecule")
        self.v_h = tk.BooleanVar(value=False)
        self.v_dis = tk.BooleanVar(value=True)
        self.v_dist = tk.BooleanVar(value=False)
        self.v_out = tk.StringVar(value="")
        self.v_status = tk.StringVar(value="Add .cif, .xyz, .pdb, .mol/.sdf, .mol2 or .res files to start.")

        self._build_ui(dnd)
        for v in (self.v_edge, self.v_s02, self.v_rmax, self.v_radius, self.v_mode,
                  self.v_h, self.v_dis, self.v_dist):
            v.trace_add("write", lambda *a: self.schedule())
        if dnd:
            try:
                from tkinterdnd2 import DND_FILES
                root.drop_target_register(DND_FILES)
                root.dnd_bind("<<Drop>>", lambda e: self.load(root.tk.splitlist(e.data)))
            except Exception:
                pass

    # ----- layout -----
    def _build_ui(self, dnd):
        root = self.root
        style = ttk.Style(root)
        if sys.platform.startswith("linux") and "clam" in style.theme_names():
            style.theme_use("clam")
        base = tkfont.nametofont("TkDefaultFont")
        self.f_bold = base.copy()
        self.f_bold.configure(weight="bold")
        self.f_head = base.copy()
        size = int(base.cget("size"))
        self.f_head.configure(weight="bold", size=size + 2 if size > 0 else size - 3)
        self.f_mono = tkfont.nametofont("TkFixedFont").copy()
        style.configure("Head.TLabel", font=self.f_head)
        style.configure("Note.TLabel", foreground="#9a3412")
        style.configure("Treeview", rowheight=int(base.metrics("linespace") * 1.45))

        status = ttk.Label(root, textvariable=self.v_status, anchor="w", padding=(10, 4))
        status.pack(side="bottom", fill="x")
        ttk.Separator(root).pack(side="bottom", fill="x")

        paned = ttk.PanedWindow(root, orient="horizontal")
        paned.pack(fill="both", expand=True)
        left = ttk.Frame(paned, padding=(10, 10, 6, 10))
        right = ttk.Frame(paned, padding=(6, 10, 10, 10))
        paned.add(left, weight=0)
        paned.add(right, weight=1)

        # bottom of the left column first, so the buttons never get pushed off-screen
        bf = ttk.Frame(left)
        bf.pack(side="bottom", fill="x", pady=(8, 0))
        ttk.Button(bf, text="Save this one", command=self.save_current).pack(
            side="left", fill="x", expand=True)
        ttk.Button(bf, text="Convert all", command=self.convert_all).pack(
            side="left", fill="x", expand=True, padx=(6, 0))
        of = ttk.LabelFrame(left, text=" Save to (blank = next to each input file) ", padding=(8, 4))
        of.pack(side="bottom", fill="x", pady=(8, 0))
        of.columnconfigure(0, weight=1)
        self.e_out = ttk.Entry(of, textvariable=self.v_out)
        self.e_out.grid(row=0, column=0, sticky="we")
        ttk.Button(of, text="...", width=3, command=self.pick_outdir).grid(row=0, column=1, padx=(4, 0))

        sf = ttk.LabelFrame(left, text=" FEFF settings ", padding=(8, 4))
        sf.pack(side="bottom", fill="x", pady=(8, 0))
        for c in (1, 3):
            sf.columnconfigure(c, weight=1)
        pad = dict(pady=2)

        def lab(text, r, c):
            ttk.Label(sf, text=text).grid(row=r, column=c, sticky="w", padx=(0 if c == 0 else 10, 6), **pad)

        lab("Absorber", 0, 0)
        self.cb_el = ttk.Combobox(sf, textvariable=self.v_el, state="readonly", width=5)
        self.cb_el.grid(row=0, column=1, sticky="we", **pad)
        self.cb_el.bind("<<ComboboxSelected>>", lambda e: self.on_element())
        lab("Site", 0, 2)
        self.cb_site = ttk.Combobox(sf, textvariable=self.v_site, state="readonly", width=12)
        self.cb_site.grid(row=0, column=3, sticky="we", **pad)
        self.cb_site.bind("<<ComboboxSelected>>", lambda e: self.on_site())
        lab("Edge", 1, 0)
        ttk.Combobox(sf, textvariable=self.v_edge, values=("K", "L1", "L2", "L3"),
                     state="readonly", width=5).grid(row=1, column=1, sticky="we", **pad)
        lab("S0²", 1, 2)
        ttk.Entry(sf, textvariable=self.v_s02, width=8).grid(row=1, column=3, sticky="we", **pad)
        lab("RMAX (Å)", 2, 0)
        ttk.Entry(sf, textvariable=self.v_rmax, width=6).grid(row=2, column=1, sticky="we", **pad)
        lab("Radius (Å)", 2, 2)
        ttk.Entry(sf, textvariable=self.v_radius, width=8).grid(row=2, column=3, sticky="we", **pad)
        ttk.Label(sf, text="Radius = how far from the absorber to keep atoms (0 = no cut)",
                  foreground="#6b7280").grid(row=3, column=0, columnspan=4, sticky="w")
        ttk.Radiobutton(sf, text="Bonded molecule only (like Mercury)", value="molecule",
                        variable=self.v_mode).grid(row=4, column=0, columnspan=4, sticky="w", pady=(6, 0))
        ttk.Radiobutton(sf, text="Everything in radius (packing, ions, solvent)", value="crystal",
                        variable=self.v_mode).grid(row=5, column=0, columnspan=4, sticky="w")
        ttk.Checkbutton(sf, text="Include H", variable=self.v_h).grid(
            row=6, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Checkbutton(sf, text="Distance column", variable=self.v_dist).grid(
            row=6, column=2, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Checkbutton(sf, text="Drop minor disorder parts", variable=self.v_dis).grid(
            row=7, column=0, columnspan=4, sticky="w")

        # files list takes whatever height is left
        ttk.Label(left, text="Structures", style="Head.TLabel").pack(anchor="w")
        hint = "Drag files onto this window, or use Add files" if dnd else "Each structure in each file gets a row"
        ttk.Label(left, text=hint, foreground="#6b7280").pack(anchor="w", pady=(0, 4))
        row = ttk.Frame(left)
        row.pack(side="bottom", fill="x", pady=(6, 0))
        ttk.Button(row, text="Add files...", command=self.add_files).pack(side="left")
        ttk.Button(row, text="Remove", command=self.remove_selected).pack(side="left", padx=6)
        ttk.Button(row, text="Clear", command=self.clear).pack(side="left")
        lbf = ttk.Frame(left)
        lbf.pack(fill="both", expand=True)
        self.lb = tk.Listbox(lbf, height=6, width=44, exportselection=False, activestyle="none",
                             borderwidth=1, relief="solid", highlightthickness=0)
        sb = ttk.Scrollbar(lbf, command=self.lb.yview)
        self.lb.configure(yscrollcommand=sb.set)
        self.lb.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self.lb.bind("<<ListboxSelect>>", lambda e: self.on_select())

        # right side
        self.lbl_title = ttk.Label(right, text="No structure loaded", style="Head.TLabel")
        self.lbl_title.pack(anchor="w")
        self.lbl_shell = ttk.Label(right, text="", justify="left")
        self.lbl_shell.pack(anchor="w", pady=(2, 0))
        self.lbl_notes = ttk.Label(right, text="", style="Note.TLabel", justify="left", wraplength=760)
        self.lbl_notes.pack(anchor="w", pady=(2, 6))
        right.bind("<Configure>", lambda e: self.lbl_notes.configure(wraplength=max(300, e.width - 30)))

        nb = ttk.Notebook(right)
        nb.pack(fill="both", expand=True)
        at = ttk.Frame(nb, padding=6)
        tt = ttk.Frame(nb, padding=6)
        nb.add(at, text="  Atoms  ")
        nb.add(tt, text="  feff.inp preview  ")

        tools = ttk.Frame(at)
        tools.pack(fill="x", pady=(0, 6))
        ttk.Button(tools, text="Exclude / include selected", command=self.toggle_selected).pack(side="left")
        ttk.Button(tools, text="Restore all", command=self.restore_all).pack(side="left", padx=6)
        ttk.Label(tools, text="Ctrl/Shift-click rows to pick atoms to leave out of the model",
                  foreground="#6b7280").pack(side="left", padx=8)
        tf = ttk.Frame(at)
        tf.pack(fill="both", expand=True)
        cols = ("n", "tag", "el", "ipot", "r", "x", "y", "z")
        heads = ("#", "Tag", "El", "ipot", "R (Å)", "x", "y", "z")
        widths = (50, 90, 50, 50, 90, 100, 100, 100)
        self.tree = ttk.Treeview(tf, columns=cols, show="headings", selectmode="extended")
        for c, h, w in zip(cols, heads, widths):
            self.tree.heading(c, text=h)
            self.tree.column(c, width=int(w * self.k), anchor="e" if c in ("r", "x", "y", "z", "n", "ipot") else "w",
                             stretch=c in ("x", "y", "z"))
        self.tree.tag_configure("abs", font=self.f_bold, background="#e0ecff")
        self.tree.tag_configure("excl", foreground="#9ca3af")
        ys = ttk.Scrollbar(tf, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=ys.set)
        self.tree.pack(side="left", fill="both", expand=True)
        ys.pack(side="right", fill="y")
        self.tree.bind("<Double-1>", lambda e: self.toggle_selected())
        self.tree.bind("<Delete>", lambda e: self.toggle_selected())

        self.txt = tk.Text(tt, wrap="none", font=self.f_mono, borderwidth=1, relief="solid",
                           highlightthickness=0, padx=8, pady=6)
        tys = ttk.Scrollbar(tt, orient="vertical", command=self.txt.yview)
        txs = ttk.Scrollbar(tt, orient="horizontal", command=self.txt.xview)
        self.txt.configure(yscrollcommand=tys.set, xscrollcommand=txs.set, state="disabled")
        txs.pack(side="bottom", fill="x")
        self.txt.pack(side="left", fill="both", expand=True)
        tys.pack(side="right", fill="y")

    # ----- helpers -----
    def _on_error(self, exc, val, tb):
        import traceback
        msg = "".join(traceback.format_exception(exc, val, tb))
        messagebox.showerror(APP_NAME, "Something went wrong:\n\n%s" % msg[-1500:])

    def current(self):
        sel = self.lb.curselection()
        if not sel:
            return None
        i = sel[0]
        return self.structs[i] if i < len(self.structs) else None

    def options(self):
        def f(var, name, lo=None):
            try:
                v = float(var.get())
            except ValueError:
                raise ValueError("%s must be a number" % name)
            if lo is not None and v < lo:
                raise ValueError("%s must be at least %s" % (name, lo))
            return v
        return Options(edge=self.v_edge.get(), s02=f(self.v_s02, "S0²", 0),
                       rmax=f(self.v_rmax, "RMAX", 0.5), radius=f(self.v_radius, "Cluster radius", 0),
                       hydrogens=self.v_h.get(), molecule=self.v_mode.get() == "molecule",
                       drop_minor=self.v_dis.get(), distances=self.v_dist.get())

    def _state(self, st):
        s = self.state.get(id(st))
        if s is None:
            el = st.default_element()
            s = {"el": el, "site": st.default_index(el), "excl": {}}
            self.state[id(st)] = s
        return s

    def _sites_for(self, st):
        """Concrete absorber indices for the current state (several if 'all')."""
        s = self._state(st)
        if s["site"] == "all":
            return st.indices_of(s["el"])
        return [s["site"]]

    # ----- file handling -----
    def add_files(self):
        pats = " ".join("*" + e for e in FILE_TYPES) + " POSCAR* CONTCAR*"
        paths = filedialog.askopenfilenames(
            title="Choose structure files",
            filetypes=[("Structure files", pats), ("All files", "*.*")])
        if paths:
            self.load(paths)

    def load(self, paths):
        errors, first, files = [], None, []
        for p in paths:
            p = str(p).strip()
            if os.path.isdir(p):          # a dropped folder: take every structure file in it
                files.extend(structure_files_in(p))
            elif os.path.isfile(p):
                files.append(p)
        for p in files:
            try:
                structs = load_structures(p)
            except Exception as exc:
                errors.append(str(exc))
                continue
            for st in structs:
                self.structs.append(st)
                self.lb.insert("end", "%s   [%s, %s]" % (st.name, os.path.splitext(st.path)[1].lstrip(".") or st.kind,
                                                       st.describe()))
                if first is None:
                    first = len(self.structs) - 1
        if first is not None:
            self.lb.selection_clear(0, "end")
            self.lb.selection_set(first)
            self.lb.see(first)
            self.on_select()
            self.v_status.set("Loaded %d structure(s)." % len(self.structs))
        if errors:
            messagebox.showwarning(APP_NAME, "Some files could not be read:\n\n" + "\n".join(errors))

    def remove_selected(self):
        sel = self.lb.curselection()
        if not sel:
            return
        i = sel[0]
        st = self.structs.pop(i)
        self.state.pop(id(st), None)
        self.lb.delete(i)
        if self.structs:
            j = min(i, len(self.structs) - 1)
            self.lb.selection_set(j)
        self.on_select()

    def clear(self):
        self.structs, self.state = [], {}
        self.lb.delete(0, "end")
        self.on_select()

    def pick_outdir(self):
        d = filedialog.askdirectory(title="Save .inp files to")
        if d:
            self.v_out.set(d)

    # ----- selection -----
    def on_select(self):
        st = self.current()
        if st is None:
            self.cb_el.configure(values=())
            self.cb_site.configure(values=())
            self.v_el.set("")
            self.v_site.set("")
            self.show(None)
            return
        s = self._state(st)
        self.cb_el.configure(values=st.elements())
        self.v_el.set(s["el"] or "")
        self._fill_sites(st)
        self.refresh()

    def _fill_sites(self, st):
        s = self._state(st)
        idx = st.indices_of(s["el"])
        seen, names, mapping = {}, [], []
        for i in idx:
            lab = st.label_of(i)
            seen[lab] = seen.get(lab, 0) + 1
            names.append(lab if seen[lab] == 1 else "%s (#%d)" % (lab, seen[lab]))
            mapping.append(i)
        if len(idx) > 1:
            names.append("All %d %s sites (one file each)" % (len(idx), s["el"]))
            mapping.append("all")
        self._site_map = mapping
        self.cb_site.configure(values=names)
        if s["site"] in mapping:
            self.v_site.set(names[mapping.index(s["site"])])
        elif names:
            s["site"] = mapping[0]
            self.v_site.set(names[0])

    def on_element(self):
        st = self.current()
        if st is None:
            return
        s = self._state(st)
        s["el"] = self.v_el.get()
        s["site"] = st.default_index(s["el"])
        self._fill_sites(st)
        self.refresh()

    def on_site(self):
        st = self.current()
        if st is None:
            return
        k = self.cb_site.current()
        if 0 <= k < len(self._site_map):
            self._state(st)["site"] = self._site_map[k]
        self.refresh()

    # ----- preview -----
    def schedule(self):
        if self._pending is not None:
            self.root.after_cancel(self._pending)
        self._pending = self.root.after(200, self.refresh)

    def refresh(self):
        self._pending = None
        st = self.current()
        if st is None:
            self.show(None)
            return
        try:
            opt = self.options()
        except ValueError as exc:
            self.lbl_notes.configure(text=str(exc))
            return
        s = self._state(st)
        site = self._sites_for(st)[0]
        res = build(st, site, opt, s["excl"].get(site, ()))
        self.result = (st, site, res)
        self.show(res, st)

    def show(self, res, st=None):
        self.tree.delete(*self.tree.get_children())
        self.txt.configure(state="normal")
        self.txt.delete("1.0", "end")
        if res is None:
            self.lbl_title.configure(text="No structure loaded")
            self.lbl_shell.configure(text="")
            self.lbl_notes.configure(text="")
            self.txt.configure(state="disabled")
            self.result = None
            return
        multi = self._state(st)["site"] == "all"
        nexcl = sum(1 for a in res.atoms if a.excluded)
        self.lbl_title.configure(text="%s  —  absorber %s%s" % (
            st.name, res.abs_label, "  (previewing first of all sites)" if multi else ""))
        info = "%d atoms, %d potentials" % (len(res.used), len(res.potentials))
        if nexcl:
            info += ", %d excluded" % nexcl
        if res.shell:
            info += "   |   nearest: " + res.shell
        self.lbl_shell.configure(text=info)
        self.lbl_notes.configure(text="\n".join("• " + n for n in res.notes))
        for k, a in enumerate(res.atoms):
            tags = ("abs",) if a.is_abs else (("excl",) if a.excluded else ())
            self.tree.insert("", "end", iid=str(k), tags=tags, values=(
                a.n if a.n else "–", a.tag, a.el, a.ipot if not a.excluded else "", "%.4f" % a.d,
                "%.5f" % a.x, "%.5f" % a.y, "%.5f" % a.z))
        self.txt.insert("1.0", res.text)
        self.txt.configure(state="disabled")

    def toggle_selected(self):
        if not self.result:
            return
        st, site, res = self.result
        s = self._state(st)
        if s["site"] == "all":
            messagebox.showinfo(APP_NAME, "Pick a single site to exclude atoms.")
            return
        excl = s["excl"].setdefault(site, set())
        for iid in self.tree.selection():
            a = res.atoms[int(iid)]
            if a.is_abs:
                continue
            if a.key in excl:
                excl.discard(a.key)
            else:
                excl.add(a.key)
        self.refresh()

    def restore_all(self):
        if not self.result:
            return
        st, site, res = self.result
        self._state(st)["excl"].pop(site, None)
        self.refresh()

    # ----- saving -----
    def _outdir(self, st):
        return self.v_out.get().strip() or os.path.dirname(st.path)

    def _write_jobs(self, structs):
        opt = self.options()
        jobs, taken = [], {}
        for st in structs:
            s = self._state(st)
            for site in self._sites_for(st):
                if site is None:
                    continue
                res = build(st, site, opt, s["excl"].get(site, ()))
                jobs.append((st, res, unique_path(st, res, self._outdir(st), taken)))
        exist = [p for _, _, p in jobs if os.path.exists(p)]
        if exist and not messagebox.askyesno(
                APP_NAME, "%d file(s) already exist, e.g.\n%s\n\nOverwrite?" % (len(exist), exist[0])):
            return []
        written = []
        for st, res, p in jobs:
            written.append(write_inp(res, os.path.dirname(p)))
        return written

    def save_current(self):
        st = self.current()
        if st is None:
            return
        try:
            written = self._write_jobs([st])
        except ValueError as exc:
            messagebox.showerror(APP_NAME, str(exc))
            return
        if written:
            self.v_status.set("Saved " + ("  |  ".join(written) if len(written) < 3 else
                                          "%d files to %s" % (len(written), os.path.dirname(written[0]))))

    def convert_all(self):
        if not self.structs:
            return
        try:
            written = self._write_jobs(self.structs)
        except ValueError as exc:
            messagebox.showerror(APP_NAME, str(exc))
            return
        if written:
            folders = sorted({os.path.dirname(p) for p in written})
            self.v_status.set("Wrote %d .inp file(s)." % len(written))
            messagebox.showinfo(APP_NAME, "Wrote %d .inp file(s) to:\n\n%s" % (len(written), "\n".join(folders)))


def run_gui(paths=()):
    if tk is None:
        sys.stderr.write("tkinter is not available in this Python. Use --batch, or add tkinter: the "
                         "python.org installers include it; on Linux install python3-tk; with Homebrew "
                         "install python-tk.\n")
        return 1
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    dnd = False
    try:
        from tkinterdnd2 import TkinterDnD
        root = TkinterDnD.Tk()
        dnd = True
    except Exception:
        root = tk.Tk()
    app = App(root, dnd=dnd)
    if paths:
        root.after(50, lambda: app.load(paths))
    root.mainloop()
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--batch" in argv:
        return run_batch(argv)
    if argv and argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    return run_gui([a for a in argv if not a.startswith("-")])


if __name__ == "__main__":
    sys.exit(main())
