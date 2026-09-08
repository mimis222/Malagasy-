#!/usr/bin/env python3
"""Genere un motif topographique (lignes de contour) en SVG, palette ALTO.
Bruit de valeur fractal + marching squares. Sans dependance externe."""
import random
import math

random.seed(7)

# ---- Champ de bruit de valeur fractal ----
def make_lattice(nx, ny):
    return [[random.random() for _ in range(nx + 1)] for _ in range(ny + 1)]

def fade(t):
    return t * t * t * (t * (t * 6 - 15) + 10)  # smootherstep

def sample_octave(lat, nx, ny, u, v):
    x = u * nx
    y = v * ny
    xi = min(int(x), nx - 1)
    yi = min(int(y), ny - 1)
    tx = fade(x - xi)
    ty = fade(y - yi)
    a = lat[yi][xi]
    b = lat[yi][xi + 1]
    c = lat[yi + 1][xi]
    d = lat[yi + 1][xi + 1]
    top = a + (b - a) * tx
    bot = c + (d - c) * tx
    return top + (bot - top) * ty

# octaves : basses frequences = grandes ondulations souples
OCTAVES = [
    (2, 3, 1.0),
    (4, 6, 0.55),
    (8, 11, 0.28),
    (16, 22, 0.14),
    (32, 44, 0.07),
]
lattices = [(make_lattice(nx, ny), nx, ny, amp) for (nx, ny, amp) in OCTAVES]
amp_total = sum(a for *_, a in lattices)

GW, GH = 176, 236  # resolution du champ echantillonne
field = [[0.0] * GW for _ in range(GH)]
mn, mx = 1e9, -1e9
for j in range(GH):
    v = j / (GH - 1)
    for i in range(GW):
        u = i / (GW - 1)
        s = 0.0
        for lat, nx, ny, amp in lattices:
            s += sample_octave(lat, nx, ny, u, v) * amp
        s /= amp_total
        field[j][i] = s
        if s < mn: mn = s
        if s > mx: mx = s

# normalise 0..1
for j in range(GH):
    for i in range(GW):
        field[j][i] = (field[j][i] - mn) / (mx - mn)

# ---- Marching squares ----
def interp(p1, p2, v1, v2, level):
    if abs(v2 - v1) < 1e-9:
        return p1
    t = (level - v1) / (v2 - v1)
    return (p1[0] + (p2[0] - p1[0]) * t, p1[1] + (p2[1] - p1[1]) * t)

def contour_segments(level):
    segs = []
    for j in range(GH - 1):
        for i in range(GW - 1):
            tl = field[j][i]
            tr = field[j][i + 1]
            br = field[j + 1][i + 1]
            bl = field[j + 1][i]
            idx = 0
            if tl > level: idx |= 8
            if tr > level: idx |= 4
            if br > level: idx |= 2
            if bl > level: idx |= 1
            if idx == 0 or idx == 15:
                continue
            ptl = (i, j); ptr = (i + 1, j); pbr = (i + 1, j + 1); pbl = (i, j + 1)
            top = lambda: interp(ptl, ptr, tl, tr, level)
            right = lambda: interp(ptr, pbr, tr, br, level)
            bottom = lambda: interp(pbl, pbr, bl, br, level)
            left = lambda: interp(ptl, pbl, tl, bl, level)
            e = {
                1: [(left, bottom)],
                2: [(bottom, right)],
                3: [(left, right)],
                4: [(top, right)],
                5: [(left, top), (bottom, right)],
                6: [(top, bottom)],
                7: [(left, top)],
                8: [(left, top)],
                9: [(top, bottom)],
                10: [(left, bottom), (top, right)],
                11: [(top, right)],
                12: [(left, right)],
                13: [(bottom, right)],
                14: [(left, bottom)],
            }[idx]
            for a, b in e:
                segs.append((a(), b()))
    return segs

def key(p):
    return (round(p[0], 4), round(p[1], 4))

def join(segs):
    """Assemble les segments en polylignes continues."""
    from collections import defaultdict
    adj = defaultdict(list)
    for a, b in segs:
        adj[key(a)].append((key(b), b))
        adj[key(b)].append((key(a), a))
    visited = set()
    lines = []
    seg_set = set()
    for a, b in segs:
        seg_set.add((key(a), key(b)))
    used = set()

    def take(ka, kb):
        used.add((ka, kb)); used.add((kb, ka))

    # index rapide vers coords
    coord = {}
    for a, b in segs:
        coord[key(a)] = a; coord[key(b)] = b

    for a, b in segs:
        ka, kb = key(a), key(b)
        if (ka, kb) in used:
            continue
        take(ka, kb)
        line = [a, b]
        # etendre vers l'avant
        cur = kb
        while True:
            nxt = None
            for nk, npt in adj[cur]:
                if (cur, nk) not in used:
                    nxt = (nk, npt); break
            if not nxt:
                break
            take(cur, nxt[0])
            line.append(nxt[1])
            cur = nxt[0]
        # etendre vers l'arriere
        cur = ka
        pts_front = []
        while True:
            nxt = None
            for nk, npt in adj[cur]:
                if (cur, nk) not in used:
                    nxt = (nk, npt); break
            if not nxt:
                break
            take(cur, nxt[0])
            pts_front.append(nxt[1])
            cur = nxt[0]
        line = list(reversed(pts_front)) + line
        lines.append(line)
    return lines

def simplify(pts, eps):
    """Douglas-Peucker."""
    if len(pts) < 3:
        return pts
    dmax, idx = 0.0, 0
    a, b = pts[0], pts[-1]
    dx, dy = b[0] - a[0], b[1] - a[1]
    nrm = math.hypot(dx, dy) or 1e-9
    for i in range(1, len(pts) - 1):
        px, py = pts[i]
        d = abs((px - a[0]) * dy - (py - a[1]) * dx) / nrm
        if d > dmax:
            dmax, idx = d, i
    if dmax > eps:
        left = simplify(pts[:idx + 1], eps)
        right = simplify(pts[idx:], eps)
        return left[:-1] + right
    return [a, b]

# ---- Construction du SVG ----
LEVELS = 42
paths = []
for k in range(1, LEVELS):
    level = k / LEVELS
    segs = contour_segments(level)
    if not segs:
        continue
    for line in join(segs):
        line = simplify(line, 0.09)
        if len(line) < 2:
            continue
        d = "M" + " ".join(f"{x:.1f},{y:.1f}" for x, y in line)
        # variation d'epaisseur/opacite selon l'altitude -> profondeur
        t = level
        paths.append((t, d))

W, H = GW - 1, GH - 1
out = []
out.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" preserveAspectRatio="xMidYMid slice" width="{W}" height="{H}">')
out.append('<defs>')
out.append(f'''<linearGradient id="bg" x1="0" y1="0" x2="1" y2="1">
  <stop offset="0" stop-color="#0a2138"/>
  <stop offset="0.5" stop-color="#071726"/>
  <stop offset="1" stop-color="#04101b"/>
</linearGradient>
<radialGradient id="glow" cx="0.32" cy="0.42" r="0.75">
  <stop offset="0" stop-color="#c38a2e" stop-opacity="0.16"/>
  <stop offset="1" stop-color="#c38a2e" stop-opacity="0"/>
</radialGradient>''')
out.append('</defs>')
out.append(f'<rect width="{W}" height="{H}" fill="url(#bg)"/>')
out.append(f'<rect width="{W}" height="{H}" fill="url(#glow)"/>')
out.append('<g fill="none" stroke-linecap="round" stroke-linejoin="round">')
for t, d in paths:
    # lignes plus claires et plus fines en altitude, plus chaudes en bas
    if t > 0.72:
        col = "#f6d792"; w = 0.36; op = 1.0
    elif t > 0.5:
        col = "#e9c473"; w = 0.32; op = 0.92
    elif t > 0.3:
        col = "#cf9435"; w = 0.3; op = 0.82
    else:
        col = "#a8752a"; w = 0.28; op = 0.72
    out.append(f'<path d="{d}" stroke="{col}" stroke-width="{w}" opacity="{op}"/>')
out.append('</g>')
out.append('</svg>')

svg = "\n".join(out)
with open("topo.svg", "w") as f:
    f.write(svg)
print("topo.svg ecrit :", len(svg), "octets,", len(paths), "chemins")
