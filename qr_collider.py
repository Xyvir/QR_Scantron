#!/usr/bin/env python3
"""
qr_collider.py -- search for "pen-flip" QR collision pairs.

Goal
----
Find two URLs in a family (e.g. lithic.uk/e and lithic.uk/f) whose QR symbols
are close enough that a human can redirect the destination by filling in a
handful of white modules with a pen.  Scan before -> A, scan after -> B.

How it actually works (this is the important part)
--------------------------------------------------
You never have to turn symbol A into symbol B exactly.  A QR symbol carries
Reed-Solomon parity, so a scanner repairs up to `t` damaged codewords
(t = parity_codewords // 2) before it even looks at the payload.  That gives
three free levers, all of which this tool exploits:

  * dark->light modules never need erasing.  Leave them dark; the scanner of the
    "after" image repairs them.  A pen that only adds black is enough, as long
    as those modules fit inside the repair budget.
  * light->dark work can be skipped entirely in up to `t` codewords.  A
    codeword's 8 modules occupy one small 2x4 block, so skipping the heaviest
    codewords deletes whole clumps of pen work for one unit of budget each.
  * --preblank splits the two budgets.  Instead of printing symbol A, print it
    with the unreachable modules already erased, so they are not errors at all
    for the after scan.  It still scans as A (the damage is inside its repair
    budget), and the after scan can now spend its whole budget on skipping fill
    work.  --margin buys headroom against print noise.

The mask is a lever too, and not for the reason you might expect
---------------------------------------------------------------
Two symbols sharing a mask are both XORed with the same pattern, so masking
cancels out of the pairwise difference: the SET of differing modules cannot be
moved by any mask, and the raw distance is mask-invariant (measured identical
for masks 0/1/2).  Masks therefore cannot slide a differing module next to
another one.

They do decide which differing modules a pen can reach.  A printed module is
fillable iff it is light, i.e. P_m(p) == C_A(p), and filling sets the *unmasked*
value to 1 XOR P_m(p) -- so where P_m(p) == 1 a light->dark stroke pushes the
decoded bit DOWN, repairing a difference that mask 0 could only repair by
erasing.  Since fillable and unreachable are complementary, each mask gives a
different partition of the same differing set: one mask may need fewer fills but
strand them in many codewords, another may need more fills but concentrate the
stranded ones where ECC can absorb them.  So the mask trades fill count against
clumping, and it changes the feasibility of every pair -- searching all 8
broadens the pool of usable pairs roughly eightfold.

The printed sheet is self-sufficient
------------------------------------
fillmap.png carries the answer: a light outline on every cell to fill, with no
numbers and no separate reference to consult.  The stroke is drawn
INSIDE the module it marks, so it never touches a neighbour and it shows the
true size of the pixel to ink; outlines rather than solid highlights also keep
each module's centre white, and the centre is exactly where a decoder samples,
so the highlighted sheet still scans as the "before" url -- here the pretty
option and the robust option are the same one.  Cells that share an edge fuse
into a single region (no line is drawn across the shared edge, so their white
interiors join up), and cells that meet only at a corner keep a square corner
there instead of a rounded one, so the two outlines meet and a diagonal pair
reads as one connected region.

The default outline colour is a light magenta (--highlight #ffd9ff), and that
is not a preference but a printing fact: it is a tint of exactly ONE ink.  A
light yellow such as #ffd633 is yellow plus magenta, so the printer has to mix
two inks to build it, and when that mix is off the sheet comes out pink or
purple instead of yellow.  A yellow-dependent tint also inherits the yellow
cartridge -- a two-ink mix turns pink the moment yellow runs low -- and yellow
has the least contrast on white paper, so it fades first on recycled stock or
through a copy.  A tint of one ink has nothing to mix, so its hue cannot drift.
--highlight takes any colour; swatches.py prints the printer-safe candidates.

Everything runs on the module matrix.  No images, no PNGs, no decoder in the
search loop: a symbol is one Python int, so the module distance between two
candidates is a popcount.  SVG/ASCII/PNG output exists only for the winner, and
a real decoder (OpenCV) is only consulted for the final pair, with --verify.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import random
import sys
from collections import deque

import segno
import segno.encoder as E
from segno import consts as C

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

ERROR_LEVEL = {"L": C.ERROR_LEVEL_L, "M": C.ERROR_LEVEL_M,
               "Q": C.ERROR_LEVEL_Q, "H": C.ERROR_LEVEL_H}

# ---------------------------------------------------------------------------
# symbol geometry: which module each message bit lives in
# ---------------------------------------------------------------------------
_GEOM = {}


def geometry(version):
    """(size, free, order) for a QR version.

    `order` is the ISO/IEC 18004 7.7.3 placement order: two-module-wide columns
    zig-zagging right to left, skipping function patterns.  Bit i of the final
    message lands on order[i], so codeword c owns order[8c : 8c+8] -- on screen
    one small 2x4 block, which is why skipping a codeword removes a clump.
    """
    if version in _GEOM:
        return _GEOM[version]
    size = E.calc_matrix_size(version)
    m = E.make_matrix(size, size)
    E.add_finder_patterns(m, size, size)
    E.add_alignment_patterns(m, size, size)
    free = [[bool(m[r][c] & 0x2) for c in range(size)] for r in range(size)]
    order, upward, right = [], True, size - 1
    while right > 0:
        if right == 6:                      # timing column
            right -= 1
        for k in range(size):
            row = (size - 1 - k) if upward else k
            for col in (right, right - 1):
                if free[row][col]:
                    order.append((row, col))
        upward = not upward
        right -= 2
    _GEOM[version] = (size, free, order)
    return _GEOM[version]


def raw_message_bits(url, version, error):
    """The final codeword stream (data + Reed-Solomon parity), as bits.

    Assembled from segno's own primitives so that comparing it against what we
    read out of a finished matrix is an independent check on `geometry()`,
    not a restatement of it.
    """
    segs = E.prepare_data(url, None, None)
    buff = E.Buffer()
    E.write_segment(buff, segs[0], None, E.version_range(version), False)
    cap = C.SYMBOL_CAPACITY[version][error]
    E.write_terminator(buff, cap, None, len(buff))
    E.write_padding_bits(buff, version, len(buff))
    E.write_pad_codewords(buff, version, cap, len(buff))
    return list(E.make_final_message(version, error, buff))


class Codec:
    """One (version, error level) plus all 8 mask patterns, as packed ints.

    Bit i of every int below is message bit i (placement slot i), so codeword c
    is bits 8c..8c+8 and a popcount is a module distance.
    """

    def __init__(self, version, error_name):
        self.version = version
        self.error_name = error_name
        self.error = ERROR_LEVEL[error_name]
        self.size, self.free, self.order = geometry(version)
        info = C.ECC[version][self.error]
        if any(e.num_blocks != 1 for e in info):
            raise SystemExit(
                f"version {version}-{error_name} uses multiple Reed-Solomon "
                "blocks; this tool models single-block symbols only (v1/v2 and "
                "other low versions). Use a shorter payload or a lower --ecc.")
        self.block = info[0]
        self.parity = self.block.num_total - self.block.num_data
        self.n_codewords = self.block.num_total
        self.correctable = self.parity // 2
        self.min_symbol_distance = self.parity + 1
        self.n_slots = len(self.order)
        self.remainder = self.n_slots - 8 * self.n_codewords
        self.n_bits = 8 * self.n_codewords
        self.bit_mask = (1 << self.n_bits) - 1
        # P[m]: which data modules mask m toggles, as a packed int.  Only the
        # codeword slots are kept: the remainder/filler modules carry the mask
        # pattern too but are identically zero in every payload, so they never
        # differ and never contribute to a distance.
        self.mask_bits = [
            sum((1 if fn(r, c) else 0) << i
                for i, (r, c) in enumerate(self.order) if i < self.n_bits)
            for fn in E.get_data_mask_functions(False)]

    def matrix(self, url, mask=0):
        q = segno.make_qr(url, version=self.version, error=self.error_name,
                          mask=mask, boost_error=False)
        return [[int(v) for v in row] for row in q.matrix]

    def printed_bits(self, rows):
        """Masked data-region bits, i.e. the modules as they appear on paper."""
        v = 0
        for i, (r, c) in enumerate(self.order):
            if rows[r][c]:
                v |= 1 << i
        return v & self.bit_mask

    def filler_bits(self, rows):
        """Remainder modules past the codeword stream (constant for all urls)."""
        v = 0
        for i in range(self.n_bits, self.n_slots):
            r, c = self.order[i]
            v |= int(rows[r][c]) << (i - self.n_bits)
        return v

    def codeword_bits(self, url_or_rows, mask=0):
        """Canonical (unmasked) codeword array: the mask is XORed back out."""
        rows = url_or_rows if isinstance(url_or_rows, list) \
            else self.matrix(url_or_rows, mask)
        return (self.printed_bits(rows) ^ self.mask_bits[mask]) & self.bit_mask

    def cell(self, bit):
        return self.order[bit]


# ---------------------------------------------------------------------------
# image-free metrics
# ---------------------------------------------------------------------------
def blobs(cells, diag=True):
    """Connected components of a cell set.

    diag=False counts only edge-sharing neighbours, which is what one
    uninterrupted pen stroke can cover.  diag=True also merges cells that
    merely touch at a corner.
    """
    steps = ([(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0),
              (1, 1)] if diag else [(-1, 0), (1, 0), (0, -1), (0, 1)])
    cells = set(cells)
    seen, out = set(), []
    for start in cells:
        if start in seen:
            continue
        q, comp = deque([start]), []
        seen.add(start)
        while q:
            cr, cc = q.popleft()
            comp.append((cr, cc))
            for dr, dc in steps:
                nxt = (cr + dr, cc + dc)
                if nxt in cells and nxt not in seen:
                    seen.add(nxt)
                    q.append(nxt)
        out.append(comp)
    return out


def clump_metrics(cells):
    """How pen-friendly the fill set is, purely geometrically.

    regions         groups one uninterrupted stroke can cover (edge-sharing)
    blobs8          groups when corner contacts also count as joined
    regions_eff     blobs8, but each corner-only merge counts as HALF a merge:
                    joining two cells along a shared edge is worth much more to
                    a hand than clipping a corner, so a diagonal pair is scored
                    as 1.5 regions rather than a full 1
    diag_links      groups held together only by a corner touch
    orth_pairs      edge-sharing neighbour pairs (the good kind)
    biggest_region  size of the largest single-stroke group
    """
    cells = set(cells)
    edge = blobs(cells, diag=False)
    corner = blobs(cells, diag=True)
    merges = len(edge) - len(corner)          # merges caused by corner contact
    return {
        "regions": len(edge),
        "blobs8": len(corner),
        "regions_eff": len(corner) + 0.5 * merges,
        "diag_links": merges,
        "orth_pairs": sum(1 for (r, c) in cells
                          if (r + 1, c) in cells or (r, c + 1) in cells),
        "biggest_region": max((len(g) for g in edge), default=0),
    }


def popcount(a):
    if hasattr(np, "bitwise_count"):
        return np.bitwise_count(a)
    lut = np.array([bin(b).count("1") for b in range(256)], dtype=np.uint8)
    return lut[a]


def _skip_saving(up_per_cw, mask_free, budget):
    """Fill modules deleted by spending `budget` codewords.

    `mask_free` is non-zero only where a codeword may be abandoned (in
    pristine-print mode a codeword already broken by an unreachable module is
    abandoned for free, being an error either way).  `budget` is per-pair.
    """
    pool = np.where(mask_free, up_per_cw, 0)
    desc = -np.sort(-pool, axis=1)
    csum = np.cumsum(desc, axis=1)
    idx = np.clip(budget, 1, up_per_cw.shape[1]) - 1
    rows = np.arange(len(csum))
    return np.where(budget > 0, csum[rows, idx], 0)


# ---------------------------------------------------------------------------
# the search
# ---------------------------------------------------------------------------
def scan(codec, urls, canon, cfg):
    """One vectorised pass per mask over every pair, sharing one encode per url.

    For each (pair, mask) exactly, and without building any image: how many
    modules a pen can reach, which codewords the unreachable ones strand, and
    how many fills survive after spending repair budget on the heaviest
    codewords.
    """
    n = len(canon)
    ncw = codec.n_codewords
    t = codec.correctable - cfg.margin
    masks = cfg.masks
    b = np.empty((n, ncw), dtype=np.uint8)
    for i, k in enumerate(canon):
        for j in range(ncw):
            b[i, j] = (k >> (8 * j)) & 0xFF
    pm = np.array([[(m >> (8 * j)) & 0xFF for j in range(ncw)]
                   for m in codec.mask_bits], dtype=np.uint8)

    passing, reach_pairs = [], set()
    stat = dict(pairs=0, reachable_combos=0, best_raw=None, best_raw_pair=None)
    described = False
    for m in masks:
        printed = b ^ pm[m]                 # this is what a scanner sees
        for i in range(n - 1):
            # differing modules are mask-invariant; the mask only decides which
            # of them are light on paper (reachable) and which are dark
            # (stranded).  Printed_B is forced to the opposite value here, so
            # "light in A and differing" already implies "dark in B".
            diff = b[i] ^ b[i + 1:]
            up_per = popcount(diff & ~printed[i]).astype(np.int32)
            st_per = popcount(diff & printed[i]).astype(np.int32)
            raw = popcount(diff).sum(axis=1).astype(np.int32)
            up_bits = up_per.sum(1)
            k_down = (st_per > 0).sum(1)

            budget = (np.full(len(k_down), t, dtype=np.int32) if cfg.preblank
                      else np.clip(t - k_down, 0, None))
            if cfg.slack:
                free = np.ones_like(up_per, dtype=bool) if cfg.preblank \
                    else (st_per == 0)
                saving = _skip_saving(up_per, free, budget)
                if not cfg.preblank:
                    saving = saving + (up_per * (st_per > 0)).sum(1)
            else:
                saving = np.zeros(len(up_per), dtype=np.int32)
            fills = up_bits - saving

            ok = k_down <= t
            if cfg.max_raw is not None:
                ok &= raw <= cfg.max_raw

            stat["reachable_combos"] += int(ok.sum())
            for j in np.nonzero(ok)[0]:
                passing.append((int(fills[j]), i, i + 1 + int(j), m))
                reach_pairs.add((i, i + 1 + int(j)))
            if not described:
                if stat["best_raw"] is None or raw.min() < stat["best_raw"]:
                    stat["best_raw"] = int(raw.min())
                    stat["best_raw_pair"] = (
                        urls[i], urls[i + 1 + int(raw.argmin())])
        if not described:
            stat["pairs"] = n * (n - 1) // 2
            described = True

    passing.sort()
    scored = []
    for _, i, j, m in passing[:cfg.refine]:
        detail = analyse(codec, urls, canon, i, j, cfg, m)
        if detail:
            scored.append(detail)
    scored.sort(key=lambda d: (d["score"], d["fills"]))
    stat["reachable_pairs"] = len(reach_pairs)
    stat["meets_target"] = sum(1 for f, _, _, _ in passing
                               if f <= cfg.max_fills)
    return (scored[0] if scored else None), scored, stat


def analyse(codec, urls, canon, i, j, cfg, mask=0):
    """Full detail for one (pair, mask), straight off the packed bit ints."""
    t = codec.correctable - cfg.margin
    ca, cb = canon[i], canon[j]
    printed = ca ^ codec.mask_bits[mask]          # what the before symbol shows
    up_cw, st_cw = {}, {}
    for c in range(codec.n_codewords):
        a = (ca >> (8 * c)) & 0xFF
        b_ = (cb >> (8 * c)) & 0xFF
        p = (printed >> (8 * c)) & 0xFF
        if a == b_:
            continue
        if ~p & (a ^ b_) & 0xFF:                  # light on paper -> pen work
            up_cw[c] = [8 * c + bit for bit in range(8)
                        if (a ^ b_) >> bit & 1 and not (p >> bit) & 1]
        if p & (a ^ b_) & 0xFF:                   # dark on paper -> unreachable
            st_cw[c] = [8 * c + bit for bit in range(8)
                        if (a ^ b_) >> bit & 1 and (p >> bit) & 1]
    k_down = len(st_cw)
    if k_down > t:
        return None

    skipped = set()
    if cfg.slack:
        if cfg.preblank:
            budget, pool = t, dict(up_cw)
        else:
            budget = t - k_down
            pool = {c: v for c, v in up_cw.items() if c not in st_cw}
            skipped.update(c for c in up_cw if c in st_cw)   # free: already broken
        for c in sorted(pool, key=lambda c: -len(pool[c]))[:budget]:
            skipped.add(c)
    fills = [bit for c, bits in up_cw.items() if c not in skipped for bit in bits]
    stranded = [bit for bits in st_cw.values() for bit in bits]

    cells = [codec.cell(bit) for bit in fills]
    m = clump_metrics(cells)
    return {
        # Cost is fills plus penalty per extra fill REGION.  A corner-only join
        # counts as half a merge (regions_eff), and each shared edge earns its
        # own bonus, so a clump a hand can shade in one movement outranks the
        # same number of fills scattered around, and outranks a diagonal chain
        # that would need separate dabs.  Ordering is orthogonal < diagonal <
        # isolated, with orthogonal contact worth several isolated regions.
        "score": len(fills) + cfg.blob_weight * (m["regions_eff"] - 1)
        - cfg.orth_weight * m["orth_pairs"],
        "before": urls[i], "after": urls[j], "mask": mask,
        "up": sum(len(v) for v in up_cw.values()),
        "down": sum(len(v) for v in st_cw.values()),
        "raw": sum(len(v) for v in up_cw.values())
        + sum(len(v) for v in st_cw.values()),
        "k_down": k_down, "k_up": len(up_cw),
        "fills": len(fills),
        **m,
        "skipped_codewords": sorted(skipped),
        "fill_bits": sorted(fills),
        "fill_cells": [list(c) for c in sorted(cells)],
        "stranded_cells": [list(codec.cell(b)) for b in sorted(stranded)],
        "before_headroom": t - k_down,
        "after_headroom": t - len(skipped),
        "matrix_before": codec.matrix(urls[i], mask),
        "matrix_after": codec.matrix(urls[j], mask),
    }


# ---------------------------------------------------------------------------
# candidates / codec setup
# ---------------------------------------------------------------------------
def build_candidates(args):
    if args.urls:
        with open(args.urls, encoding="utf-8") as fh:
            urls = [ln.strip() for ln in fh
                    if ln.strip() and not ln.startswith("#")]
    else:
        urls = [args.base + "".join(p) for p in
                itertools.product(args.suffix_chars, repeat=args.suffix_len)]
    urls = list(dict.fromkeys(urls))
    if len(urls) > args.max_candidates:
        random.Random(args.seed).shuffle(urls)
        urls = urls[:args.max_candidates]
        print(f"note: space truncated to {len(urls)} candidates "
              f"(raise --max-candidates for a wider search)")
    return urls


def pick_version(urls, error_name):
    longest = max(urls, key=len)
    for v in range(1, 41):
        try:
            segno.make_qr(longest, version=v, error=error_name,
                          boost_error=False)
            return v
        except E.DataOverflowError:
            continue
    raise SystemExit("payload too large for any QR version")


# ---------------------------------------------------------------------------
# output for the winner only
# ---------------------------------------------------------------------------
def printed_matrix(res, cfg):
    """What actually gets printed: the before symbol, unreachable bits removed."""
    p = [row[:] for row in res["matrix_before"]]
    if cfg.preblank:
        for r, c in res["stranded_cells"]:
            p[r][c] = 0
    return p


def final_matrix(res, cfg):
    f = printed_matrix(res, cfg)
    for r, c in res["fill_cells"]:
        f[r][c] = 1
    return f


def ascii_map(codec, res, cfg):
    fills = {tuple(c) for c in res["fill_cells"]}
    stranded = {tuple(c) for c in res["stranded_cells"]}
    skipped = set(res["skipped_codewords"])
    lines = ["    " + "".join("%x" % (i % 16) for i in range(codec.size))]
    for r in range(codec.size):
        row = []
        for c in range(codec.size):
            if (r, c) in fills:
                row.append("X")
            elif (r, c) in stranded:
                row.append("o" if cfg.preblank else "D")
            elif res["matrix_before"][r][c] != res["matrix_after"][r][c]:
                row.append("-")
            else:
                row.append("#" if res["matrix_before"][r][c] else ".")
        lines.append(f"{r:3d} " + "".join(row))
    legend = ["  #  dark now, dark after -- unchanged",
              "  .  light now, light after -- unchanged",
              "  X  FILL THIS: light now, dark after (one pen stroke)"]
    if cfg.preblank:
        legend.append("  o  printed light on purpose (dark in a plain symbol); "
                      "the before-scan's ECC repairs it")
    else:
        legend.append("  D  left dark: unreachable with a pen, so the "
                      "after-scan's ECC repairs it")
    legend.append("  -  left light: after-scan's ECC repairs this codeword")
    return "\n".join(lines), legend


def svg(codec, res, cfg, path):
    scale, border, gap = 16, 4, 3
    panel = (codec.size + 2 * border) * scale
    w = panel * 2 + gap * scale
    h = panel + 60
    fills = sorted({tuple(c) for c in res["fill_cells"]})
    printed, final = printed_matrix(res, cfg), final_matrix(res, cfg)

    hl = "#" + cfg.highlight.lstrip("#")

    def panel_svg(ox, grid, title, emphasise=()):
        out = [f'<text x="{ox}" y="18" class="t">{title}</text>']
        for r in range(codec.size):
            for c in range(codec.size):
                x, y = ox + (c + border) * scale, 28 + (r + border) * scale
                out.append(f'<rect x="{x}" y="{y}" width="{scale}" '
                           f'height="{scale}" '
                           f'fill="{"#111" if grid[r][c] else "#fff"}"/>')
        # Same outline as the PNG: the stroke is kept INSIDE the module and no
        # line is drawn across an edge shared with another fill cell.  Corners
        # here are square throughout -- this is the schematic panel; the
        # printable PNG rounds every corner except the ones where two cells
        # meet diagonally.
        lw = 2
        ins = lw / 2                 # keeps the stroke inside the module edge
        marked = set(emphasise)

        def pull(r, c, orth):
            """Extend an edge end to close a convex corner, or leave it for the
            neighbouring cell's own edge to cover."""
            return ins if orth not in marked else 0

        for r, c in emphasise:
            x0 = ox + (c + border) * scale
            y0 = 28 + (r + border) * scale
            x1, y1 = x0 + scale, y0 + scale
            if (r - 1, c) not in marked:
                a = x0 + pull(r, c, (r, c - 1))
                b = x1 - pull(r, c, (r, c + 1))
                out.append(f'<line x1="{a}" y1="{y0 + ins}" x2="{b}" '
                           f'y2="{y0 + ins}" stroke="{hl}" '
                           f'stroke-width="{lw}"/>')
            if (r + 1, c) not in marked:
                a = x0 + pull(r, c, (r, c - 1))
                b = x1 - pull(r, c, (r, c + 1))
                out.append(f'<line x1="{a}" y1="{y1 - ins}" x2="{b}" '
                           f'y2="{y1 - ins}" stroke="{hl}" '
                           f'stroke-width="{lw}"/>')
            if (r, c - 1) not in marked:
                a = y0 + pull(r, c, (r - 1, c))
                b = y1 - pull(r, c, (r + 1, c))
                out.append(f'<line x1="{x0 + ins}" y1="{a}" '
                           f'x2="{x0 + ins}" y2="{b}" stroke="{hl}" '
                           f'stroke-width="{lw}"/>')
            if (r, c + 1) not in marked:
                a = y0 + pull(r, c, (r - 1, c))
                b = y1 - pull(r, c, (r + 1, c))
                out.append(f'<line x1="{x1 - ins}" y1="{a}" '
                           f'x2="{x1 - ins}" y2="{b}" stroke="{hl}" '
                           f'stroke-width="{lw}"/>')
        return "".join(out)

    body = [
        panel_svg(0, printed, "PRINT THIS (outlined = cells to ink in)",
                  emphasise=fills),
        panel_svg(panel + gap * scale, final,
                  "AFTER FILLING (scans as the 'after' url)"),
    ]
    doc = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
           f'viewBox="0 0 {w} {h}">'
           f'<style>.t{{font:bold 13px sans-serif;fill:#111}}</style>'
           f'<rect width="{w}" height="{h}" fill="#fff"/>' + "".join(body) +
           f'<text x="0" y="{panel + 44}" class="t">'
           f'{res["before"]} &#8594; fill {res["fills"]} cells (mask '
           f'{res["mask"]}) &#8594; {res["after"]}</text></svg>')
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)


def _outline_mask(side, cells, scale, border, width=2):
    """Outline for the fill set, drawn INSIDE the pixel it marks.

    The ring is the union of the cell squares minus that union eroded by
    `width`, so its outer edge is exactly the module boundary and the stroke
    never leaves the pixel it marks.  It cannot reach a neighbour, and it shows
    the true size of the pixel to ink.  Two consequences worth knowing:

      * cells sharing an edge get no line between them, because that edge is
        interior to the union, so their white interiors merge into one region;
      * cells meeting only at a corner keep a SQUARE corner there instead of a
        rounded one, so the two outlines meet and the pair reads as one
        connected region.  A rounded corner would pull both sides back and
        leave a diagonal pair looking like two separate cells.  Rounding is
        applied to the solid union and those corners are restored out of it, so
        every restored pixel is a pixel of a fill cell -- nothing can leak into
        a module that is not being filled.
    """
    from PIL import Image, ImageChops, ImageDraw, ImageFilter

    union = Image.new("L", (side, side), 0)
    du = ImageDraw.Draw(union)
    for (r, c) in cells:
        x, y = (c + border) * scale, (r + border) * scale
        du.rectangle([x, y, x + scale - 1, y + scale - 1], fill=255)

    radius = max(1.0, scale * 0.12)
    rounded = union.filter(ImageFilter.GaussianBlur(radius))
    rounded = rounded.point(lambda v: 255 if v >= 128 else 0)
    # Blur-then-threshold rounds the corners but can also tip a pixel just
    # outside the union over the line at a shallow step.  Clipping back to the
    # union makes "never leaves the pixel it marks" structural rather than a
    # property to hope for: ring = rounded - core is then always inside a fill
    # cell, whatever the rounding did.
    rounded = ImageChops.multiply(rounded, union)

    # Big enough to undo the rounding at a contact corner; clipping it to the
    # union below is what keeps it off every other module.
    half = max(2, int(round(scale * 0.18)))
    cells = set(cells)
    for (r, c) in cells:
        for dr, dc in ((1, 1), (1, -1)):
            if (r + dr, c + dc) not in cells:
                continue
            if (r, c + dc) in cells or (r + dr, c) in cells:
                continue            # edge-connected already, union joins it
            px = (c + (1 if dc > 0 else 0) + border) * scale
            py = (r + 1 + border) * scale
            box = (px - half, py - half, px + half, py + half)
            # `union` crops the sharp corner patch down to fill-cell pixels
            rounded.paste(ImageChops.lighter(rounded.crop(box),
                                             union.crop(box)), box[:2])

    core = rounded.filter(ImageFilter.MinFilter(2 * width + 1))
    return ImageChops.subtract(rounded, core)


def pngs(codec, res, cfg, directory):
    """Write the printable and its companions.

    fillmap.png is the artefact to print and fill: the symbol with a light
    outline on each cell to fill, and nothing else -- no numbers, no
    separate sheet to consult.  Outlines rather than solid blocks keep the
    module centre white, and the centre is exactly where a decoder samples, so
    the highlighted sheet still scans as the 'before' url while reading like a
    scantron: here the pretty option and the robust option are the same one.
    Modules stay exact integer squares -- no anti-aliasing and no resampling
    anywhere -- so the printer has nothing to interpolate.  The quiet zone is
    the ISO 18004 4-module minimum and the DPI is embedded so print dialogs
    show the true physical size.
    """
    from PIL import Image, ImageDraw

    printed = printed_matrix(res, cfg)
    fills = sorted({tuple(c) for c in res["fill_cells"]})
    n, border = codec.size, 4
    hl = cfg.highlight.lstrip("#")
    ink = tuple(int(hl[i:i + 2], 16) for i in (0, 2, 4))

    def render(path, scale, highlight=False, filled=False):
        side = (n + 2 * border) * scale
        img = Image.new("RGB", (side, side), "white")
        d = ImageDraw.Draw(img)
        pad = max(1, scale // 10)
        for r in range(n):
            for c in range(n):
                if printed[r][c]:
                    x, y = (c + border) * scale, (r + border) * scale
                    d.rectangle([x, y, x + scale - 1, y + scale - 1], fill=0)
        if highlight and fills:
            mask = _outline_mask(side, fills, scale, border,
                                 max(2, int(round(scale * 0.12))))
            img.paste(Image.new("RGB", (side, side), ink), (0, 0), mask)
        if filled:
            # what the sheet looks like once the outlined cells are inked in
            for r, c in fills:
                x, y = (c + border) * scale, (r + border) * scale
                d.rectangle([x - pad, y - pad, x + scale - 1 + pad,
                             y + scale - 1 + pad], fill=0)
        img.save(path, dpi=(cfg.dpi, cfg.dpi))
        return side

    # fillmap.png IS the printable: one artefact, no numbers, no separate sheet
    # to consult.  after-filling.png shows the same sheet once inked in.
    made = [(os.path.join(directory, "fillmap.png"), res["before"]),
            (os.path.join(directory, "after-filling.png"), res["after"])]
    render(made[0][0], cfg.png_scale, highlight=True)
    render(made[1][0], cfg.png_scale, highlight=True, filled=True)
    return made


def verify(codec, res, cfg, files=None):
    """Decode the winning pair. Only ever called on the winner.

    With `files` (path, expected url) pairs it decodes the PNGs actually written
    to disk, which is the artefact that goes to the printer; otherwise it builds
    the two images in memory.
    """
    try:
        import cv2
        import numpy as _np
    except ImportError:
        return "verify skipped: opencv not available"

    det = cv2.QRCodeDetector()

    def decode_rows(rows):
        a = (_np.array(rows, dtype=_np.uint8) ^ 1) * 255
        n = a.shape[0]
        img = _np.full(((n + 8) * 10,) * 2, 255, dtype=_np.uint8)
        img[40:40 + n * 10, 40:40 + n * 10] = _np.kron(
            a, _np.ones((10, 10), dtype=_np.uint8))
        return det.detectAndDecode(img)[0]

    if files:
        checks = [(os.path.basename(p),
                   det.detectAndDecode(cv2.imread(p, cv2.IMREAD_GRAYSCALE))[0],
                   want) for p, want in files]
    else:
        checks = [("printed (in memory)", decode_rows(printed_matrix(res, cfg)),
                   res["before"]),
                  ("filled (in memory)", decode_rows(final_matrix(res, cfg)),
                   res["after"])]
    good = all(got == want for _, got, want in checks)
    detail = "; ".join(f"{name} reads {got!r} (want {want!r})"
                       for name, got, want in checks)
    return f"decoder check [{'OK' if good else 'MISMATCH'}]: {detail}"


def selftest():
    """Check our bit->module map and mask machinery against segno.

    The finished matrix carries the data mask, so XOR it back out before
    comparing -- otherwise this checks the wrong thing.
    """
    print("selftest: does our bit->module map agree with segno's own bit stream?")
    ok = True
    for version in (1, 2, 3, 4):
        for name in "LMQH":
            if any(e.num_blocks != 1 for e in C.ECC[version][ERROR_LEVEL[name]]):
                print(f"  v{version}-{name}: multi-block, out of scope")
                continue
            url = "lithic.uk/" + "x" * version
            codec = Codec(version, name)
            try:
                bits = raw_message_bits(url, version, ERROR_LEVEL[name])
                canon = codec.codeword_bits(url, 0)
            except E.DataOverflowError:
                print(f"  v{version}-{name}: payload does not fit, skipped")
                continue
            n = codec.n_bits
            want = sum(bit << i for i, bit in enumerate(bits[:n]))
            same = canon == want
            masks_ok = all(
                canon ^ codec.mask_bits[m] == codec.printed_bits(
                    codec.matrix(url, m)) for m in range(8))
            ok &= same and masks_ok
            print(f"  v{version}-{name}: {n} codeword bits + {codec.remainder} "
                  f"filler over {codec.n_slots} slots -> "
                  f"{'MATCH' if same else 'MISMATCH'}, 8/8 masks "
                  f"{'consistent' if masks_ok else 'INCONSISTENT'} "
                  f"(repairs {codec.correctable}, d>={codec.min_symbol_distance})")
    print("selftest:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
def forced_pair(codec, urls, canon, cfg):
    """--pair BEFORE,AFTER: detail exactly this pair, without searching.

    Comparing two candidates -- or re-printing one found earlier -- should not
    cost a full scan.  Everything downstream (flip map, printable, --verify) is
    the ordinary code path, so a forced pair is emitted exactly like a winner.
    """
    parts = [s.strip() for s in cfg.pair.split(",")]
    if len(parts) != 2 or parts[0] not in urls or parts[1] not in urls:
        raise SystemExit("--pair needs 'BEFORE,AFTER', both from the candidate "
                         f"family; got {cfg.pair!r}")
    i, j = urls.index(parts[0]), urls.index(parts[1])
    best = None
    for m in cfg.masks:
        detail = analyse(codec, urls, canon, i, j, cfg, m)
        if detail and (best is None or detail["score"] < best["score"]):
            best = detail
    if best is None:
        raise SystemExit(f"{parts[0]} -> {parts[1]} is not pen-reachable under "
                         f"mask(s) {cfg.masks}")
    print(f"\nforced pair (--pair): {parts[0]} -> {parts[1]}  mask "
          f"{best['mask']}, {best['fills']} fills in {best['regions']} "
          f"region(s)")
    stat = dict(pairs=1, reachable_pairs=1, reachable_combos=1,
                best_raw=best["raw"], best_raw_pair=(parts[0], parts[1]),
                meets_target=int(best["fills"] <= cfg.max_fills),
                requested=cfg.pair)
    return best, [best], stat


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base", default="lithic.uk/",
                   help="prefix used verbatim; include any trailing slash")
    p.add_argument("--suffix-chars",
                   default="abcdefghijklmnopqrstuvwxyz0123456789")
    p.add_argument("--suffix-len", type=int, default=2)
    p.add_argument("--urls", help="file of candidate urls, one per line")
    p.add_argument("--ecc", default="L", choices=list(ERROR_LEVEL),
                   help="error correction level; L gives the closest pairs")
    p.add_argument("--version", type=int, help="force a QR version")
    p.add_argument("--mask", default="auto",
                   help="'auto' tries all 8 masks per pair (recommended: the "
                        "mask cannot move differing modules but it decides "
                        "which of them a pen can reach, and how they clump), "
                        "or pin one mask 0-7")
    p.add_argument("--max-fills", type=int, default=4,
                   help="how many fills you are willing to do by hand "
                        "(the collision threshold)")
    p.add_argument("--blob-weight", type=float, default=1.0,
                   help="cost of one extra fill region; raise it to prefer "
                        "fewer, larger regions")
    p.add_argument("--orth-weight", type=float, default=1.5,
                   help="bonus per edge-sharing (orthogonal) fill pair.  A "
                        "corner-only join already halves a region, this makes "
                        "sharing an edge worth more still, so an orthogonal "
                        "clump outranks extra scattered fills")
    p.add_argument("--highlight", default="#ffd9ff",
                   help="outline colour for the printable. The default is a "
                        "tint of ONE ink (light magenta), because a printer "
                        "has to MIX a light yellow from yellow + magenta and "
                        "that mix is what drifts pink; a single-ink tint has "
                        "nothing to mix. Run swatches.py to test candidates "
                        "on your own printer")
    p.add_argument("--margin", type=int, default=0,
                   help="repair budget to leave unused (1 = one symbol of "
                        "slack against print noise)")
    p.add_argument("--no-preblank", action="store_true",
                   help="print a pristine 'before' symbol instead of "
                        "pre-erasing the unreachable modules")
    p.add_argument("--no-slack", action="store_true",
                   help="never spend repair budget to skip fill work")
    p.add_argument("--max-candidates", type=int, default=2500)
    p.add_argument("--max-raw", type=int,
                   help="ignore pairs whose raw module distance exceeds this")
    p.add_argument("--refine", type=int, default=400,
                   help="how many cheapest (pair, mask) results get the "
                        "clumping analysis")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--show", type=int, default=5,
                   help="how many alternatives to list")
    p.add_argument("--emit", metavar="DIR",
                   help="write report.json, flipmap.svg and print-ready PNGs")
    p.add_argument("--png-scale", type=int, default=16,
                   help="pixels per module in the emitted PNGs")
    p.add_argument("--dpi", type=int, default=300,
                   help="dpi recorded in the PNGs so print dialogs show the "
                        "true physical size")
    p.add_argument("--pair", metavar="BEFORE,AFTER",
                   help="skip the search and emit exactly this pair from the "
                        "candidate family (with --mask N, or all 8 masks)")
    p.add_argument("--verify", action="store_true",
                   help="decode the winning pair with OpenCV")
    p.add_argument("--selftest", action="store_true")
    cfg = p.parse_args(argv)
    cfg.preblank = not cfg.no_preblank
    cfg.slack = not cfg.no_slack

    if cfg.selftest:
        return selftest()

    if cfg.mask == "auto":
        cfg.masks = list(range(8))
    else:
        m = int(cfg.mask)
        if not 0 <= m <= 7:
            raise SystemExit("--mask must be 'auto' or 0-7")
        cfg.masks = [m]

    urls = build_candidates(cfg)
    if len(urls) < 2:
        raise SystemExit("need at least two candidates")
    version = cfg.version or pick_version(urls, cfg.ecc)
    codec = Codec(version, cfg.ecc)
    print(f"{len(urls)} candidates, {len(urls)*(len(urls)-1)//2:,} pairs | "
          f"v{version}-{cfg.ecc} | {codec.n_codewords} codewords "
          f"({codec.parity} parity) | repairs {codec.correctable} symbol "
          f"errors | min symbol distance {codec.min_symbol_distance}")
    print(f"masks: {'all 8' if len(cfg.masks) == 8 else cfg.masks} | "
          f"mode: {'pre-blank' if cfg.preblank else 'pristine print'} | "
          f"slack {'on' if cfg.slack else 'off'} | margin {cfg.margin} fills")

    canon, kept, fillers = [], [], set()
    for u in urls:
        try:
            rows = codec.matrix(u, 0)
        except E.DataOverflowError:
            continue
        canon.append(codec.printed_bits(rows) ^ codec.mask_bits[0])
        fillers.add(codec.filler_bits(rows))
        kept.append(u)
    urls = kept

    if cfg.pair:
        best, scored, stat = forced_pair(codec, urls, canon, cfg)
    else:
        best, scored, stat = scan(codec, urls, canon, cfg)
        print(f"\nscanned {stat['pairs']:,} pairs x {len(cfg.masks)} mask(s); "
              f"{stat['reachable_pairs']:,} pairs are pen-reachable under some "
              f"mask ({stat['reachable_combos']:,} pair+mask combinations)")
        print(f"closest pair by raw module distance: {stat['best_raw']} modules "
              f"({stat['best_raw_pair'][0]} vs {stat['best_raw_pair'][1]}) "
              f"-- invariant across masks")
        print(f"results reaching --max-fills {cfg.max_fills}: "
              f"{stat['meets_target']}")

    if best is None:
        print("\nNo pen-reachable pair in this family. Widen the search "
              "(--suffix-len 3, a bigger --max-candidates, or --urls FILE).")
        return 0

    line = "=" * 74
    print(f"\n{line}\nBEST PAIR:  {best['before']}  ->  {best['after']}"
          f"   (mask {best['mask']})\n{line}")
    print(f"  differing modules      {best['raw']}  (up {best['up']} reachable "
          f"with a pen, down {best['down']} not)")
    print(f"  codewords touched      stranded {best['k_down']}, fillable "
          f"{best['k_up']} of a repair budget of {codec.correctable}")
    print(f"  FILL TARGETS           {best['fills']} modules in "
          f"{best['regions']} region(s) a single stroke can cover "
          f"({best['orth_pairs']} orthogonal link(s), "
          f"{best['diag_links']} diagonal-only), largest region "
          f"{best['biggest_region']}")
    print(f"  skipped codewords      {best['skipped_codewords']}")
    print(f"  score                  {best['score']:.2f} = {best['fills']} + "
          f"{cfg.blob_weight}*({best['regions_eff']:.1f}-1) - "
          f"{cfg.orth_weight}*{best['orth_pairs']} orthogonal pair(s)")
    print(f"  pen regions            {best['regions_eff']:.1f} effective "
          f"({best['regions']} separate strokes, {best['diag_links']} of the "
          f"joins only a corner touch)")
    print(f"  repair headroom        before-scan {best['before_headroom']}, "
          f"after-scan {best['after_headroom']} "
          f"(of {codec.correctable} symbols each)")
    if min(best["before_headroom"], best["after_headroom"]) == 0:
        if cfg.margin:
            print(f"  !!! budget fully spent: this pair uses all "
                  f"{codec.correctable - cfg.margin} repair(s) it is allowed, "
                  f"so the {cfg.margin} reserved symbol(s) are the only "
                  f"protection against print noise.")
        else:
            print("  !!! no headroom: this pair touches the correctability "
                  "limit, so one stray speck\n      breaks it. Use --margin 1 "
                  "for a version that survives print noise.")
    print(f"  verdict                "
          + ("MEETS --max-fills" if best["fills"] <= cfg.max_fills else
             f"exceeds --max-fills {cfg.max_fills}; this is the best this "
             f"family offers"))

    art, legend = ascii_map(codec, best, cfg)
    print("\nflip map (row 0 at top):")
    print(art)
    print("\nlegend:\n" + "\n".join(legend))
    print(f"\nthe {best['fills']} cells to fill, in order (1-based):")
    for n, bit in enumerate(best["fill_bits"], 1):
        r, c = codec.cell(bit)
        print(f"  {n:3d}. r{r+1:>3} c{c+1:>3}   codeword {bit//8}")

    if len(scored) > 1:
        print(f"\ntop {min(cfg.show, len(scored))} by score:")
        for d in scored[:cfg.show]:
            print(f"  {d['before']} -> {d['after']}  mask {d['mask']}   "
                  f"fills {d['fills']:2d}  regions {d['regions']:2d} "
                  f"(pen-eff {d['regions_eff']:.1f}, diag-only {d['diag_links']}, "
                  f"orth {d['orth_pairs']})  raw {d['raw']:3d}  "
                  f"score {d['score']:.2f}")

    png_list = []
    if cfg.emit:
        os.makedirs(cfg.emit, exist_ok=True)
        png_list = pngs(codec, best, cfg, cfg.emit)
        svg(codec, best, cfg, os.path.join(cfg.emit, "flipmap.svg"))
        payload = {k: v for k, v in best.items() if not k.startswith("matrix")}
        payload["codec"] = dict(version=version, ecc=cfg.ecc, mask=best["mask"],
                                size=codec.size,
                                correctable=codec.correctable,
                                parity=codec.parity,
                                mode="preblank" if cfg.preblank else "pristine")
        payload["stats"] = stat
        with open(os.path.join(cfg.emit, "report.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        side = (codec.size + 8) * cfg.png_scale
        print(f"\nwrote {cfg.emit}/report.json, flipmap.svg and "
              f"{len(png_list) + 1} pngs")
        print(f"PRINT THIS: {cfg.emit}/fillmap.png -- {side}x{side}px = "
              f"{side / cfg.dpi * 25.4:.1f} mm square at {cfg.dpi} dpi "
              f"({cfg.png_scale}px per module, 4-module quiet zone)")
        print(f"{cfg.highlight} outlines mark the cells to ink in; no numbers "
              "on the sheet")

    if cfg.verify:
        print("\n" + verify(codec, best, cfg, png_list or None))
    return 0


if __name__ == "__main__":
    sys.exit(main())
