# QR pen-flip collision searcher

Finds two URLs in a family — `lithic.uk/e4` and `lithic.uk/zt`, say — whose QR
symbols are close enough that a **human with a pen can change where the code
points**. Print one symbol, scan it, ink in a few outlined cells, scan again, and
you land somewhere else.

The output is a printable sheet: the "before" symbol with a light outline on
every cell you need to fill in, in the cheapest position the search could find.

```
python3 qr_collider.py --emit out --verify
```

---

## Why this is possible at all

You never have to turn symbol A into symbol B exactly. A QR symbol carries
Reed–Solomon parity, so a scanner **repairs up to `t` damaged codewords**
(`t = parity_codewords // 2`) before it even reads the payload. At v1-L that is
RS(26,19): 26 codewords, 7 parity, **3 repairable symbol errors**, minimum symbol
distance 8. Three consequences, and the tool exploits all three:

- **Dark→light modules never need erasing.** Leave them dark. The *after* scan
  repairs them, as long as they fall inside its repair budget. A pen that only
  adds black is enough.
- **Light→dark work can be skipped in up to `t` codewords.** A codeword's 8
  modules sit in one small 2×4 block, so skipping the *heaviest* codewords
  deletes whole clumps of pen work for one unit of budget each.
- **Pre-blanking splits the two budgets.** Instead of printing A, print A with
  the unreachable modules already erased. It still scans as A, and the after
  scan can then spend its whole budget skipping fill work. `--margin` buys
  extra headroom against print noise.

**What is impossible** is the naive framing: "flip 3–4 pixels and A becomes B".
Measured over all 839,160 pairs in the worked example below, the closest two
symbols differ in **20 modules**. Pen strokes are not codeword edits, which is
exactly why the trick works.

### The mask is a lever, but not for the usual reason

Two symbols sharing a mask are XORed with the same pattern, so masking cancels
out of the pairwise difference: **the set of differing modules cannot be moved by
any mask.** What the mask *does* decide is which of those modules a pen can
reach — a printed module is fillable iff it is light on paper, and filling it
drives the decoded bit to `1 XOR P_m(p)`. So the mask re-partitions the same
differing modules into fixable and stranded, changing the fill set and how it
clumps. Searching all 8 masks widens reachability by ~3.5× (258 → 894 pairs in
the worked example). See `masks.py` for the measurements.

---

## Requirements

| need | tested with | note |
|---|---|---|
| Python | 3.14.5 | no syntax newer than f-strings |
| `numpy` | 2.5.4 | uses `np.bitwise_count` when present, with a fallback |
| `segno` | 1.6.6 | pure-Python QR encoder; the search reads its symbol matrices |
| `Pillow` | 12.3.0 | 1-bit PNG output |
| `opencv-python` | 5.0.0 | **optional**, only for `--verify` |

```bash
pip install numpy segno pillow opencv-python
```

## Quick start

```bash
python3 qr_collider.py --selftest            # checks our bit -> module map
python3 qr_collider.py                       # search the default family
python3 qr_collider.py --emit out --verify   # write the printable + decode it
```

Print `out/fillmap.png` at **100% / actual size** (not "fit to page"). Scan it
before you touch it, then fill the outlined cells and scan again.

## What `--emit DIR` writes

| file | what it is |
|---|---|
| `fillmap.png` | **the printable.** The before-symbol with a light outline on each cell to fill. No numbers, no second sheet to consult. |
| `after-filling.png` | the same sheet once those cells are inked — the state you are aiming for |
| `report.json` | the pair, mask, fill cells, repair headroom, clump metrics (matrices stripped) |
| `flipmap.svg` | two-panel vector schematic. Square corners; the PNG is the artefact that matters |

Only the winning pair is ever rasterised. The search itself runs on packed
integers — one symbol is one Python `int`, so "how different are these two
symbols" is a popcount. No images, no PNGs, and no decoder in the search loop.

## Why the printable looks the way it does

- The stroke is drawn **inside** the module it marks, so it can never touch a
  neighbouring module and it shows the true size of the pixel to ink.
- **No line is drawn across an edge shared by two fill cells** — their interiors
  merge into one region, which is what you want to shade in.
- Corners are **rounded, except where two cells meet diagonally**. Those corners
  stay square so the two outlines actually meet and a diagonal pair reads as one
  region instead of two boxes that happen to be near each other.
- The colour is a **tint of a single ink** (default `#ffd9ff`, light magenta).
  This is not a preference: a light yellow such as `#ffd633` is yellow *plus
  magenta*, so the printer must mix two inks, and when that mix is off the sheet
  comes out pink. A single-ink tint has nothing to mix, so its hue cannot drift.
  Run `swatches.py` to test candidates on your own printer.
- Output is **1-bit PNG** at an exact integer module size with embedded DPI, so
  the printer has nothing to interpolate. Default: 16 px/module at 300 dpi =
  464×464 px = **39.3 mm** square, with the ISO 18004 4-module quiet zone.

## How candidates are ranked

```
score = fills + blob_weight * (pen_regions - 1) - orth_weight * orth_pairs
        (defaults: blob_weight 1.0, orth_weight 1.5)
```

`pen_regions` counts groups a **single uninterrupted stroke** can cover
(edge-sharing), but a **corner-only join counts as half a merge** — so a
diagonal pair is scored as *1.5 regions*, not a full 1. Combining that with a
per-edge bonus gives the ordering a hand actually wants:

**orthogonal clump < diagonal chain < scattered fills.**

The trade-off is deliberate: the tool will happily take one or two more fills if
they arrive as a clump you can shade in one movement. `--blob-weight` and
`--orth-weight` move that line; `--max-fills` is your hard willingness threshold.

## Options worth knowing

| flag | why you'd touch it |
|---|---|
| `--base`, `--suffix-len`, `--suffix-chars` | the URL family. Default: `lithic.uk/` + 36 suffixes → 1,296 URLs |
| `--urls FILE` | search an explicit list instead |
| `--ecc {L,M,Q,H}` | `L` gives the closest pairs (most repair budget per payload) |
| `--mask auto\|0-7` | `auto` tries all 8; pin one only to reproduce an older result |
| `--max-fills N` | how many cells you are willing to fill by hand |
| `--margin N` | reserve N repairs against print noise — **the difference between a sheet that works and one that dies at the first speck** |
| `--blob-weight`, `--orth-weight` | clumping vs count, as above |
| `--pair BEFORE,AFTER` | skip the search and emit one specific pair (with `--mask N`) |
| `--highlight #rrggbb` | outline colour; keep it a single-ink tint |
| `--png-scale`, `--dpi` | print size. 16 px/module at 300 dpi is 39.3 mm |
| `--no-preblank`, `--no-slack` | turn off the two budget tricks, for comparison |
| `--emit DIR`, `--verify`, `--show N` | write artefacts, decode the winner, list alternatives |
| `--selftest` | verifies the bit→module map against segno's own stream |

## Helper scripts

| script | what it does |
|---|---|
| `swatches.py` | writes `printer-swatches.png`, a calibration strip: single-ink colour candidates, each shown as a solid patch **and** as the real fill outline, so you can pick one that survives your printer |
| `verify_pair.py out/report.json` | independent re-check of a pair (its own code path) plus **print-noise tolerance**: how many stray specks the sheet survives, 20 trials each |
| `masks.py` | the evidence for the mask claims: a mask-0 regression against the shipped tool, and a mask-by-mask table for one pair |

---

## Worked example: the `lithic.uk/` search

Setup: `--base lithic.uk/`, suffixes `a–z0–9`, 2 characters → **1,296 candidates,
839,160 pairs**, v1-L, all 8 masks, pre-blank on.

**Scan**

- **894** of 839,160 pairs are pen-reachable under some mask (920 pair+mask
  combinations). With `--mask 0` only: **258**.
- Closest pair by raw module distance: **20 modules** (`lithic.uk/aa` vs
  `lithic.uk/dh`) — invariant across masks.
- Pairs reaching `--max-fills 4`: **3**.

**Top 5 by the default score**

| before | after | mask | fills | pen-regions | shared edges | corner-only joins | raw | score |
|---|---|---|---|---|---|---|---|---|
| `lithic.uk/3u` | `lithic.uk/5g` | 6 | **5** | 3.0 | 2 | 0 | 22 | **4.00** |
| `lithic.uk/j9` | `lithic.uk/za` | 0 | 8 | 3.0 | 4 | 2 | 28 | 4.00 |
| `lithic.uk/e4` | `lithic.uk/zt` | 5 | **3** | 2.5 | 0 | 1 | 24 | 4.50 |
| `lithic.uk/j8` | `lithic.uk/zl` | 0 | 8 | 3.5 | 4 | 1 | 26 | 4.50 |
| `lithic.uk/dx` | `lithic.uk/q8` | 0 | **3** | 3.0 | 0 | 0 | 22 | 5.00 |

Note what the ranking does here. The cheapest sheet by raw cell count, `e4 → zt`
at 3 fills, **does not win**: one of its joins is diagonal-only, so it scores as
2.5 regions. `3u → 5g` costs two more fills but they arrive as one 3-cell run
plus a pair of clumps, worth 3.0 regions with 2 shared edges — the easier thing
for a hand to shade, and it ranks first.

**What each result actually does**

`lithic.uk/3u → lithic.uk/5g` (mask 6) — the top pick
5 fills in 3 regions: one vertical 3-cell run plus two isolated cells, with 2
shared edges. Decoder-verified from the written PNGs: the outlined sheet reads
`lithic.uk/3u`, the inked sheet reads `lithic.uk/5g` (also at 20 px/module). It
spends the *entire* 3-symbol repair budget on both scans, so print-noise
tolerance is thin: 0 stray specks 20/20, but 1 speck only 4/20 before / 1/20
after.

`lithic.uk/mf → lithic.uk/mj` (mask 3) — the one to actually print

```bash
python3 qr_collider.py --margin 1 --max-fills 12 --blob-weight 2 \
                       --emit out-margin1 --verify
```

13 fills, but in 5 regions with 8 shared edges and a 5-cell run: one long scribble
plus a few dabs. Because `--margin 1` reserves a repair, it **survives a stray
speck 20/20** (and still 2 specks 2/20). This is the trade the sheet exists to
make: more ink, far more forgiving.

`lithic.uk/e4 → lithic.uk/zt` (mask 5) — provably the closest possible

3 fills in 3 regions, one of them diagonal-only. Its two symbols differ in
**exactly 8 codewords**, which is the minimum distance of RS(26,19) — you cannot
find a closer pair at this version and ECC level, so this is a floor, not a
lucky hit. It is also the most fragile sheet in the set: 1 stray speck takes it
to 1/20 before / 2/20 after.

Reproduce any of them with `--pair`:

```bash
python3 qr_collider.py --pair lithic.uk/e4,lithic.uk/zt --mask 5 --emit out-e4 --verify
```

**Why mask 5, specifically?** The differing modules never move, but reachability
does. For `e4 → zt`, `masks.py` reports:

| mask | stranded codewords | fills | 8-connected clumps | reachable? |
|---|---|---|---|---|
| 0 | 7 | 2 | 2 | no (over budget) |
| 2 | 6 | 4 | 4 | no |
| **5** | **3** | **3** | **2** | **yes — the only one** |
| 7 | 7 | 2 | 2 | no |

Raw distance is 24 at every mask; only mask 5 lands inside the repair budget.
At mask 5 those 2 clumps are 3 separate strokes joined by one corner contact,
which is the 2.5 pen-regions the ranking table shows.

---

## Limits, stated plainly

- **Single-block ECC versions only.** Multi-block symbols (e.g. v3-Q) raise a
  clear error rather than guessing at a codeword layout.
- **OpenCV is the only decoder oracle.** It agrees across scale 16 and 20, but a
  real phone camera may differ. Nothing here has been verified on paper.
- The search is exhaustive up to `--max-candidates` (2,500); beyond that it
  samples, and `--seed` makes the sample reproducible.
- The cheapest pairs spend the whole repair budget, so they are brittle by
  construction. `verify_pair.py` measures exactly how brittle; `--margin 1` is
  usually the right answer for a sheet you intend to hand someone.
- Printer colour is environment-specific. The single-ink rule is general; the
  exact tint is not — run `swatches.py` once.

## Files

Commit: `qr_collider.py`, `masks.py`, `verify_pair.py`, `swatches.py`,
`README.md`, `LICENSE`, `.gitignore`.

Generated, safe to ignore: `out*/` (sheets and reports from a run),
`printer-swatches.png`, `__pycache__/` — all covered by `.gitignore`, so a plain
`git add .` picks up exactly the files above.

Console-only runs write nothing at all; PNGs, SVG and JSON appear only with
`--emit`.

## License

MIT — see [LICENSE](LICENSE).
