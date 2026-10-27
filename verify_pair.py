"""Independent check of the headline result, using report.json and segno/cv2 only.

Confirms:
  * the printed 'before' symbol really is A with only its dark->light modules
    erased, and really decodes as the 'before' url
  * filling exactly the listed cells yields a symbol that decodes as 'after'
  * how much print noise that zero-margin solution survives
"""
import json
import random
import sys

import cv2
import numpy as np
import segno

CFG = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "out/report.json"))
VERSION = CFG["codec"]["version"]
ERROR = CFG["codec"]["ecc"]
MASK = CFG["codec"]["mask"]
FILLS = [tuple(c) for c in CFG["fill_cells"]]
BEFORE, AFTER = CFG["before"], CFG["after"]
SIZE = CFG["codec"]["size"]


def matrix(url):
    q = segno.make_qr(url, version=VERSION, error=ERROR, mask=MASK,
                      boost_error=False)
    return [[int(v) for v in row] for row in q.matrix]


def decode(rows, scale=10):
    a = (np.array(rows, dtype=np.uint8) ^ 1) * 255
    n = a.shape[0]
    img = np.full(((n + 8) * scale,) * 2, 255, dtype=np.uint8)
    img[4 * scale:4 * scale + n * scale, 4 * scale:4 * scale + n * scale] = \
        np.kron(a, np.ones((scale, scale), dtype=np.uint8))
    return cv2.QRCodeDetector().detectAndDecode(img)[0]


A, B = matrix(BEFORE), matrix(AFTER)
P = [[A[r][c] & B[r][c] for c in range(SIZE)] for r in range(SIZE)]
F = [row[:] for row in P]
for r, c in FILLS:
    F[r][c] = 1


def diff(x, y):
    return sum(1 for r in range(SIZE) for c in range(SIZE) if x[r][c] != y[r][c])


print(f"claim: {BEFORE} -> {AFTER}, {len(FILLS)} fills, v{VERSION}-{ERROR} "
      f"mask {MASK}")
print(f"  printed vs A (dark->light erased): {diff(P, A)} modules   "
      f"(vs B: {diff(P, B)})")
print(f"  filled  vs B (left for ECC):       {diff(F, B)} modules   "
      f"(vs A: {diff(F, A)})")
print(f"  printed decodes as: {decode(P)!r}  want {BEFORE!r}")
print(f"  filled  decodes as: {decode(F)!r}  want {AFTER!r}")
print(f"  also at scale 20:   {decode(F, 20)!r}")

# function modules for version 1: finders+separators, format areas, timing
def is_function(r, c):
    if r < 9 and c < 9:
        return True
    if r < 9 and c > SIZE - 10:
        return True
    if r > SIZE - 10 and c < 9:
        return True
    return r == 6 or c == 6


free = [(r, c) for r in range(SIZE) for c in range(SIZE)
        if not is_function(r, c) and (r, c) not in FILLS]
print(f"\nprint-noise tolerance (random extra specks in {len(free)} data modules), "
      f"20 trials each:")
rng = random.Random(7)
for k in range(0, 5):
    ok_p = ok_f = 0
    for _ in range(20):
        spots = rng.sample(free, k)
        p = [row[:] for row in P]
        f = [row[:] for row in F]
        for r, c in spots:
            p[r][c] ^= 1
            f[r][c] ^= 1
        ok_p += decode(p) == BEFORE
        ok_f += decode(f) == AFTER
    print(f"  {k} extra speck(s): before-scan correct {ok_p:2d}/20, "
          f"after-scan correct {ok_f:2d}/20")
