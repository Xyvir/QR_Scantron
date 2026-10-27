"""Is the mask choice really irrelevant?

Raw distance: yes, provably.  Both symbols are XORed with the same pattern, so
the set of differing modules cannot move.

What a pen can fix: no.  Position p is fillable iff the printed module is light,
i.e. P_m(p) == C_A(p); filling then drives the *unmasked* bit to 1 XOR P_m(p).
So where P_m(p) == 1 a light->dark stroke pushes the decoded bit DOWN, repairing
a difference that mask 0 can only repair by erasing.  Masks therefore re-partition
the same differing positions into fixable / unfixable, which changes the fill set
and hence its clumping.

This checks the mask-0 case reproduces the tool exactly, then measures the gain.
"""
import itertools
import types

from qr_collider import Codec, analyse, blobs

c = Codec(1, "L")
NB = c.n_bits
MSK = c.bit_mask
P = c.mask_bits
assert len(P) == 8 and P[0] != 0, "expected 8 distinct mask patterns"

URLS = ["lithic.uk/" + x for x in "abcdefghijklmnopqrstuvwxyz0123456789"]


def codeword(url):
    """Canonical unmasked codeword array (mask 0 XORed back out)."""
    return c.codeword_bits(url, 0)


C = {u: codeword(u) for u in URLS}
CFG = types.SimpleNamespace(preblank=True, slack=True, margin=0, blob_weight=1.0,
                            orth_weight=1.5)
KEYS = [c.codeword_bits(u, 0) for u in URLS]


def evaluate(a, b, m):
    """Same model as the tool, but for an arbitrary shared mask."""
    A, B = C[a], C[b]
    d = A ^ B
    u0 = d & ~A & MSK               # C_A=0, C_B=1
    u1 = d & A                      # C_A=1, C_B=0
    fill = ((u0 & ~P[m]) | (u1 & P[m])) & MSK
    unfix = d & ~fill & MSK
    kd = sum(1 for i in range(c.n_codewords) if (unfix >> (8 * i)) & 0xFF)
    up = [bin((fill >> (8 * i)) & 0xFF).count("1") for i in range(c.n_codewords)]
    skip = set(sorted(range(c.n_codewords), key=lambda i: -up[i])[:c.correctable])
    cells = [c.order[8 * i + bit] for i in range(c.n_codewords)
             if i not in skip for bit in range(8) if (fill >> (8 * i + bit)) & 1]
    comps = blobs(cells)                     # corner contacts count as joined
    comps4 = blobs(cells, diag=False)        # one stroke can cover these
    cs = set(cells)
    orth = sum(1 for (r, cc) in cs for dr, dc in ((1, 0), (0, 1))
               if (r + dr, cc + dc) in cs)
    return dict(raw=bin(d).count("1"), kd=kd, fills=len(cells), blobs=len(comps),
                regions=len(comps4), orth=orth,
                regions_eff=len(comps) + 0.5 * (len(comps4) - len(comps)),
                biggest=max((len(x) for x in comps), default=0))


pairs = list(itertools.combinations(URLS, 2))
print(f"=== 1. regression: does mask 0 here match the shipped tool exactly? ===")
bad = 0
for i, (a, b) in enumerate(pairs):
    ref = analyse(c, URLS, KEYS, URLS.index(a), URLS.index(b), CFG, 0)
    mine = evaluate(a, b, 0)
    if (ref is None) != (mine["kd"] > c.correctable):
        bad += 1
        print(f"  MISMATCH {a}->{b}: reachability disagrees")
    elif ref is not None and (ref["fills"], ref["blobs8"], ref["k_down"]) != \
            (mine["fills"], mine["blobs"], mine["kd"]):
        bad += 1
        if bad < 4:
            print(f"  MISMATCH {a}->{b}: tool {ref['fills']}/{ref['blobs']}/"
                  f"{ref['k_down']} vs {mine['fills']}/{mine['blobs']}/{mine['kd']}")
print(f"  {len(pairs) - bad}/{len(pairs)} pairs identical")
assert bad == 0, "formula disagrees with the tool; do not trust the rest"

print("\n=== 2. does the fill set (and its clumping) depend on the mask? ===")
moved = 0
for a, b in pairs:
    fillsets = {m: evaluate(a, b, m)["fills"] for m in range(8)}
    if len(set(fillsets.values())) > 1:
        moved += 1
print(f"  {moved}/{len(pairs)} pairs change their fill COUNT under some mask")
blobmoved = sum(1 for a, b in pairs
                if len({evaluate(a, b, m)["blobs"] for m in range(8)}) > 1)
print(f"  {blobmoved}/{len(pairs)} pairs change their BLOB count under some mask")

print("\n=== 3. is any pair fixable with a pen only under a nonzero mask? ===")
only_nonzero = []
for a, b in pairs:
    reach = [m for m in range(8) if evaluate(a, b, m)["kd"] <= c.correctable]
    if reach and 0 not in reach:
        only_nonzero.append((a, b, reach, evaluate(a, b, reach[0])["fills"]))
print(f"  {len(only_nonzero)} pairs need a nonzero mask to be pen-reachable")
for a, b, reach, f in only_nonzero[:3]:
    print(f"    {a} -> {b}: reachable with mask(s) {reach}, {f} fills")

print("\n=== 4. best available result, mask 0 only vs best of all 8 masks ===")
print("  (the tool's default score: fills + 1.0*(pen_regions-1) - 1.5*orth_pairs)")
for label, masks in (("mask 0 only", [0]), ("best of 8 masks", range(8))):
    best = None
    for a, b in pairs:
        for m in masks:
            r = evaluate(a, b, m)
            if r["kd"] > c.correctable:
                continue
            score = (r["fills"] + 1.0 * (r["regions_eff"] - 1)
                     - 1.5 * r["orth"])
            if best is None or score < best[0]:
                best = (score, a, b, m, r)
    score, a, b, m, r = best
    print(f"  {label:16s}: {a} -> {b} mask {m}  fills {r['fills']:2d}  "
          f"pen_regions {r['regions_eff']:.1f}  biggest {r['biggest']}  "
          f"score {score:.1f}  raw {r['raw']}")
print("  (in this small 1-char space mask 0 already holds the best result; what")
print("   the mask search buys is reachability -- section 3 -- and it buys much")
print("   more of it in the larger 2-char space the tool actually searches)")

print("\n=== 5. the 3-fill pair from the big search, mask by mask ===")
A, B = "lithic.uk/e4", "lithic.uk/zt"       # cheapest pair in the 2-char space
for _u in (A, B):
    C[_u] = codeword(_u)
print(f"  differing modules (invariant) = {evaluate(A, B, 0)['raw']}")
print(f"  {'mask':>4} {'k_down':>6} {'fills':>5} {'blobs':>5} {'biggest':>7}  reachable")
for m in range(8):
    r = evaluate(A, B, m)
    print(f"  {m:>4} {r['kd']:>6} {r['fills']:>5} {r['blobs']:>5} "
          f"{r['biggest']:>7}  {r['kd'] <= c.correctable}")
