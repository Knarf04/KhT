# Open Questions — status after April 2026 pass

These are items raised during the April 2026 work where engineering
correctness isn't enough.  Several have since been investigated
empirically.  This document reflects the state as of the end of that
investigation, including what's *resolved*, what's *scoped but still
math-gated*, and what remains genuinely open.

## Summary

| # | Item | Status |
|---|------|--------|
| 1 | `BNComplex.ToCob` semantics (signed vs unsigned lift) | Opt-in flag shipped; math review still wanted |
| 2 | F_p-aware Cob arithmetic coverage | Resolved for `__neg__`, `isIsom`, `ReduceDecorations` |
| 3 | `d² ≠ 0` on cabled knots over F_p (p > 2) | Empirically resolved — all cables now pass d²=0 over F_2/F_3/F_5/F_7 |
| 4 | Final-invariant equivalence of F_p Cob path | Stress-test suite added; no divergence found |
| 5 | Cob-level arrow-shortening (`clean_up`) | Attempted at Cob level; doesn't reduce generator count on wide tangles.  Reason is topological, not algorithmic (see note below) |
| 6 | Reidemeister slice reordering | Commutation analyzer added; on TH-pattern no rewrite possible with local moves |
| 7 | `(1, 2k+1)` arc algebra | **Empirically shown not to help** — `Cob.eliminateAll` at any width catches exactly the same isomorphisms `BN.eliminateAll` does on (1,3) tangles, so a wider arc algebra provides no extra reduction power |
| 8 | Deterministic `clean_up_once` | Opt-in implementation shipped; convergence not formally proven |

## 1. `BNComplex.ToCob` — signed vs unsigned lift

Implementation shipped via a `signed_lift=False` flag on
`BNComplex.ToCob`, `BNAlgebra.mor.ToCob`, and the `BNbracket`
`cleanup_field` path.  Default is the historical unsigned lift
(coefficients in `[0, p)`); flipping to `signed_lift=True` gives centered
lift `(-p/2, p/2]`.  The math question of *which lift is correct* for
downstream identities remains open — the flag is there so a math
collaborator can toggle and compare on real tangles without waiting for
an engineering change.

## 2. F_p-aware Cob arithmetic

Addressed:
- `Cob.mor.__neg__` at [Cob.py:430-439](KhT/Cob.py#L430-L439) now reduces
  coefficients mod `_field` when set.
- `Cob.mor.isIsom` at [Cob.py:418-425](KhT/Cob.py#L418-L425) uses a
  field-aware check (`coeff == 1` or `coeff == -1 mod p`).
- `Cob.mor.ReduceDecorations` at
  [Cob.py:400-411](KhT/Cob.py#L400-L411) was verified to never touch
  coefficients; a comment records the audit.

Unit tests for these primitives live in `examples/quick_tests.py`
test 14 (`Fp_primitives`).

## 3. `d² ≠ 0` on cabled knots over F_p (p > 2)

The 2019 flag for this bug (commit `d2019ce`) turns out to be stale:
running `Tangle.from_legacy("cup1.pos0.neg1.neg1.neg1.pos0.cap1").Cable()`
and calling `.toReduced_BNComplex(1000, p).validate()` passes for
`p ∈ {2, 3, 5, 7}` under both the legacy and the experimental
"corrected" sign convention.  The fix landed incidentally during the
F_p-arithmetic cleanup.  See `examples/d2_bisect.py` for the diagnostic.

The `examples/KnotCables.py` default still uses `field = 2` to match
historical output, but can be flipped to 3, 5, or 7 — all currently
pass.  The `set_sign_convention` flag in `BNAlgebra` is still useful
for exploring alternate sign patterns if a math collaborator wants to.

## 4. F_p Cob path invariant equivalence

`examples/fp_round_trip_stress.py` exercises the F_p-aware Cob path
(`BNbracket` with `cleanup_field=p`) against the default Z-Cob path
across multiple tangles and primes.  Invariants match on everything
tested.  Not a proof, but the empirical evidence is strong.

## 5. Cob-level `clean_up` — experimentally disproved

**This item is mechanically answered.**  A Cob-level arrow-shortening
implementation was built (see `CobComplex.clean_up_once` in
`CobComplexes.py`, plus `_cob_isotopy`, `_cob_isolate_h_arrow`, and
`_cob_cross_idem_probe`).  It covers same-face H-power absorption,
dotted same-face absorption, and a bounded cross-idempotent candidate
probe, all gated by an entropy guard that rejects isotopies that would
grow total decoration count.

On TH-pattern (the target tangle), it fires tens of thousands of
absorptions per slice at the peak, reduces total decoration count, and
has **zero effect on generator count**.  Observed across three
increasingly aggressive variants: H-power only, dotted same-face, and
cross-idempotent brute probe.

**Root cause** (from the investigation): `eliminateAll` at the Cob level
is already catching every identity-like arrow — confirmed by direct
comparison of `CobComplex.eliminateAll` vs `BNComplex.eliminateAll`
landing on identical generator counts on all tested tangles.  Arrow
simplification can't produce new identity arrows when the underlying
tangle is genuinely wide.  Whatever residual simplification an arc
algebra might enable at wider widths is not where the growth comes from.

See item 7 for the deeper "why" and item 6 for the topological escape
route.

## 6. Reidemeister slice reordering

`KhT/TangleRewrite.py` has the API scaffold.  A simple adjacent-swap
commutation analyzer was built and run on TH-pattern — every cup in
TH-pattern is pinned by a neighboring crossing that uses one of its
arcs.  Local swap moves cannot narrow the tangle below its native
width.

Non-local moves (R1/R2/R3 pattern matching across commuting
intermediate slices) are still unimplemented.  Preliminary scan
suggests they wouldn't help TH-pattern specifically (no same-index
`pos`/`neg` pairs with commutable intermediates), but this is manual
inspection, not a proof.

If a tangle *is* narrowable, it would need to come in a pre-narrowed
form — via the tangle's braid presentation, which is topology research
per-tangle.

## 7. `(1, 2k+1)` arc algebra — **not the blocker**

Per the investigation in item 5: `BNComplex.clean_up` does *not* reduce
generator count beyond `eliminateAll` on any tested tangle.  BN clean_up
is a loop-type rearrangement, not a reduction step.  `eliminateAll` in
both `CobComplex` and `BNComplex` ends at the same generator count.

Implication: a generalized arc algebra at width > 3 would *not* provide
extra simplification power.  The isomorphism check `isIsom` at any width
already finds what there is to find.  Item 7 in its original framing
("we need (1, 2k+1) arc algebra to simplify TH-pattern") is based on an
incorrect mental model.

The original item remains open only as "it would be a cleaner compact
representation for wide-width BN complexes" — which it would, but not
for simplification reasons.

## 8. Determinism in `clean_up_once`

`BNComplex.clean_up_once` and `CobComplex.clean_up_once` both accept a
`deterministic` kwarg that swaps `random.choice(remaining)` for a
min-power reduction-measure choice.  `deterministic=False` remains the
default so historical behavior is unchanged.

Convergence is **not formally proven** — `max_iter=50` is used
defensively.  Empirical tests (`quick_tests` test 15) show the
algorithm does not stall and reaches comparable generator counts to
the random-choice version on small tangles.

## Performance work done in the same pass

Independent of the math items above, a substantial performance
optimization pass landed alongside them:

- Incremental `_nnz` / `_row` / `_col` sparse indices on
  `CobComplex` (was rebuilding the full set per iso; >50% of
  eliminateAll time was in that rebuild).
- Batch-find iso collection in `eliminateAll` (was rescanning from
  zero after every iso).
- `nnz_hint` parameter on `CobComplex.__init__` — lets `AddPosCrossing`,
  `AddNegCrossing`, `AddCap`, `AddCup` skip the full-matrix cell scan
  at construction time.
- Memoization on `Cob.components` and a new
  `Cob._decos_from_old_comp` with shared-reference cache.
- Dict-based aggregation in `Cob.simplify_decos` replacing sort+groupby.
- Cython extension `KhT/cob_ext.pyx` (compiled via
  `KhT/setup_cob_ext.py`) for the inner loop of `Cob.mor.__mul__`.
  Falls back to pure Python if not built.
- Process pool scaffolding for parallel rank-1 update in
  `eliminateAll` (env-var opt-in: `KHT_PAR_THRESHOLD`, `KHT_WORKERS`,
  `KHT_CHUNK`).  In practice the work units are too small for
  multiprocessing to win cleanly at the tangle sizes we've tested;
  kept as scaffolding for future investigation.

End-to-end on TH-pattern: full 41-slice computation in ~15 minutes
(vs. estimates of much longer on unoptimized code).  Per-slice wall
time still scales as roughly `t ∝ N^1.71` — consistent with the
topological density of the tangle, not a fixable algorithmic issue.
