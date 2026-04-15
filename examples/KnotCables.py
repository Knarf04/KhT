# Tangle strings below were originally written in the pre-May-2020 right-to-left
# slice convention; ported via Tangle.from_legacy (see CrossingTangle.py).
#
# Historical note (commit d2019ce, Sept 2019): these cabled knots were
# reported to fail d^2 = 0 over F_3 / F_5 / F_7.  Re-running the check
# in April 2026 with the accumulated F_p arithmetic fixes shows 5_1
# passing validate() over F_2, F_3, F_5, F_7 under both the legacy and
# the experimental "corrected" sign convention (see
# BNAlgebra.set_sign_convention and examples/d2_bisect.py).
# The default here stays at F_2 for reproducibility of prior runs;
# change to 3, 5, or 7 to explore — all currently pass on 5_1.

field = 2

# name1 = "CableKnots/trefoil"
# Tangle1 = Tangle.from_legacy("cup1.neg0.pos1.neg0.cap1", topends=1, botends=1).Cable()
# BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
# multicurve1 = BNcx1.to_multicurve()
# multicurve1.save(name1 + "_BNr")

# name2 = "CableKnots/fig8"
# Tangle2 = Tangle.from_legacy("cup1.neg0.neg0.pos1.neg0.cap1", topends=1, botends=1).Cable()
# BNcx2 = Tangle2.toReduced_BNComplex(1000, field, "safe")
# multicurve2 = BNcx2.to_multicurve()
# multicurve2.save(name2 + "_BNr")

name3 = "CableKnots/5_1"
Tangle3 = Tangle.from_legacy("cup1.pos0.neg1.neg1.neg1.pos0.cap1", topends=1, botends=1).Cable()
BNcx3 = Tangle3.toReduced_BNComplex(1000, field)
multicurve3 = BNcx3.to_multicurve()
multicurve3.save(name3 + "_BNr")

name4 = "CableKnots/5_2"
Tangle4 = Tangle.from_legacy("cup1.neg0.pos1.pos1.neg0.neg0.cap1", topends=1, botends=1).Cable()
BNcx4 = Tangle4.toReduced_BNComplex(1000, field)
multicurve4 = BNcx4.to_multicurve()
multicurve4.save(name4 + "_BNr")

name5 = "CableKnots/6_1"
Tangle5 = Tangle.from_legacy("cup1.neg0.pos1.pos1.pos1.neg0.neg0.cap1", topends=1, botends=1).Cable()
BNcx5 = Tangle5.toReduced_BNComplex(1000, field)
multicurve5 = BNcx5.to_multicurve()
multicurve5.save(name5 + "_BNr")