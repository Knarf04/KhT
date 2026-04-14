# Tangle strings below were originally written in the pre-May-2020 right-to-left
# slice convention; ported via Tangle.from_legacy (see CrossingTangle.py).
#
# Known pre-existing library issue (commit d2019ce, Sept 2019): the cabled
# forms of 5_1, 5_2, 6_1 fail the d^2 = 0 check over F_3 and F_5.  They do
# compute cleanly over F_2, so this example uses field = 2.

field = 2

# name1 = "CableKnots/trefoil"
# Tangle1 = Tangle.from_legacy("cup1.neg0.pos1.neg0.cap1", topends=1, botends=1).Cable()
# BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
# multicurve1 = BNcx1.to_multicurve()
# multicurve1.drawpng(name1+"_BNr","hdelta",Tangle1.slices)

# name2 = "CableKnots/fig8"
# Tangle2 = Tangle.from_legacy("cup1.neg0.neg0.pos1.neg0.cap1", topends=1, botends=1).Cable()
# BNcx2 = Tangle2.toReduced_BNComplex(1000, field, "safe")
# multicurve2 = BNcx2.to_multicurve()
# multicurve2.drawpng(name2+"_BNr","hdelta",Tangle2.slices)

name3 = "CableKnots/5_1"
Tangle3 = Tangle.from_legacy("cup1.pos0.neg1.neg1.neg1.pos0.cap1", topends=1, botends=1).Cable()
BNcx3 = Tangle3.toReduced_BNComplex(1000, field)
multicurve3 = BNcx3.to_multicurve()
multicurve3.drawpng(name3+"_BNr","hdelta",Tangle3.slices)

name4 = "CableKnots/5_2"
Tangle4 = Tangle.from_legacy("cup1.neg0.pos1.pos1.neg0.neg0.cap1", topends=1, botends=1).Cable()
BNcx4 = Tangle4.toReduced_BNComplex(1000, field)
multicurve4 = BNcx4.to_multicurve()
multicurve4.drawpng(name4+"_BNr","hdelta",Tangle4.slices)

name5 = "CableKnots/6_1"
Tangle5 = Tangle.from_legacy("cup1.neg0.pos1.pos1.pos1.neg0.neg0.cap1", topends=1, botends=1).Cable()
BNcx5 = Tangle5.toReduced_BNComplex(1000, field)
multicurve5 = BNcx5.to_multicurve()
multicurve5.drawpng(name5+"_BNr","hdelta",Tangle5.slices)