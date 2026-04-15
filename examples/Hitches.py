field = 5

name1 = "TwoTwistHitch/NoCrossings"
Tangle1 = Tangle.two_twist_hitch(0)
BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
multicurve1 = BNcx1.to_multicurve()
multicurve1.save(name1 + "_BNr")

name2 = "TwoTwistHitch/+1Crossing"
Tangle2 = Tangle.two_twist_hitch(1)
BNcx2 = Tangle2.toReduced_BNComplex(1000, field)
multicurve2 = BNcx2.to_multicurve()
multicurve2.save(name2 + "_BNr")

name3 = "TwoTwistHitch/+2Crossing"
Tangle3 = Tangle.two_twist_hitch(2)
BNcx3 = Tangle3.toReduced_BNComplex(1000, field)
multicurve3 = BNcx3.to_multicurve()
multicurve3.save(name3 + "_BNr")

name4 = "TwoTwistHitch/+3Crossing"
Tangle4 = Tangle.two_twist_hitch(3)
BNcx4 = Tangle4.toReduced_BNComplex(1000, field)
multicurve4 = BNcx4.to_multicurve()
multicurve4.save(name4 + "_BNr")

name5 = "TwoTwistHitch/+4Crossing"
Tangle5 = Tangle.two_twist_hitch(4)
BNcx5 = Tangle5.toReduced_BNComplex(1000, field)
multicurve5 = BNcx5.to_multicurve()
multicurve5.save(name5 + "_BNr")

name6 = "TwoTwistHitch/-1Crossing"
Tangle6 = Tangle.two_twist_hitch(-1)
BNcx6 = Tangle6.toReduced_BNComplex(1000, field)
multicurve6 = BNcx6.to_multicurve()
multicurve6.save(name6 + "_BNr")

name7 = "TwoTwistHitch/-2Crossing"
Tangle7 = Tangle.two_twist_hitch(-2)
BNcx7 = Tangle7.toReduced_BNComplex(1000, field)
multicurve7 = BNcx7.to_multicurve()
multicurve7.save(name7 + "_BNr")