# Pre-existing library issue (commit d2019ce, Sept 2019): the +3 and +4
# twist quotients of (2,-3)-pretzel tangle fail d^2 = 0 over F_3 and F_5.
# They compute cleanly over F_2, so this example uses field = 2.
field = 2

name1 = "Quotients_of_(2,-3)_PretzelTangle/NoTwists"
Tangle1 = Tangle.quotient_of_2_m3_pretzel_tangle(0)
BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
multicurve1 = BNcx1.to_multicurve()
multicurve1.save(name1 + "_BNr")

name2 = "Quotients_of_(2,-3)_PretzelTangle/+1Twist"
Tangle2 = Tangle.quotient_of_2_m3_pretzel_tangle(1)
BNcx2 = Tangle2.toReduced_BNComplex(1000, field)
multicurve2 = BNcx2.to_multicurve()
multicurve2.save(name2 + "_BNr")

name3 = "Quotients_of_(2,-3)_PretzelTangle/+2Twists"
Tangle3 = Tangle.quotient_of_2_m3_pretzel_tangle(2)
BNcx3 = Tangle3.toReduced_BNComplex(1000, field)
multicurve3 = BNcx3.to_multicurve()
multicurve3.save(name3 + "_BNr")

name4 = "Quotients_of_(2,-3)_PretzelTangle/+3Twists"
Tangle4 = Tangle.quotient_of_2_m3_pretzel_tangle(3)
BNcx4 = Tangle4.toReduced_BNComplex(1000, field)
multicurve4 = BNcx4.to_multicurve()
multicurve4.save(name4 + "_BNr")

name5 = "Quotients_of_(2,-3)_PretzelTangle/+4Twists"
Tangle5 = Tangle.quotient_of_2_m3_pretzel_tangle(4)
BNcx5 = Tangle5.toReduced_BNComplex(1000, field)
multicurve5 = BNcx5.to_multicurve()
multicurve5.save(name5 + "_BNr")

name6 = "Quotients_of_(2,-3)_PretzelTangle/-1Twist"
Tangle6 = Tangle.quotient_of_2_m3_pretzel_tangle(-1)
BNcx6 = Tangle6.toReduced_BNComplex(1000, field)
multicurve6 = BNcx6.to_multicurve()
multicurve6.save(name6 + "_BNr")

name7 = "Quotients_of_(2,-3)_PretzelTangle/-2Twists"
Tangle7 = Tangle.quotient_of_2_m3_pretzel_tangle(-2)
BNcx7 = Tangle7.toReduced_BNComplex(1000, field)
multicurve7 = BNcx7.to_multicurve()
multicurve7.save(name7 + "_BNr")