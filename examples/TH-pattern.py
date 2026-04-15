name = "TH-pattern"
tangle = "cap1.cap2.cap3.neg5.neg4.neg1.neg2.pos0.pos1.pos0.pos5.cap3.neg2.cap7.cap8.pos6.pos7.neg5.neg6.cup7.cap11.neg10.neg9.neg8.neg7.neg7.neg8.neg9.neg10.neg1.neg0.neg4.neg3.neg5.neg4.pos1.cup2.cup1.cup2.cup1.cup2"
cx = BNbracket(tangle, 0, 0, 1)
BNr = cx.ToBNAlgebra(2)
BNr.eliminateAll()
BNr.clean_up()
print("RESULT: gens={} looptype={} comps={}".format(
    len(BNr.gens), BNr.is_looptype(), len(BNr.to_multicurve().comps) if BNr.is_looptype() else "n/a"))
