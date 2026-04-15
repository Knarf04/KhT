name = "TH-pattern"
tangle_str = "cap1.cap2.cap3.neg5.neg4.neg1.neg2.pos0.pos1.pos0.pos5.cap3.neg2.cap7.cap8.pos6.pos7.neg5.neg6.cup7.cap11.neg10.neg9.neg8.neg7.neg7.neg8.neg9.neg10.neg1.neg0.neg4.neg3.neg5.neg4.pos1.cup2.cup1.cup2.cup1.cup2"
Tangle = Tangle(tangle_str)

cx = BNbracket(tangle_str, 0, 0, 1)
BNr = cx.ToBNAlgebra(2)
BNr.eliminateAll()
BNr.clean_up()
print("RESULT: gens={} looptype={} comps={}".format(
    len(BNr.gens), BNr.is_looptype(), len(BNr.to_multicurve().comps) if BNr.is_looptype() else "n/a"))

if BNr.is_looptype():
    multicurve = BNr.to_multicurve()
    multicurve.save(name)
    # Produce the PDF + PNG visual output of the multicurve invariant.
    html_row = multicurve.html(name, "_BNr2", "hdelta", Tangle)
    header = """---
title: TH-pattern
layout: default
filename: TH-pattern
---
<link rel="stylesheet" href="css/main.css">
<h4>Arc invariant over \\(\\mathbb{F}_2\\)</h4>
"""
    with open("examples/" + filename + ".html", "w") as f:
        print(header + html_row, file=f)
    print("wrote examples/{}.html plus PDFs in examples/{}/".format(name, name))
