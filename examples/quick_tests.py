# Fast regression suite for KhT.  Runs in well under 5 minutes.
#
# Usage:   python3 KhT quick_tests
#
# The suite computes a handful of tangle invariants that mirror the slower
# examples (the 2-cable trefoil from the README, a pretzel tangle, a scaled
# torus knot, a scaled 3-braid) and compares four deterministic structural
# invariants against values captured from a known-good run:
#
#     * validate()        - d^2 = 0 check passes without raising
#     * is_looptype()     - True after clean-up, so the complex is a multicurve
#     * len(mc.comps)     - number of loop-type components
#     * sorted list of    - (q, h, delta) generator gradings after clean-up
#       generator gradings  (order is random due to clean_up's choice(remaining),
#                            but the sorted multiset is an invariant)
#
# Baseline captured 2026-04-14 against the HEAD of the KhT repository.
# All computations are over F_2 unless noted.  Observed wall time: ~1.5 s.
#
# If a refactor changes any frozen value below, either:
#   (a) you broke correctness - investigate before proceeding, or
#   (b) the change is mathematically equivalent but renormalizes something -
#       re-capture the baseline and document why in a commit message.


from time import time

# BNbracket, Tangle, BNComplex etc. are already in scope via the star-imports
# performed by KhT/__main__.py before exec()'ing this file.

# Using a class to hold shared state:  KhT/__main__.py executes this file
# via exec() inside a function, so module-level variables are not visible to
# helper functions defined alongside them.  A class namespace sidesteps that.


class Suite:
    failures = []
    t_start = time()

    @staticmethod
    def grad(cx):
        return sorted((g.q, g.h, g.delta) for g in cx.gens)

    @staticmethod
    def compsizes(cx):
        return sorted(len(c.gens) for c in cx.to_multicurve().comps)

    @staticmethod
    def validate_ok(cx):
        cx.validate()   # raises on failure, returns None on success
        return True

    @classmethod
    def check(cls, label, actual, expected):
        ok = actual == expected
        print("[{}] {}".format("PASS" if ok else "FAIL", label))
        if not ok:
            print("       expected: {}".format(expected))
            print("       actual:   {}".format(actual))
            cls.failures.append(label)

    @classmethod
    def check_truthy(cls, label, cond):
        ok = bool(cond)
        print("[{}] {}".format("PASS" if ok else "FAIL", label))
        if not ok:
            cls.failures.append(label)


check = Suite.check
check_truthy = Suite.check_truthy
_grad = Suite.grad
_compsizes = Suite.compsizes
validate_ok = Suite.validate_ok


# ---------- 1) 2-cable-trefoil (the README demo), F_2 ---------------------
#     Covers:  BNbracket -> ToBNAlgebra(2) -> eliminateAll -> clean_up
#              plus cone(1) + clean_up for the figure-8 invariant

t1 = time()
tangle_2cable = ("cap1.cap2.cap3.neg1.neg2.neg0.neg1.pos3.pos2.pos4.pos3."
                 "neg1.neg0.neg2.neg1.cup3.cup2")

cx = BNbracket(tangle_2cable, 0, 0, 1)
BNr = cx.ToBNAlgebra(2)
BNr.eliminateAll()
BNr.clean_up()

check_truthy("2cable_trefoil_BNr_F2 validate",   validate_ok(BNr))
check_truthy("2cable_trefoil_BNr_F2 looptype",   BNr.is_looptype())
check("2cable_trefoil_BNr_F2 ngens",             len(BNr.gens), 27)
check("2cable_trefoil_BNr_F2 compsizes",         _compsizes(BNr), [12, 15])
check("2cable_trefoil_BNr_F2 gradings",          _grad(BNr), [
    (-2, 0, -1.0), (-1, 1, -1.5), (1, 2, -1.5), (3, 3, -1.5),
    (4, 4, -2.0), (5, 4, -1.5), (5, 5, -2.5), (6, 5, -2.0),
    (6, 6, -3.0), (7, 6, -2.5), (7, 6, -2.5), (7, 7, -3.5),
    (8, 7, -3.0), (9, 7, -2.5), (9, 7, -2.5), (9, 8, -3.5),
    (9, 8, -3.5), (10, 8, -3.0), (11, 8, -2.5), (11, 9, -3.5),
    (11, 9, -3.5), (12, 9, -3.0), (13, 10, -3.5), (13, 10, -3.5),
    (14, 11, -4.0), (15, 11, -3.5), (16, 12, -4.0)])

Khr = BNr.cone(1)
Khr.clean_up()
check_truthy("2cable_trefoil_Khr_F2 validate",   validate_ok(Khr))
check_truthy("2cable_trefoil_Khr_F2 looptype",   Khr.is_looptype())
check("2cable_trefoil_Khr_F2 ngens",             len(Khr.gens), 54)
check("2cable_trefoil_Khr_F2 compsizes",         _compsizes(Khr), [6, 12, 12, 12, 12])
check("2cable_trefoil_Khr_F2 gradings",          _grad(Khr), [
    (-2, 0, -1.0), (-1, 1, -1.5), (0, 1, -1.0), (1, 2, -1.5),
    (1, 2, -1.5), (3, 3, -1.5), (3, 3, -1.5), (4, 4, -2.0),
    (5, 4, -1.5), (5, 4, -1.5), (5, 5, -2.5), (6, 5, -2.0),
    (6, 5, -2.0), (6, 6, -3.0), (7, 5, -1.5), (7, 6, -2.5),
    (7, 6, -2.5), (7, 6, -2.5), (7, 7, -3.5), (8, 6, -2.0),
    (8, 7, -3.0), (8, 7, -3.0), (9, 7, -2.5), (9, 7, -2.5),
    (9, 7, -2.5), (9, 7, -2.5), (9, 8, -3.5), (9, 8, -3.5),
    (9, 8, -3.5), (10, 8, -3.0), (10, 8, -3.0), (11, 8, -2.5),
    (11, 8, -2.5), (11, 8, -2.5), (11, 9, -3.5), (11, 9, -3.5),
    (11, 9, -3.5), (11, 9, -3.5), (12, 9, -3.0), (12, 9, -3.0),
    (13, 9, -2.5), (13, 10, -3.5), (13, 10, -3.5), (13, 10, -3.5),
    (13, 10, -3.5), (14, 10, -3.0), (14, 11, -4.0), (15, 11, -3.5),
    (15, 11, -3.5), (15, 11, -3.5), (16, 12, -4.0), (16, 12, -4.0),
    (17, 12, -3.5), (18, 13, -4.0)])

print("#  2cable_trefoil_F2 wall={:.2f}s".format(time() - t1))


# ---------- 2) 2-cable-trefoil over the default field (F_7) ---------------
#     Same tangle, ToBNAlgebra() with its default field.  Exercises the
#     default-arg codepath.

t2 = time()
cx = BNbracket(tangle_2cable, 0, 0, 1)
BNr_default = cx.ToBNAlgebra()
BNr_default.eliminateAll()
BNr_default.clean_up()
check_truthy("2cable_trefoil_BNr_default validate", validate_ok(BNr_default))
check_truthy("2cable_trefoil_BNr_default looptype", BNr_default.is_looptype())
check("2cable_trefoil_BNr_default ngens",           len(BNr_default.gens), 27)
check("2cable_trefoil_BNr_default compsizes",       _compsizes(BNr_default), [12, 15])
# Gradings come from the cobordism structure and are independent of the field,
# so the default-field result must match the F_2 result exactly.
check("2cable_trefoil_BNr_default gradings == F_2", _grad(BNr_default), _grad(BNr))
print("#  2cable_trefoil_default wall={:.2f}s".format(time() - t2))


# ---------- 3) Pretzel tangle (2, -3), F_2 --------------------------------
#     Covers:  Tangle.PretzelTangle() -> toReduced_BNComplex -> clean_up
#              plus cone(1).  Same pattern as examples/PretzelTangles.py.

t3 = time()
T = Tangle.PretzelTangle(2, -3)
BNcx = T.toReduced_BNComplex(10000, 2)
BNcx.clean_up(10000)
check_truthy("Pretzel(2,-3)_BNr_F2 validate",    validate_ok(BNcx))
check_truthy("Pretzel(2,-3)_BNr_F2 looptype",    BNcx.is_looptype())
check("Pretzel(2,-3)_BNr_F2 ngens",              len(BNcx.gens), 9)
check("Pretzel(2,-3)_BNr_F2 compsizes",          _compsizes(BNcx), [9])
check("Pretzel(2,-3)_BNr_F2 gradings",           _grad(BNcx), [
    (-2, 0, -1.0), (-1, 1, -1.5), (0, 1, -1.0), (1, 2, -1.5),
    (1, 2, -1.5), (3, 3, -1.5), (3, 3, -1.5), (5, 4, -1.5),
    (6, 5, -2.0)])

Khr = BNcx.cone(1)
Khr.clean_up(1000)
check_truthy("Pretzel(2,-3)_Khr_F2 validate",    validate_ok(Khr))
check_truthy("Pretzel(2,-3)_Khr_F2 looptype",    Khr.is_looptype())
check("Pretzel(2,-3)_Khr_F2 ngens",              len(Khr.gens), 18)
check("Pretzel(2,-3)_Khr_F2 compsizes",          _compsizes(Khr), [6, 12])
check("Pretzel(2,-3)_Khr_F2 gradings",           _grad(Khr), [
    (-2, 0, -1.0), (-1, 1, -1.5), (0, 1, -1.0), (0, 1, -1.0),
    (1, 2, -1.5), (1, 2, -1.5), (1, 2, -1.5), (2, 2, -1.0),
    (3, 3, -1.5), (3, 3, -1.5), (3, 3, -1.5), (3, 3, -1.5),
    (5, 4, -1.5), (5, 4, -1.5), (5, 4, -1.5), (6, 5, -2.0),
    (7, 5, -1.5), (8, 6, -2.0)])
print("#  Pretzel(2,-3) wall={:.2f}s".format(time() - t3))


# ---------- 4) Pretzel tangle (2, -2), F_2 --------------------------------
#     Smallest member of the UnlinkTangles family.

t4 = time()
T = Tangle.PretzelTangle(2, -2)
BNcx = T.toReduced_BNComplex(10000, 2)
BNcx.clean_up(10000)
check_truthy("Pretzel(2,-2)_BNr_F2 validate",    validate_ok(BNcx))
check_truthy("Pretzel(2,-2)_BNr_F2 looptype",    BNcx.is_looptype())
check("Pretzel(2,-2)_BNr_F2 ngens",              len(BNcx.gens), 6)
check("Pretzel(2,-2)_BNr_F2 compsizes",          _compsizes(BNcx), [3, 3])
check("Pretzel(2,-2)_BNr_F2 gradings",           _grad(BNcx), [
    (-1, 0, -0.5), (0, 1, -1.0), (2, 2, -1.0),
    (2, 2, -1.0), (4, 3, -1.0), (5, 4, -1.5)])
print("#  Pretzel(2,-2) wall={:.2f}s".format(time() - t4))


# ---------- 5) Scaled torus knot, F_2 -------------------------------------
#     examples/TorusKnot.py uses 6 repetitions of neg2.neg1.neg0 (18 negs).
#     We use 3 repetitions; the full version runs in ~0.9 s but this keeps
#     the suite well under budget on slower machines.

t5 = time()
torus_str = ("cap1.cap2.cap3." + ("neg2.neg1.neg0." * 3) + "cup3.cup2")
T = Tangle(torus_str)
BNcx = T.toReduced_BNComplex(10000, 2)
BNcx.clean_up(10000)
check_truthy("TorusKnot_scaled_BNr_F2 validate", validate_ok(BNcx))
check_truthy("TorusKnot_scaled_BNr_F2 looptype", BNcx.is_looptype())
check("TorusKnot_scaled_BNr_F2 ngens",           len(BNcx.gens), 8)
check("TorusKnot_scaled_BNr_F2 compsizes",       _compsizes(BNcx), [8])
check("TorusKnot_scaled_BNr_F2 gradings",        _grad(BNcx), [
    (3, 4, -2.5), (4, 4, -2.0), (4, 5, -3.0), (5, 5, -2.5),
    (6, 6, -3.0), (8, 7, -3.0), (10, 8, -3.0), (11, 9, -3.5)])
print("#  TorusKnot_scaled wall={:.2f}s".format(time() - t5))


# ---------- 6) Scaled 3-braid, F_2 ----------------------------------------
#     examples/3braid.py uses (neg2 * 21) + (pos0 * 20) and takes ~100s.
#     We cut both arms to 8; the scaled variant runs in under 1 s and still
#     exercises the same compose/eliminate/clean-up machinery on a
#     mid-sized (~78-gen) complex.

t6 = time()
braid_str = ("cap1.cap3." + ("neg2." * 8) + ("pos0." * 8) + "cup1")
T = Tangle(braid_str)
BNcx = T.toReduced_BNComplex(10000, 2)
BNcx.clean_up(10000)
check_truthy("3braid_scaled_F2 validate",        validate_ok(BNcx))
check_truthy("3braid_scaled_F2 looptype",        BNcx.is_looptype())
check("3braid_scaled_F2 ngens",                  len(BNcx.gens), 78)
check("3braid_scaled_F2 compsizes",              _compsizes(BNcx), [9, 9, 12, 20, 28])
check("3braid_scaled_F2 gradings", _grad(BNcx), [
    (-7, 0, -3.5), (-6, 1, -4.0), (-5, 1, -3.5), (-4, 2, -4.0),
    (-4, 2, -4.0), (-3, 2, -3.5), (-2, 3, -4.0), (-2, 3, -4.0),
    (-2, 3, -4.0), (-1, 3, -3.5), (0, 4, -4.0), (0, 4, -4.0),
    (0, 4, -4.0), (0, 4, -4.0), (1, 4, -3.5), (2, 5, -4.0),
    (2, 5, -4.0), (2, 5, -4.0), (2, 5, -4.0), (2, 5, -4.0),
    (3, 5, -3.5), (4, 6, -4.0), (4, 6, -4.0), (4, 6, -4.0),
    (4, 6, -4.0), (4, 6, -4.0), (4, 6, -4.0), (5, 6, -3.5),
    (6, 7, -4.0), (6, 7, -4.0), (6, 7, -4.0), (6, 7, -4.0),
    (6, 7, -4.0), (6, 7, -4.0), (6, 7, -4.0), (8, 8, -4.0),
    (8, 8, -4.0), (8, 8, -4.0), (8, 8, -4.0), (8, 8, -4.0),
    (8, 8, -4.0), (8, 8, -4.0), (8, 8, -4.0), (10, 9, -4.0),
    (10, 9, -4.0), (10, 9, -4.0), (10, 9, -4.0), (10, 9, -4.0),
    (10, 9, -4.0), (10, 9, -4.0), (11, 10, -4.5), (12, 10, -4.0),
    (12, 10, -4.0), (12, 10, -4.0), (12, 10, -4.0), (12, 10, -4.0),
    (12, 10, -4.0), (13, 11, -4.5), (14, 11, -4.0), (14, 11, -4.0),
    (14, 11, -4.0), (14, 11, -4.0), (14, 11, -4.0), (15, 12, -4.5),
    (16, 12, -4.0), (16, 12, -4.0), (16, 12, -4.0), (16, 12, -4.0),
    (17, 13, -4.5), (18, 13, -4.0), (18, 13, -4.0), (18, 13, -4.0),
    (19, 14, -4.5), (20, 14, -4.0), (20, 14, -4.0), (21, 15, -4.5),
    (22, 15, -4.0), (23, 16, -4.5)])
print("#  3braid_scaled wall={:.2f}s".format(time() - t6))


# ---------- 7) two_twist_hitch(0), F_5 -------------------------------------
#     Covers Hitches.py and the two_twist_hitch classmethod, which builds its
#     string in the pre-May-2020 right-to-left convention and reverses it.

t7 = time()
T = Tangle.two_twist_hitch(0)
BNcx = T.toReduced_BNComplex(1000, 5)
BNcx.clean_up(1000)
check_truthy("two_twist_hitch(0)_F5 validate",   validate_ok(BNcx))
check_truthy("two_twist_hitch(0)_F5 looptype",   BNcx.is_looptype())
check("two_twist_hitch(0)_F5 ngens",             len(BNcx.gens), 37)
check("two_twist_hitch(0)_F5 compsizes",         _compsizes(BNcx), [1, 12, 12, 12])
check("two_twist_hitch(0)_F5 gradings",          _grad(BNcx), [
    (-3, 0, -1.5), (-2, 1, -2.0), (-1, 1, -1.5), (-1, 1, -1.5),
    (0, 2, -2.0), (0, 2, -2.0), (0, 2, -2.0), (1, 2, -1.5),
    (1, 2, -1.5), (2, 3, -2.0), (2, 3, -2.0), (2, 3, -2.0),
    (2, 3, -2.0), (2, 3, -2.0), (3, 3, -1.5), (4, 4, -2.0),
    (4, 4, -2.0), (4, 4, -2.0), (4, 4, -2.0), (4, 4, -2.0),
    (4, 4, -2.0), (4, 4, -2.0), (5, 5, -2.5), (6, 5, -2.0),
    (6, 5, -2.0), (6, 5, -2.0), (6, 5, -2.0), (6, 5, -2.0),
    (7, 6, -2.5), (7, 6, -2.5), (8, 6, -2.0), (8, 6, -2.0),
    (8, 6, -2.0), (9, 7, -2.5), (9, 7, -2.5), (10, 7, -2.0),
    (11, 8, -2.5)])
print("#  two_twist_hitch(0)_F5 wall={:.2f}s".format(time() - t7))


# ---------- 8) LiamsTangle(1, [0]), F_2 ------------------------------------
#     Covers Computations.py and the LiamsTangle classmethod.

t8 = time()
T = Tangle.LiamsTangle(1, [0])
BNcx = T.toReduced_BNComplex(1000, 2)
BNcx.clean_up(1000)
check_truthy("LiamsTangle(1,[0])_F2 validate",   validate_ok(BNcx))
check_truthy("LiamsTangle(1,[0])_F2 looptype",   BNcx.is_looptype())
check("LiamsTangle(1,[0])_F2 ngens",             len(BNcx.gens), 17)
check("LiamsTangle(1,[0])_F2 compsizes",         _compsizes(BNcx), [17])
check("LiamsTangle(1,[0])_F2 gradings",          _grad(BNcx), [
    (-4, 0, -2.0), (-3, 1, -2.5), (-2, 1, -2.0), (-1, 2, -2.5),
    (-1, 2, -2.5), (1, 3, -2.5), (1, 3, -2.5), (2, 4, -3.0),
    (3, 4, -2.5), (3, 4, -2.5), (3, 5, -3.5), (4, 5, -3.0),
    (5, 5, -2.5), (5, 6, -3.5), (7, 7, -3.5), (9, 8, -3.5),
    (10, 9, -4.0)])
print("#  LiamsTangle(1,[0])_F2 wall={:.2f}s".format(time() - t8))


# ---------- 9) quotient_of_2_m3_pretzel_tangle(0), F_5 ---------------------
#     Covers Quotients.py and the quotient_of_2_m3_pretzel_tangle classmethod.

t9 = time()
T = Tangle.quotient_of_2_m3_pretzel_tangle(0)
BNcx = T.toReduced_BNComplex(1000, 5)
BNcx.clean_up(1000)
check_truthy("quotient(0)_F5 validate",          validate_ok(BNcx))
check_truthy("quotient(0)_F5 looptype",          BNcx.is_looptype())
check("quotient(0)_F5 ngens",                    len(BNcx.gens), 15)
check("quotient(0)_F5 compsizes",                _compsizes(BNcx), [3, 12])
check("quotient(0)_F5 gradings",                 _grad(BNcx), [
    (5, 6, -3.5), (6, 6, -3.0), (6, 7, -4.0), (7, 7, -3.5),
    (8, 7, -3.0), (8, 8, -4.0), (8, 8, -4.0), (9, 8, -3.5),
    (10, 9, -4.0), (10, 9, -4.0), (12, 10, -4.0), (12, 10, -4.0),
    (13, 11, -4.5), (14, 11, -4.0), (15, 12, -4.5)])
print("#  quotient(0)_F5 wall={:.2f}s".format(time() - t9))


# ---------- 10) from_legacy on a Computations2 tangle, F_2 -----------------
#     Covers Computations2.py and the from_legacy classmethod path.

t10 = time()
T = Tangle.from_legacy("cup0.pos1.neg2.pos1.neg0.neg1.pos2.cap3.cap1")
BNcx = T.toReduced_BNComplex(1000, 2)
BNcx.clean_up(1000)
check_truthy("trefoil_linked_F2 validate",       validate_ok(BNcx))
check_truthy("trefoil_linked_F2 looptype",       BNcx.is_looptype())
check("trefoil_linked_F2 ngens",                 len(BNcx.gens), 10)
check("trefoil_linked_F2 compsizes",             _compsizes(BNcx), [10])
check("trefoil_linked_F2 gradings",              _grad(BNcx), [
    (-2, 0, -1.0), (-1, 1, -1.5), (0, 1, -1.0), (1, 2, -1.5),
    (1, 2, -1.5), (3, 3, -1.5), (3, 3, -1.5), (5, 4, -1.5),
    (5, 4, -1.5), (6, 5, -2.0)])
print("#  trefoil_linked_F2 wall={:.2f}s".format(time() - t10))


# ---------- 11) Cable() of from_legacy, F_2 --------------------------------
#     Covers KnotCables.py:  a 1-1 legacy tangle cabled, then reduced.
#     (F_3/F_5/F_7 hit a known pre-existing d^2 != 0 library bug on this
#     tangle; F_2 is clean.)

t11 = time()
T = Tangle.from_legacy("cup1.pos0.neg1.neg1.neg1.pos0.cap1", topends=1, botends=1).Cable()
BNcx = T.toReduced_BNComplex(1000, 2)
BNcx.clean_up(1000)
check_truthy("cabled_5_1_F2 validate",           validate_ok(BNcx))
check_truthy("cabled_5_1_F2 looptype",           BNcx.is_looptype())
check("cabled_5_1_F2 ngens",                     len(BNcx.gens), 79)
check("cabled_5_1_F2 compsizes",                 _compsizes(BNcx), [12, 12, 12, 12, 31])
print("#  cabled_5_1_F2 wall={:.2f}s".format(time() - t11))


# ---------- 12) F_p-aware Cob path (intermediate (1,3) cleanup) ------------
#     Runs the same 2-cable-trefoil and cabled_5_1 computations but with
#     cleanup_field=2 passed through to BNbracket.  This exercises
#     BNComplex.ToCob (previously broken), the F_p simplify_decos branch in
#     Cob, and the intermediate cleanup heuristics.  Final invariants must
#     still match the Z-Cob path.

t12 = time()
tangle_2cable = ("cap1.cap2.cap3.neg1.neg2.neg0.neg1.pos3.pos2.pos4.pos3."
                 "neg1.neg0.neg2.neg1.cup3.cup2")
cx = BNbracket(tangle_2cable, 0, 0, 1, cleanup_field=2)
BNr = cx.ToBNAlgebra(2)
BNr.eliminateAll()
BNr.clean_up()
check_truthy("Fp_path_2cable validate", validate_ok(BNr))
check_truthy("Fp_path_2cable looptype", BNr.is_looptype())
check("Fp_path_2cable ngens",           len(BNr.gens), 27)
check("Fp_path_2cable compsizes",       _compsizes(BNr), [12, 15])
check("Fp_path_2cable gradings",        _grad(BNr), [
    (-2, 0, -1.0), (-1, 1, -1.5), (1, 2, -1.5), (3, 3, -1.5),
    (4, 4, -2.0), (5, 4, -1.5), (5, 5, -2.5), (6, 5, -2.0),
    (6, 6, -3.0), (7, 6, -2.5), (7, 6, -2.5), (7, 7, -3.5),
    (8, 7, -3.0), (9, 7, -2.5), (9, 7, -2.5), (9, 8, -3.5),
    (9, 8, -3.5), (10, 8, -3.0), (11, 8, -2.5), (11, 9, -3.5),
    (11, 9, -3.5), (12, 9, -3.0), (13, 10, -3.5), (13, 10, -3.5),
    (14, 11, -4.0), (15, 11, -3.5), (16, 12, -4.0)])
print("#  Fp_path_2cable wall={:.2f}s".format(time() - t12))

t13 = time()
T = Tangle.from_legacy("cup1.pos0.neg1.neg1.neg1.pos0.cap1", topends=1, botends=1).Cable()
# Same computation as test 11 but with intermediate cleanup enabled.
cx = BNbracket(T.slices, T.pos, T.neg, T.top, cleanup_field=2)
BNcx = cx.ToBNAlgebra(2)
BNcx.eliminateAll()
BNcx.clean_up()
check_truthy("Fp_path_cabled_5_1 validate", validate_ok(BNcx))
check_truthy("Fp_path_cabled_5_1 looptype", BNcx.is_looptype())
check("Fp_path_cabled_5_1 ngens",           len(BNcx.gens), 79)
check("Fp_path_cabled_5_1 compsizes",       _compsizes(BNcx), [12, 12, 12, 12, 31])
print("#  Fp_path_cabled_5_1 wall={:.2f}s".format(time() - t13))


# ---------- 14) F_p-aware Cob.mor primitives (isIsom, __neg__) -------------
#     Directly exercises the field-aware branches of isIsom and __neg__
#     under Cob.set_field(2).  Ensures -1 and 1 collapse correctly mod 2
#     and the isIsom predicate still recognises an identity cobordism.

t14 = time()

def _fp_primitives_tests(Suite=Suite):
    import BNAlgebra as _BNAlg
    import Cob as _Cob
    # Build a trivial (1,3) identity cobordism via the BN-algebra route.
    obj = _BNAlg.obj(0, 0, 0)
    unit = _BNAlg.mor([[0, 1]], 2)
    clt = obj.ToCob()

    prev = _Cob.set_field(2)
    try:
        idm = unit.ToCob(clt, clt)
        Suite.check_truthy("Fp_isIsom identity under F_2", idm.isIsom())

        neg_idm = -idm
        # Under F_2, coeff -1 must lift to 1; __neg__ should reduce mod p.
        Suite.check_truthy("Fp_neg coeff in [0, p)",
                     all(0 <= d[-1] < 2 for d in neg_idm.decos))
        Suite.check_truthy("Fp_isIsom neg identity under F_2", neg_idm.isIsom())

        dneg = -(-idm)
        Suite.check_truthy("Fp_double_neg identity under F_2", dneg.isIsom())
    finally:
        _Cob.set_field(prev)

_fp_primitives_tests()
print("#  Fp_primitives wall={:.2f}s".format(time() - t14))


# ---------- 15) Deterministic clean_up_once ---------------------------------
#     Exercises the deterministic=True path through clean_up on the
#     2-cable-trefoil.  The reduction-measure choice is not proven to
#     converge to loop-type (OPEN_QUESTIONS.md item 8, residual), so we
#     only assert: (a) the algorithm does not stall / crash, (b) d^2 = 0
#     is preserved, (c) generator count matches the randomised path, so
#     the deterministic variant at least gets eliminateAll-level
#     reduction for free.

t15 = time()
cx = BNbracket(tangle_2cable, 0, 0, 1)
BNr = cx.ToBNAlgebra(2)
BNr.eliminateAll()
BNr.clean_up(max_iter=50, deterministic=True)
check_truthy("deterministic_2cable validate", validate_ok(BNr))
check("deterministic_2cable ngens",           len(BNr.gens), 27)
print("#  deterministic_2cable wall={:.2f}s".format(time() - t15))


# ---------- 16) signed-lift opt-in path ------------------------------------
#     Runs the 2-cable-trefoil through the F_p intermediate-cleanup path
#     with signed_lift=True.  Final invariants must match the unsigned
#     path (both are reduced mod p downstream, so the set of mod-p
#     classes is the same).  Exercises the new signed-lift branch in
#     BNAlgebra.mor.ToCob + BNComplex.ToCob.

t16 = time()
cx = BNbracket(tangle_2cable, 0, 0, 1, cleanup_field=2, signed_lift=True)
BNr = cx.ToBNAlgebra(2)
BNr.eliminateAll()
BNr.clean_up()
check_truthy("signed_lift_2cable validate", validate_ok(BNr))
check_truthy("signed_lift_2cable looptype", BNr.is_looptype())
check("signed_lift_2cable ngens",           len(BNr.gens), 27)
check("signed_lift_2cable compsizes",       _compsizes(BNr), [12, 15])
print("#  signed_lift_2cable wall={:.2f}s".format(time() - t16))


# ---------- 17) Cob.mor.arrow_length + CobComplex.clean_up scaffold -------
#     Smoke test for the Stage-4 scaffolding from OPEN_QUESTIONS.md
#     item 5.  arrow_length returns min H-power; clean_up currently
#     delegates to eliminateAll and must be idempotent on an already
#     cleaned complex.

t17 = time()
cx = BNbracket(tangle_2cable, 0, 0, 1)
before = len(cx.gens)
cx.clean_up(max_iter=5)
after = len(cx.gens)
check_truthy("cob_clean_up non-increasing", after <= before)
# arrow_length exposed on Cob.mor
import Cob as _Cob_smoke
sample = next((m for row in cx.diff for m in row if m != 0), None)
check_truthy("arrow_length numeric",
             sample is None or isinstance(sample.arrow_length(), (int, float)))
print("#  cob_clean_up wall={:.2f}s".format(time() - t17))


# ---------- Summary --------------------------------------------------------

elapsed = time() - Suite.t_start
print()
print("=" * 60)
if Suite.failures:
    print("FAIL: {} regression(s):".format(len(Suite.failures)))
    for name in Suite.failures:
        print("  - " + name)
    print("=" * 60)
    raise Exception("quick_tests regressed on {} check(s)".format(len(Suite.failures)))
else:
    print("PASS: all checks green ({:.2f}s wall)".format(elapsed))
    print("=" * 60)
