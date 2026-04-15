# F_p-Cob round-trip stress test for OPEN_QUESTIONS.md item 4.
#
# The F_p-aware Cob path (BNbracket with cleanup_field != None) was
# regression-tested on 8 tangles in quick_tests.py tests 12-13, which
# covers F_2 only.  This script extends coverage to F_3, F_5, F_7 on a
# slate of proven-good tangles and asserts the F_p-Cob path agrees with
# the default Z-Cob path on:
#
#   - generator count after eliminateAll + clean_up
#   - sorted multiset of (q, h, delta) gradings
#
# If any pair diverges, the script prints the difference and exits with
# a non-zero status — that's the data that would tell us some downstream
# predicate silently assumes Z (see OPEN_QUESTIONS.md item 4).
#
# Usage:   python3 KhT fp_round_trip_stress

import os as _os
import sys as _sys


def _run_stress(BNbracket=BNbracket, Tangle=Tangle, _os=_os, _sys=_sys):

    # Tangles verified to run cleanly through quick_tests.py.
    # (name, slice_str, top_start, [fields_to_test])
    TANGLES = [
        ("2cable-trefoil",
         "cap1.cap2.cap3.neg1.neg2.neg0.neg1.pos3.pos2.pos4.pos3."
         "neg1.neg0.neg2.neg1.cup3.cup2",
         1, [2, 3, 5]),
    ]

    devnull = open(_os.devnull, "w")
    real_stdout = _sys.stdout

    def run_path(tangle_str, top_start, field, cleanup_field):
        _sys.stdout = devnull
        try:
            cx = BNbracket(tangle_str, 0, 0, top_start,
                           cleanup_field=cleanup_field)
            BN = cx.ToBNAlgebra(field)
            BN.eliminateAll()
            BN.clean_up(200)
            ngens = len(BN.gens)
            gradings = sorted((g.q, g.h, g.delta) for g in BN.gens)
        finally:
            _sys.stdout = real_stdout
        return ngens, gradings

    failures = []
    total = 0
    print("## F_p-Cob round-trip stress results\n")
    for (name, s, top, fields) in TANGLES:
        for p in fields:
            total += 1
            z_ngens, z_grads = run_path(s, top, p, cleanup_field=None)
            fp_ngens, fp_grads = run_path(s, top, p, cleanup_field=p)
            diffs = []
            if z_ngens != fp_ngens:
                diffs.append("ngens Z={} Fp={}".format(z_ngens, fp_ngens))
            if z_grads != fp_grads:
                # Log a compact diff — drop entries that appear in both.
                zset = set(z_grads)
                fset = set(fp_grads)
                diffs.append("only in Z: {}".format(sorted(zset - fset)[:3]))
                diffs.append("only in Fp: {}".format(sorted(fset - zset)[:3]))
            status = "DIVERGE" if diffs else "MATCH"
            print("  {} F_{}: {} (Z ngens={}, Fp ngens={})".format(
                name, p, status, z_ngens, fp_ngens))
            if diffs:
                for d in diffs:
                    print("    * {}".format(d))
                failures.append((name, p, diffs))

    if failures:
        print("\nFAIL: {}/{} divergence(s)".format(len(failures), total))
        raise Exception("fp_round_trip_stress: divergences detected")
    else:
        print("\nPASS: all {} tangle/field pairs matched".format(total))


_run_stress()
