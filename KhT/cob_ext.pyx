# cython: language_level=3, boundscheck=False, wraparound=False
# Cython acceleration of the hot inner loop of Cob.mor.__mul__.
#
# Ports:
#   _compose_decos_inner: the "for deco1 in decos1: for deco2 in decos2: ..."
#       double loop plus combine_decos + partial-decoration collection.
#
# Keeps Python-level concerns (caching, Cob.mor construction) in Cob.py.
# The decos passed in are standard Python list-of-list-of-int.

from itertools import product as _product


def compose_decos_inner(list decos1,
                        list decos2,
                        list genus,
                        list comps1_x,
                        list comps2_x,
                        list old_comps_x,
                        decos_from_old_comp):
    """Mirror of the inner body of Cob.mor.__mul__.

    Returns a flat list of new decorations (each a Python list:
    [H_power, dot_1, ..., dot_n, coeff]) corresponding to the product
    of ``decos1`` and ``decos2`` under Bar-Natan neck-cutting.

    Intended to be called from Python; strongly typed locals let Cython
    eliminate most of the interpreter overhead in the loop.
    """
    cdef Py_ssize_t n_old = len(genus)
    cdef list new_decos = []
    cdef Py_ssize_t i_d1, i_d2, oi, idx, r
    cdef int hpow, combined_h
    cdef long coeff, coeff1
    cdef list deco1, deco2, c1_idx, c2_idx, oc_idx
    cdef list partial_decos
    cdef list lengths = [len(oc) for oc in old_comps_x]
    cdef Py_ssize_t n_d1 = len(decos1)
    cdef Py_ssize_t n_d2 = len(decos2)

    for i_d1 in range(n_d1):
        deco1 = <list>decos1[i_d1]
        hpow = <int>deco1[0]
        coeff1 = <long>deco1[len(deco1) - 1]
        for i_d2 in range(n_d2):
            deco2 = <list>decos2[i_d2]
            combined_h = hpow + <int>deco2[0]
            coeff = coeff1 * <long>deco2[len(deco2) - 1]

            partial_decos = []
            for oi in range(n_old):
                c1_idx = <list>comps1_x[oi]
                c2_idx = <list>comps2_x[oi]
                r = 0
                for idx in c1_idx:
                    r += <int>deco1[idx + 1]
                for idx in c2_idx:
                    r += <int>deco2[idx + 1]
                partial_decos.append(decos_from_old_comp(<int>genus[oi], r, lengths[oi]))

            # Combine decorations across the cartesian product of partials.
            for combo in _product(*partial_decos):
                new_h = combined_h
                new_coeff = coeff
                middle = []
                for sub in combo:
                    new_h += <int>(<list>sub)[0]
                    new_coeff *= <long>(<list>sub)[len(<list>sub) - 1]
                    middle.extend((<list>sub)[1:len(<list>sub) - 1])
                new_decos.append([new_h] + middle + [new_coeff])

    return new_decos
