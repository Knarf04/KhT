# -*- coding: utf-8 -*-
# COPYRIGHT 2019 Gurkeerat Chhina, Claudius Zibrowius
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd
from tabulate import tabulate
from time import time
import os
import multiprocessing

import Cob
import BNComplexes


# Module-level process pool for parallelizing the rank-1 update in
# eliminateAll.  Created lazily on first use; shared across all
# eliminateAll calls to amortize fork + startup cost.
#
# The pool's workers execute ``_parallel_mul_chunk`` — a simple wrapper
# around Cob.mor.__mul__ that processes a batch of (i, j, iv, ov) tuples
# to amortize per-task pickling overhead.  Threshold below which we fall
# back to the serial path (parallelism overhead would dominate):
#   ``nonzero_out * nonzero_in >= _PARALLEL_THRESHOLD``
_POOL = None
_POOL_SIZE = min(int(os.environ.get("KHT_WORKERS", os.cpu_count() or 1)),
                 os.cpu_count() or 1)
_PARALLEL_THRESHOLD = int(os.environ.get("KHT_PAR_THRESHOLD", "400"))
_CHUNK_SIZE = int(os.environ.get("KHT_CHUNK", "32"))


def _parallel_mul_chunk(chunk):
    """Worker entry point: compute (i, j, iv * ov) for each tuple in chunk.

    Workers run in forked processes, so they inherit the loaded Cob
    module (and its cached tables).  Returns only nonzero deltas to
    keep the pickle-back payload small.
    """
    out = []
    for (i, j, iv, ov) in chunk:
        delta = iv * ov
        if delta != 0:
            out.append((i, j, delta))
    return out


def _get_pool():
    """Lazily create the module-level pool.  Fork start method keeps
    the child processes cheap (no Python re-initialization)."""
    global _POOL
    if _POOL is None and _POOL_SIZE > 1:
        ctx = multiprocessing.get_context("fork")
        _POOL = ctx.Pool(_POOL_SIZE)
    return _POOL

class CobComplex(object):
    """ A chain complex is a directed graph, consisting of 
        - A list of CLTs (Cob.obj) as labels on the vertices
        - A a matrix of cobordisms as the adjacency matrix.
        These should satisfy the usual rules of a chain complex, ie that the differential squared = 0
        Note that the matrix's rows and columns depend on the order the CLTS are given in the list 
        We assume that all entries of 'diff' are reduced in the sense that 'ReduceDecorations()' does not change them.
        """
    __slots__ = 'gens','diff','field','_nnz','_row','_col'

    def __init__(self,gens,diff,field=1,nnz_hint=None):
        """Construct a chain complex.

        ``nnz_hint``, if provided, must be an iterable of ``(target, source)``
        index pairs listing every non-zero position in ``diff``.  When set,
        the full matrix scan is skipped — constructors like
        :func:`AddPosCrossing`, :func:`AddNegCrossing`, :func:`AddCup`, and
        :func:`AddCap` already know exactly which cells they filled, so
        passing a hint shaves O(N^2) per construction on wide tangles.
        """
        self.gens = gens
        self.diff = np.array(diff)
        self.field = field # not implemented; assuming integer coefficients throughout
        self._nnz = set()
        self._row = {}
        self._col = {}
        if nnz_hint is not None:
            for (i, j) in nnz_hint:
                self._nnz.add((i, j))
                self._row.setdefault(i, set()).add(j)
                self._col.setdefault(j, set()).add(i)
        else:
            for i, row in enumerate(self.diff):
                for j, cob in enumerate(row):
                    if cob != 0:
                        self._nnz.add((i, j))
                        self._row.setdefault(i, set()).add(j)
                        self._col.setdefault(j, set()).add(i)

    def _nnz_add(self, i, j):
        if (i, j) in self._nnz:
            return
        self._nnz.add((i, j))
        self._row.setdefault(i, set()).add(j)
        self._col.setdefault(j, set()).add(i)

    def _nnz_discard(self, i, j):
        if (i, j) not in self._nnz:
            return
        self._nnz.discard((i, j))
        r = self._row.get(i)
        if r is not None:
            r.discard(j)
            if not r:
                del self._row[i]
        c = self._col.get(j)
        if c is not None:
            c.discard(i)
            if not c:
                del self._col[j]
    
    def __repr__(self):
        return "CobComplex({},{},{})".format(self.gens,[list(row) for row in self.diff],self.field)
    
    def save(self,filename):
        with open("examples/data/CobComplexes/"+filename, "w") as text_file:
            print(repr(self), file=text_file)
    
    def print(self,switch="long"):
        """Print a complex in human readable form. The optional parameter should be one of the following strings: 
        - 'short' (default) prints only the length of cobordisms.
        - 'long' prints all cobordism data in a nice table.
        - 'old long' prints all cobordism data, but as a list of lists.
        """
        print("The generators:")
        print(pd.DataFrame({\
            "clt.pairs": [clt.pairs for clt in self.gens],\
            "q": [clt.q for clt in self.gens],\
            "h": [clt.h for clt in self.gens]
            },columns=["clt.pairs","q","h"]))
        print("The differential: ("+switch+" form)")
        #print(pd.DataFrame([[print(entry,switch)  for entry in row] for row in self.diff]))
        def prt(entry):
            if entry ==0:
                return ""
            return entry.print(switch)
        print(tabulate(pd.DataFrame([[prt(entry)  for entry in row] for row in self.diff]),range(len(self.diff)),tablefmt="fancy_grid"))
    
    def validate(self): #checks that the differential squares to 0, has no self loops, and is a matrix of the correct size
        length = len(self.gens)
        if len(self.diff) != length:
            raise Exception('Differential does not have n rows (where n is the number of gens in chain complex)')
                
        for i in self.diff:
            if len(i) != length:
                raise Exception('Differential does not have n columns (where n is the number of gens in chain complex)')
        for i in range(length):
            if self.diff[i][i] != 0:
                raise Exception('Differential has self loops')
        
        for i,row in enumerate(self.diff):
            for j,cob in enumerate(row):   
                if cob != 0:
                    if cob.homogeneousQ() == False:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: The component of the differential in row "+str(i)+" and column "+str(j)+" is not homoegenous:")
                        print(cob.print("long"))
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Non-homogeneous morphism in differential!')
                    if (self.gens[i]).h-(self.gens[j]).h!=1:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: The homological grading along the component of the differential in row "+str(i)+" and column "+str(j)+" does not increase by 1:")
                        print(cob.print("long"))
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Something is wrong with the homological grading!')
                    if (self.gens[i]).q-(self.gens[j]).q+cob.deg()!=0:
                        print("!!!!!!!!!!!!!!!!!!")
                        print("ERROR: The quantum grading is not preserved along the component of the differential in row "+str(i)+"(q:"+str((self.gens[i]).q)+") and column "+str(j)+"(q:"+str((self.gens[j]).q)+"):")
                        print(cob.print("long"), " degree:", cob.deg())
                        print("!!!!!!!!!!!!!!!!!!")
                        raise Exception('Something is wrong with the quantum grading!')
        # Computing diff squared
        for source in range(length):
            for target in range(length):
                sum = 0
                for i, row in enumerate(self.diff):
                    sum += row[source] * self.diff[target][i]
                if (sum!=0) and (sum.ReduceDecorations() != []):
                    raise Exception("diff does not square to 0")
        # squared = np.tensordot(self.diff,self.diff, axes=(-2,-1))
        # print("got this far 2")
        # for i,row in enumerate(squared):
            # for j,cob in enumerate(row):
                # if j == 0:
                    # print("got this far 4")
                # if (cob != 0) and (cob.ReduceDecorations() != []):
                    # print("got this far 5")
                    # print("!!!!!!!!!!!!!!!!!!")
                    # print("ERROR: Found non-zero term in d² in row "+str(i)+" and column "+str(j)+":")
                    # print(cob.print("long"))
                    # print("!!!!!!!!!!!!!!!!!!")
                    # raise Exception('Differential does not square to 0')
    
    def findIsom(self):
        """Returns the location of the first isomorphism it finds
           If no isomorphism is found, returns None"""
        # Iterate only the non-zero cells via the incrementally-maintained
        # sparse index.  Set-iteration order is arbitrary, but elimination
        # order does not affect the final invariants (d^2 = 0 is preserved
        # regardless of which iso we cancel first) - we explicitly verify
        # this via the quick_tests regression suite.
        for ti, si in self._nnz:
            if self.diff[ti, si].isIsom():
                return [si, ti]
        return None
    
    def eliminateIsom(self, sourceindex, targetindex):
        """Mutates self by eliminating the isomorphism at the specified location, via the gaussian elimination lemma
           Note that this does not check that the cobodism specified is actually an isomorphism.
        """

        Max=max(targetindex,sourceindex)
        Min=min(targetindex,sourceindex)

        del self.gens[Max] # eliminate source and target from list of generators
        del self.gens[Min] # eliminate source and target from list of generators

        out_source = np.delete(self.diff[:,sourceindex],[Min,Max],0) #arrows starting at the source, omiting indices targetindex and sourceindex
        in_target = np.delete(self.diff[targetindex],[Min,Max],0) #arrows ending at the target, omiting indices targetindex and sourceindex

        flip_sign = (self.diff[targetindex,sourceindex]).decos[0][-1]==1

        # Drop rows + cols at {Min, Max} in one allocation via bool-mask
        # indexing.  Faster than two sequential np.delete calls.
        N = self.diff.shape[0]
        keep_mask = np.ones(N, dtype=bool)
        keep_mask[Min] = False
        keep_mask[Max] = False
        self.diff = self.diff[keep_mask][:, keep_mask]

        # Re-index _nnz / _row / _col: drop positions involving the deleted
        # row/col and shift remaining indices down by 1 or 2 as appropriate.
        def shift(k):
            if k < Min:  return k
            if k < Max:  return k - 1
            return k - 2
        self._nnz = {(shift(i), shift(j)) for (i, j) in self._nnz
                     if i != Min and i != Max and j != Min and j != Max}
        self._row = {shift(t): {shift(s) for s in ss
                                if s != Min and s != Max}
                     for (t, ss) in self._row.items()
                     if t != Min and t != Max}
        self._row = {t: ss for (t, ss) in self._row.items() if ss}
        self._col = {shift(s): {shift(t) for t in ts
                                if t != Min and t != Max}
                     for (s, ts) in self._col.items()
                     if s != Min and s != Max}
        self._col = {s: ts for (s, ts) in self._col.items() if ts}

        # Sparse rank-1 update:  diff[i, j] += in_target[j] * out_source[i].
        # Profiling a cabled-knot tangle shows in_target and out_source are
        # typically ~3 non-zero entries out of ~290, so building the full
        # outer product via np.tensordot and iterating every cell of the
        # matrix-add wastes essentially all of the work.  Iterating only the
        # non-zero pairs drops eliminateIsom from O(N^2) mor-adds to
        # O(nnz_in * nnz_out).
        nonzero_out = [(i, v) for i, v in enumerate(out_source) if v != 0]
        if not nonzero_out:
            return
        nonzero_in  = [(j, -v if flip_sign else v) for j, v in enumerate(in_target) if v != 0]
        if not nonzero_in:
            return
        for i, ov in nonzero_out:
            for j, iv in nonzero_in:
                delta = iv * ov
                if delta == 0:
                    continue
                cur = self.diff[i, j]
                if cur == 0:
                    self.diff[i, j] = delta
                    self._nnz_add(i, j)
                else:
                    new = cur + delta
                    self.diff[i, j] = new
                    if new == 0:
                        self._nnz_discard(i, j)

    def eliminateAll(self): #mutates self by eliminating isomorphisms as long as it can find one
        """Eliminate all rank-1 isomorphism arrows from the diff matrix.

        Performance-critical path.  Two optimizations vs a naive scan:

        1. **Batch-find per pass**: each outer pass scans ``_nnz`` ONCE
           and collects every iso that exists, rather than restarting
           the scan after each iso.  Non-conflicting isos (disjoint
           row/col sets) get applied in sequence; conflicts queue for
           the next pass.  On TH-pattern at width 13, reduces the find
           phase from ~50% of runtime to near zero on later passes.

        2. **Deferred matrix compaction**: row/col deletions are just
           marked in the ``deleted`` set; the actual ``np.delete`` and
           remapping happens once at the end.
        """
        deleted = set()
        while True:
            # One full scan of _nnz: collect every iso position.
            candidates = []
            for (ti, si) in self._nnz:
                if ti in deleted or si in deleted:
                    continue
                if self.diff[ti, si].isIsom():
                    candidates.append((si, ti))
            if not candidates:
                break

            # Apply non-conflicting isos in sequence.  Two isos conflict
            # if they share Min/Max with any already-applied iso in this
            # pass — applying one invalidates the "is iso" status of
            # entries in those rows/cols.  Conflicts get deferred to
            # the next pass.
            applied = 0
            touched = set()  # row/col indices whose arrows were modified this pass
            for (si, ti) in candidates:
                if si in deleted or ti in deleted:
                    continue
                # Conflict check: if this iso's pivot is in a row/col that
                # another iso already modified this pass, re-verify.
                if si in touched or ti in touched:
                    cur_pivot = self.diff[ti, si]
                    if cur_pivot == 0 or not cur_pivot.isIsom():
                        continue  # no longer an iso after earlier updates

                Min = si if si < ti else ti
                Max = si if si > ti else ti
                flip_sign = self.diff[ti, si].decos[0][-1] == 1

                out_source_col = self.diff[:, si]
                in_target_row  = self.diff[ti]
                nonzero_out = [(i, v) for i, v in enumerate(out_source_col)
                               if v != 0 and i != si and i != ti and i not in deleted]
                nonzero_in = [(j, -v if flip_sign else v) for j, v in enumerate(in_target_row)
                              if v != 0 and j != si and j != ti and j not in deleted]

                deleted.add(Min); deleted.add(Max)
                touched.add(Min); touched.add(Max)
                if nonzero_out and nonzero_in:
                    # Track touched indices for conflict-detection on the
                    # next iso in this pass.
                    for i, _ in nonzero_out: touched.add(i)
                    for j, _ in nonzero_in:  touched.add(j)

                    total_muls = len(nonzero_out) * len(nonzero_in)
                    pool = _get_pool() if total_muls >= _PARALLEL_THRESHOLD else None

                    if pool is None:
                        # Serial path (small isos — parallel overhead would dominate).
                        for i, ov in nonzero_out:
                            for j, iv in nonzero_in:
                                delta = iv * ov
                                if delta == 0:
                                    continue
                                cur = self.diff[i, j]
                                if cur == 0:
                                    self.diff[i, j] = delta
                                    self._nnz_add(i, j)
                                else:
                                    new = cur + delta
                                    self.diff[i, j] = new
                                    if new == 0:
                                        self._nnz_discard(i, j)
                    else:
                        # Parallel path: batch mul computations into chunks,
                        # farm out to workers, apply results serially.
                        tasks = [(i, j, iv, ov)
                                 for (i, ov) in nonzero_out
                                 for (j, iv) in nonzero_in]
                        chunks = [tasks[k:k+_CHUNK_SIZE]
                                  for k in range(0, len(tasks), _CHUNK_SIZE)]
                        for chunk_result in pool.imap_unordered(
                                _parallel_mul_chunk, chunks):
                            for (i, j, delta) in chunk_result:
                                cur = self.diff[i, j]
                                if cur == 0:
                                    self.diff[i, j] = delta
                                    self._nnz_add(i, j)
                                else:
                                    new = cur + delta
                                    self.diff[i, j] = new
                                    if new == 0:
                                        self._nnz_discard(i, j)
                # Prune pivot row/col from _nnz.
                for t in (Min, Max):
                    for j in list(self._row.get(t, ())):
                        self._nnz_discard(t, j)
                    for i in list(self._col.get(t, ())):
                        self._nnz_discard(i, t)
                applied += 1

            if applied == 0:
                break

        # Compact: drop all deleted rows/cols in one numpy allocation.
        if deleted:
            N = self.diff.shape[0]
            keep_mask = np.ones(N, dtype=bool)
            for k in deleted:
                keep_mask[k] = False
            self.diff = self.diff[keep_mask][:, keep_mask]
            self.gens = [g for i, g in enumerate(self.gens) if i not in deleted]
            # Remap _nnz / _row / _col indices to the compacted matrix.
            deleted_sorted = sorted(deleted)
            import bisect
            def shift(k):
                return k - bisect.bisect_left(deleted_sorted, k)
            self._nnz = {(shift(i), shift(j)) for (i, j) in self._nnz}
            self._row = {shift(t): {shift(s) for s in ss} for (t, ss) in self._row.items()}
            self._col = {shift(s): {shift(t) for t in ts} for (s, ts) in self._col.items()}
    
    def _nnz_by_col(self):
        """Dict source -> sorted list of targets for nonzero diff[target, source]."""
        d = {}
        for (t, s) in self._nnz:
            d.setdefault(s, []).append(t)
        for v in d.values():
            v.sort()
        return d

    def _nnz_by_row(self):
        """Dict target -> sorted list of sources for nonzero diff[target, source]."""
        d = {}
        for (t, s) in self._nnz:
            d.setdefault(t, []).append(s)
        for v in d.values():
            v.sort()
        return d

    def _nnz_targets_of(self, source):
        """Set of targets t with diff[t, source] != 0."""
        return {t for (t, s) in self._nnz if s == source}

    def _cob_isotopy(self, start, end, alg, col_idx, row_idx, entropy_guard=True):
        """Apply a chain-homotopy-equivalence isotopy along the arrow
        ``start → end`` labeled by Cob.mor ``alg``, iterating only over
        the pre-computed non-zero row/column indices.

        ``col_idx`` maps source -> list of targets (rows nonzero in a column).

        If ``entropy_guard=True`` (default), the method first computes all
        affected entries tentatively, then commits only if the total
        decoration count of the touched entries does not increase.
        This turns the absorption into a *monotone* reduction of deco
        count, preventing the super-linear blowup observed when blindly
        applying every candidate isotopy (see OPEN_QUESTIONS item 5 —
        the BN-algebra version converges without this guard because of
        the (1,3) algebra's finite structure; at wider widths the
        isotopy's "collateral decorations" can explode).

        Returns True if the isotopy was committed, False if skipped.
        """
        alg_neg = -alg
        def _deco_count(entry):
            return 0 if entry == 0 else len(entry.decos)

        # Read iteration sets from self._nnz LIVE (not from the passed
        # col_idx/row_idx), so that earlier isotopies in the same pass
        # don't leave us with stale iteration sets — that would miss
        # updates and break d^2 = 0.
        sources_into_start = [s for (t, s) in self._nnz if t == start]
        targets_from_end = [t for (t, s) in self._nnz if s == end]

        # Phase A analogue: diff[end, s] += diff[start, s] * alg_neg
        # for each source s with an arrow s -> start.
        new_row = {}  # (end, s) -> new mor
        for s in sources_into_start:
            entry = self.diff[start, s]
            if entry == 0:
                continue
            addend = entry * alg_neg
            if addend == 0:
                continue
            cur = self.diff[end, s]
            new_row[(end, s)] = addend if cur == 0 else cur + addend

        # Phase B analogue: diff[t, start] += alg * diff[t, end]
        # for each target t with an arrow end -> t.
        new_col = {}  # (t, start) -> new mor
        for t in targets_from_end:
            entry = self.diff[t, end]
            if entry == 0:
                continue
            addend = alg * entry
            if addend == 0:
                continue
            cur = self.diff[t, start]
            new_col[(t, start)] = addend if cur == 0 else cur + addend

        if entropy_guard:
            before = 0
            after = 0
            for (i, j) in new_row:
                before += _deco_count(self.diff[i, j])
                after += _deco_count(new_row[(i, j)])
            for (i, j) in new_col:
                before += _deco_count(self.diff[i, j])
                after += _deco_count(new_col[(i, j)])
            if after > before:
                return False

        # Commit.
        for (i, j), new in new_row.items():
            self.diff[i, j] = new
            if new == 0:
                self._nnz_discard(i, j)
            else:
                self._nnz_add(i, j)
        for (i, j), new in new_col.items():
            self.diff[i, j] = new
            if new == 0:
                self._nnz_discard(i, j)
            else:
                self._nnz_add(i, j)
        return True

    def _cob_isolate_h_arrow(self, start, end, col_idx, row_idx):
        """Attempt Cob-level arrow-shortening from ``start → end``.

        Generalization of the (1,3) BN-algebra ``isolate_arrow`` to Cob
        at arbitrary width, covering **same-face absorption** including
        dotted decorations (OPEN_QUESTIONS item 5, step A).

        For each dot-pattern ``P`` appearing on a ``±1``-coeff deco of
        ``self.diff[end, start]``, the deco with minimum H-power becomes
        the "canonical short" arrow for that pattern.  Then, for each
        other arrow ``diff[index, start]`` (index ≠ end, with
        ``diff[index, end] == 0`` for isotopy safety), any deco sharing
        dot-pattern ``P`` and H-power ≥ short's gets a tentative
        absorption via a pure-H-power isotopy label.

        Cross-idempotent saddle shortening (dot pattern shifts across
        CLTs) is step B and not yet implemented.  The entropy guard
        inside :meth:`_cob_isotopy` handles the cases where composition
        through neck-cutting would produce more decoration than it
        removes — those isotopies get skipped.

        Returns the number of committed absorptions.
        """
        arrow = self.diff[end, start]
        if arrow == 0:
            return 0
        field = Cob.get_field()
        if field is not None and field > 1:
            minus_one = (-1) % field
        else:
            minus_one = -1

        # Group arrow's ±1-coeff decos by dot pattern; keep the shortest H.
        by_pattern = {}  # tuple(dots) -> (short_H, short_coeff)
        for d in arrow.decos:
            if d[-1] != 1 and d[-1] != minus_one:
                continue
            pat = tuple(d[1:-1])
            if pat not in by_pattern or d[0] < by_pattern[pat][0]:
                by_pattern[pat] = (d[0], d[-1])
        if not by_pattern:
            return 0

        end_clt = self.gens[end]
        absorbed = 0
        targets_of_end = set(col_idx.get(end, ()))

        for index in list(col_idx.get(start, ())):
            if index == end or index == start:
                continue
            if index in targets_of_end:
                continue  # isotopy not safe
            m = self.diff[index, start]
            if m == 0:
                continue
            # Snapshot candidates: decos whose pattern matches some entry
            # in by_pattern, with H-power ≥ the short entry.
            cand = []
            for d in m.decos:
                pat = tuple(d[1:-1])
                if pat not in by_pattern:
                    continue
                short_H, short_coeff = by_pattern[pat]
                if d[0] < short_H:
                    continue
                cand.append((d, short_H, short_coeff))
            for d, short_H, short_coeff in cand:
                excess_H = d[0] - short_H
                if field is not None and field > 1:
                    coeff_adj = (d[-1] * short_coeff) % field
                else:
                    coeff_adj = d[-1] * short_coeff
                alg_iso = Cob.canonical_h_mor(
                    end_clt, self.gens[index], excess_H, coeff_adj)
                if self._cob_isotopy(end, index, alg_iso, col_idx, row_idx):
                    absorbed += 1
        return absorbed

    def _cob_cross_idem_probe(self, start, end, col_idx, row_idx, max_h=2, max_dots=1):
        """Step B of OPEN_QUESTIONS item 5: cross-idempotent absorption.

        For each ``index ∈ col_idx[start]`` (with ``index ≠ end``, ``≠ start``,
        and ``diff[index, end] == 0`` for isotopy safety), enumerate a
        bounded family of canonical single-decoration cobordisms from
        ``gens[end]`` to ``gens[index]`` and try each as an isotopy
        label.  The entropy guard in :meth:`_cob_isotopy` rejects any
        label whose application would increase total decoration count,
        so non-productive candidates cost only the composition evaluation.

        Enumeration bounds:
          - H-powers ``0 .. max_h``
          - Dot subsets of size ``0 .. max_dots`` (over components of
            the specific ``gens[end] → gens[index]`` cobordism).
          - Coefficient ``1`` (we only try the unit case; signed
            coefficients fall out of the arithmetic).

        Returns the number of committed absorptions.
        """
        if self.diff[end, start] == 0:
            return 0
        field = Cob.get_field()
        committed = 0
        targets_of_end = set(col_idx.get(end, ()))
        end_clt = self.gens[end]
        for index in list(col_idx.get(start, ())):
            if index == end or index == start:
                continue
            if index in targets_of_end:
                continue
            if self.diff[index, start] == 0:
                continue
            # Enumerate candidate alg_iso: end → index.
            index_clt = self.gens[index]
            comps = Cob.components(end_clt, index_clt)
            n_comp = len(comps)
            from itertools import combinations
            for h in range(max_h + 1):
                for dot_size in range(max_dots + 1):
                    for dot_subset in combinations(range(n_comp), dot_size):
                        alg_iso = Cob.canonical_dotted_mor(
                            end_clt, index_clt, h, dot_subset, 1)
                        if self._cob_isotopy(end, index, alg_iso, col_idx, row_idx):
                            committed += 1
                            # After a successful isotopy, diff[index, start]
                            # changes — break out to revisit candidates.
                            if self.diff[index, start] == 0:
                                break
                    else:
                        continue
                    break
                else:
                    continue
                break
        return committed

    def clean_up_once(self, verbose=False, cross_idem=False):
        """One pass of Cob-level arrow-shortening (OPEN_QUESTIONS item 5).

        Step A (always on): same-face absorption via
        :meth:`_cob_isolate_h_arrow`.  Polynomial-in-H reduction within
        a fixed dot pattern.

        Step B (``cross_idem=True``): additionally, try a bounded set of
        canonical cross-idempotent cobordisms via
        :meth:`_cob_cross_idem_probe`.  This is the experimental
        generalization for item 5's cross-idempotent case; expensive,
        and the entropy guard filters non-productive candidates.

        Returns total absorptions committed.
        """
        candidates = [(s, t) for (t, s) in self._nnz if s != t]
        if not candidates:
            return 0
        col_idx = self._nnz_by_col()
        row_idx = self._nnz_by_row()
        count = 0
        for (start, end) in candidates:
            count += self._cob_isolate_h_arrow(start, end, col_idx, row_idx)
        if cross_idem:
            col_idx = self._nnz_by_col()
            row_idx = self._nnz_by_row()
            candidates_b = [(s, t) for (t, s) in self._nnz if s != t]
            for (start, end) in candidates_b:
                count += self._cob_cross_idem_probe(start, end, col_idx, row_idx)
        if verbose:
            print("  clean_up_once: {} absorptions on {} candidates (cross_idem={})".format(
                count, len(candidates), cross_idem))
        return count
        return count

    def clean_up(self, max_iter=50):
        """Iterate :meth:`clean_up_once` and :meth:`eliminateAll`
        alternately until no further reduction occurs.

        Breaks the super-linear growth from wide-width tangles (e.g.,
        TH-pattern at (1, 13)) if the main source of that growth is
        arrows whose labels differ only by a power of H — which was
        invisible to :meth:`eliminateAll` before.  See OPEN_QUESTIONS.md
        item 5 for the scope of what this covers vs. the open research
        questions (dotted / cross-idempotent shortening).
        """
        last = None
        for _ in range(max_iter):
            self.eliminateAll()
            self.clean_up_once()
            self.eliminateAll()
            current = len(self.gens)
            if current == last:
                break
            last = current
        return self

    def shift_qhd(self,q,h,delta):
        self.gens=[clt.shift_qhd(q,h,delta) for clt in self.gens]

    def ToBNAlgebra(self,field=2):
        gens=[clt.ToBNAlgebra() for clt in self.gens]
        def convert(cob):
            if cob == 0:
                return 0
            return cob.ToBNAlgebra(field)
        diff=[[convert(cob) for cob in row] for row in self.diff]
        BNcx = BNComplexes.BNComplex(gens,diff,field)
        BNcx.eliminateAll()
        # BNcx.validate()
        return BNcx

def AddCapToCLT(clt, i, grshift = "false"): 
    """creates a new CLT which is clt with a cap added to it at index i
       Here 0 <= i <= tangle.bot
       If grshift is set to "true", then it also shifts the homological and quantum grading of the clt by 1 each
       This should only be used when adding a crossing"""
    def incrementby2(j): #increment TEI by 2 if greater than where cap is to be inserted
        if j >= clt.top+i:
            return j+2
        else:
            return j
    newarcs = [incrementby2(x) for x in clt.arcs[:clt.top+i]] + [clt.top+i +1, clt.top+i] + [incrementby2(x) for x in clt.arcs[clt.top+i:]]
    if grshift == "true": # shift grading
        newh = clt.h +1
        newq = clt.q +1
        newdelta = clt.delta -0.5
    else: # dont shift grading
        newh = clt.h
        newq = clt.q
        newdelta = clt.delta
    return Cob.obj(clt.top, clt.bot+2, newarcs, newh, newq, newdelta)

def AddCap(Complex, i, grshift = "false"):
    """ Creates a new complex by adding a cap to every tangle and every cobordism in Complex, at index i
        Here 0 <= i <= tangle.bot
        If grshift is set to "true", then it applies a grading shift to every tangle (as above)
        Furthermore, it will flip the sign on every cobordism, which is the convention used so that the differential will square to 0 when adding a crossing"""
    Newgens = [AddCapToCLT(clt, i, grshift) for clt in Complex.gens]
    N = len(Complex.gens)
    Newdiff = np.zeros((N, N), dtype=object)
    nnz_hint = []
    for target, source in Complex._nnz:
        cob = Complex.diff[target, source]
        top_plus_i = cob.front.top + i
        def incrementby2(j, _base=top_plus_i):
            return j + 2 if j >= _base else j
        newcomps = [[incrementby2(x) for x in comp] for comp in cob.comps] + [[top_plus_i, top_plus_i + 1]]
        if grshift == "true":
            newDecos = [NewDeco[:-1] + [0] + [NewDeco[-1]*-1] for NewDeco in cob.decos]
        else:
            newDecos = [NewDeco[:-1] + [0] + NewDeco[-1:] for NewDeco in cob.decos]
        NewCob = Cob.mor(Newgens[source], Newgens[target], Cob.simplify_decos(newDecos), newcomps)
        reduced = NewCob.ReduceDecorations()
        if reduced != 0:
            Newdiff[target, source] = reduced
            nnz_hint.append((target, source))
    return CobComplex(Newgens, Newdiff, nnz_hint=nnz_hint)

def AddCupToCLT(clt, i):
    """ Adds a cup to the clt at index i, where 0 <= i <= clt.bot -2
        Returns a list of 2 clt if there is a closed component, obtained by neckcutting
        Otherwise returns a list of 1 clt """
    newgens = []
    def decrementby2(j):#decrement TEI by 2 if greater than where cup is to be inserted
            if j >= clt.top +i:
                return j-2
            else:
                return j
    if clt.arcs[clt.top +i] == clt.top+i+1: # adding the cup makes a closed component
        newarcs = [decrementby2(x) for x in clt.arcs if x != clt.top+i and x!= clt.top+i+1] #removes the closed component from the arcs, shifts remaining TEI
        newgens.append(Cob.obj(clt.top, clt.bot -2, newarcs, clt.h, clt.q+1, clt.delta+0.5)) #neckcutting
        newgens.append(Cob.obj(clt.top, clt.bot -2, newarcs, clt.h, clt.q-1, clt.delta-0.5)) #neckcutting
    else: # adding the cup doesnt make closed components
        leftend = clt.arcs[clt.top + i] #the endpoint of the arc which connects at i
        rightend = clt.arcs[clt.top+i+1] #the endpoint of the arc which connects at i+1
        newarcs = clt.arcs.copy()
        newarcs[leftend] = rightend # connecting the two arcs via the cup
        newarcs[rightend] = leftend
        newarcs1 = [decrementby2(entry) for x, entry in enumerate(newarcs) if x != clt.top+i and x != clt.top+i+1] #removes the ends which don't exist anymore and shifts remaining TEI
        newgens.append(Cob.obj(clt.top, clt.bot -2, newarcs1, clt.h, clt.q, clt.delta))
    return newgens

def AddCup(Complex, i): # TODO: reduce decorations
    """ Here 0 <= i <= tangle.bot -2"""
    # Precompute per-gen: does this CLT become "closed" under the cup, and
    # what's the starting new-index?  Closed gens contribute 2 new gens, open
    # gens contribute 1.  This lets us pre-size the output matrix and place
    # each cob at its correct (row, col) without the N*N Python iteration.
    old_gens = Complex.gens
    N = len(old_gens)
    closed = [g.arcs[g.top + i] == g.top + i + 1 for g in old_gens]
    new_start = [0] * N
    acc = 0
    for k in range(N):
        new_start[k] = acc
        acc += 2 if closed[k] else 1
    new_N = acc

    newgens = []
    for clt in old_gens:
        newgens.extend(AddCupToCLT(clt, i))

    # Cache AddCupToCLT results (used per old-gen up to twice in the original,
    # once for source and once for target).
    new_clts_per_gen = [None] * N

    def get_new_clts(k):
        if new_clts_per_gen[k] is None:
            new_clts_per_gen[k] = AddCupToCLT(old_gens[k], i)
        return new_clts_per_gen[k]

    Newdiff = np.zeros((new_N, new_N), dtype=object)
    nnz_hint = []

    # We only need to visit non-zero entries of the old diff.  The inner body
    # below is the original 4-branch logic, modified to *write* into the
    # pre-allocated matrix at the correct offsets rather than appending to
    # row lists.
    for target, source in Complex._nnz:
        cob = Complex.diff[target, source]
        base = cob.front.top + i

        def decrementby2(j, _b=base):
            return j - 2 if j >= _b else j

        def compcup(component):
            return [decrementby2(j) for j in component if j != base and j != base + 1]

        def incrementH(decoration):
            return [decoration[0] + 1] + decoration[1:]

        def negativeH(decoration):
            return [decoration[0] + 1] + decoration[1:-1] + [decoration[-1] * -1]

        def adddotonindex(decoration, idx):
            return decoration[:idx + 1] + [1] + decoration[idx + 2:]

        src_closed = closed[source]
        tgt_closed = closed[target]
        row0 = new_start[target]
        col0 = new_start[source]

        if src_closed and tgt_closed:
            # both closed: produce Cob1 at (r0,c0), Cob3 at (r0+1,c0), Cob4 at (r0+1,c0+1).
            # Cob2 (at (r0,c0+1)) is 0 and is already zero in Newdiff.
            magic_index = -1
            for x, comp in enumerate(cob.comps):
                if base in comp:
                    if len(comp) != 2:
                        raise Exception("closed component doesn't have 2 TEI")
                    magic_index = x
                    break
            if magic_index == -1:
                raise Exception("no magic_index")
            newcomps = [[decrementby2(j) for j in comp] for comp in cob.comps if base not in comp]

            def delclosedcomp(deco, _mi=magic_index):
                return deco[:_mi + 1] + deco[_mi + 2:]

            def computedeco4(deco, _mi=magic_index):
                if deco[_mi + 1] == 0:
                    return delclosedcomp(deco)
                else:
                    return delclosedcomp(incrementH(deco))

            newDecos1 = [delclosedcomp(deco) for deco in cob.decos if deco[magic_index + 1] == 0]
            newDecos3 = [delclosedcomp(deco) for deco in cob.decos if deco[magic_index + 1] == 1]
            newDecos4 = [computedeco4(deco) for deco in cob.decos]
            nsc = get_new_clts(source)
            ntc = get_new_clts(target)
            C1 = Cob.mor(nsc[0], ntc[0], Cob.simplify_decos(newDecos1), newcomps).ReduceDecorations()
            C3 = Cob.mor(nsc[0], ntc[1], Cob.simplify_decos(newDecos3), newcomps).ReduceDecorations()
            C4 = Cob.mor(nsc[1], ntc[1], Cob.simplify_decos(newDecos4), newcomps).ReduceDecorations()
            if C1 != 0: Newdiff[row0,     col0]     = C1; nnz_hint.append((row0, col0))
            if C3 != 0: Newdiff[row0 + 1, col0]     = C3; nnz_hint.append((row0 + 1, col0))
            if C4 != 0: Newdiff[row0 + 1, col0 + 1] = C4; nnz_hint.append((row0 + 1, col0 + 1))
        elif tgt_closed:  # target closed, source open
            magic_index = -1
            for x, comp in enumerate(cob.comps):
                if base in comp:
                    magic_index = x
                    break
            if magic_index == -1:
                raise Exception("no magic_index")
            newcomps = [compcup(component) for component in cob.comps]
            newDecos1 = []
            for deco in cob.decos:
                if deco[magic_index + 1] == 0:
                    newDecos1.extend([adddotonindex(deco, magic_index), negativeH(deco)])
            newDecos2 = [deco for deco in cob.decos]
            nsrc = get_new_clts(source)[0]
            ntc = get_new_clts(target)
            C1 = Cob.mor(nsrc, ntc[0], Cob.simplify_decos(newDecos1), newcomps).ReduceDecorations()
            C2 = Cob.mor(nsrc, ntc[1], Cob.simplify_decos(newDecos2), newcomps).ReduceDecorations()
            if C1 != 0: Newdiff[row0,     col0] = C1; nnz_hint.append((row0, col0))
            if C2 != 0: Newdiff[row0 + 1, col0] = C2; nnz_hint.append((row0 + 1, col0))
        elif src_closed:  # source closed, target open
            magic_index = -1
            for x, comp in enumerate(cob.comps):
                if base in comp:
                    magic_index = x
                    break
            if magic_index == -1:
                raise Exception("no magic_index")
            newcomps = [compcup(component) for component in cob.comps]

            def computedeco2(deco, _mi=magic_index):
                if deco[_mi + 1] == 1:
                    return incrementH(deco)
                else:
                    return adddotonindex(deco, _mi)

            newDecos1 = [deco for deco in cob.decos]
            newDecos2 = [computedeco2(deco) for deco in cob.decos]
            nsc = get_new_clts(source)
            ntgt = get_new_clts(target)[0]
            C1 = Cob.mor(nsc[0], ntgt, Cob.simplify_decos(newDecos1), newcomps).ReduceDecorations()
            C2 = Cob.mor(nsc[1], ntgt, Cob.simplify_decos(newDecos2), newcomps).ReduceDecorations()
            if C1 != 0: Newdiff[row0, col0]     = C1; nnz_hint.append((row0, col0))
            if C2 != 0: Newdiff[row0, col0 + 1] = C2; nnz_hint.append((row0, col0 + 1))
        else:  # both open
            magic_index = -1
            for x1, comp in enumerate(cob.comps):
                if base in comp:
                    magic_index = x1
                    break
            if magic_index == -1:
                raise Exception("no magic_index")
            magic_index_2 = -1
            for x2, comp in enumerate(cob.comps):
                if base + 1 in comp:
                    magic_index_2 = x2
                    break
            if magic_index_2 == -1:
                raise Exception("no magic_index_2")

            newDecos1 = []
            if magic_index == magic_index_2:
                comp = cob.comps[magic_index]
                x1 = comp.index(base)
                x2 = comp.index(base + 1)
                comp1 = comp[:min(x1, x2)] + comp[max(x1, x2) + 1:]
                comp2 = comp[min(x1, x2) + 1:max(x1, x2)]
                newcomps = [[decrementby2(j) for j in comp] for comp in cob.comps[:magic_index] + [comp1, comp2] + cob.comps[magic_index + 1:]]

                def computedeco1comps(deco, _mi=magic_index):
                    if deco[_mi + 1] == 1:
                        return [deco[:_mi + 1] + [1] + deco[_mi + 1:]]
                    else:
                        return [deco[:_mi + 1] + [1] + deco[_mi + 1:],
                                deco[:_mi + 2] + [1] + deco[_mi + 2:],
                                [deco[0] + 1] + deco[1:_mi + 1] + [0] + deco[_mi + 1:-1] + [deco[-1] * -1]]

                for deco in cob.decos:
                    newDecos1.extend(computedeco1comps(deco))
            else:
                comp1 = cob.comps[magic_index]
                comp2 = cob.comps[magic_index_2]
                location1 = comp1.index(base)
                location2 = comp2.index(base + 1)
                comp2 = (comp2[location2 + 1:] + comp2[:location2])
                if location1 % 2 == location2 % 2:
                    comp2.reverse()
                comp1 = comp1[:location1] + comp2 + comp1[location1 + 1:]
                newcomps = [[decrementby2(j) for j in comp] for comp in cob.comps[:magic_index] + [comp1] + cob.comps[magic_index + 1:]]
                del newcomps[magic_index_2]

                def computedeco2comps(deco, _mi=magic_index, _mi2=magic_index_2):
                    if deco[_mi + 1] == 1 and deco[_mi2 + 1] == 1:
                        return incrementH(deco)[:_mi2 + 1] + incrementH(deco)[_mi2 + 2:]
                    elif deco[_mi + 1] == 1 or deco[_mi2 + 1] == 1:
                        return adddotonindex(deco, _mi)[:_mi2 + 1] + adddotonindex(deco, _mi)[_mi2 + 2:]
                    else:
                        return deco[:_mi2 + 1] + deco[_mi2 + 2:]

                newDecos1 = [computedeco2comps(deco) for deco in cob.decos]

            nsrc = get_new_clts(source)[0]
            ntgt = get_new_clts(target)[0]
            C1 = Cob.mor(nsrc, ntgt, Cob.simplify_decos(newDecos1), newcomps).ReduceDecorations()
            if C1 != 0: Newdiff[row0, col0] = C1; nnz_hint.append((row0, col0))

    return CobComplex(newgens, Newdiff, nnz_hint=nnz_hint)
  
def AddPosCrossing(Complex, i):
    CapCup = AddCap(AddCup(Complex, i), i, "true")
    sourcegens = Complex.gens
    targetgens = CapCup.gens
    N = len(sourcegens)
    M = len(targetgens)

    TopLeft = Complex.diff
    BottomRight = CapCup.diff
    TopRight = np.zeros((N, M), dtype=object)

    # BottomLeft is an M x N matrix whose only non-zero cells lie on the
    # "x == y" diagonal (1 or 2 cells per x depending on whether targetclt
    # is closed).  Iterate x directly instead of the full M*N grid.
    BottomLeft = np.zeros((M, N), dtype=object)
    bl_nnz = []  # (row, col) indices written into BottomLeft
    row_cursor = 0
    for x, targetclt in enumerate(sourcegens):
        sourceclt = sourcegens[x]  # x == y case of the original
        if targetclt.arcs[targetclt.top + i] == targetclt.top + i + 1:
            newTarget1 = AddCapToCLT(AddCupToCLT(targetclt, i)[0], i, "true")
            newTarget2 = AddCapToCLT(AddCupToCLT(targetclt, i)[1], i, "true")
            newcomps = Cob.components(sourceclt, newTarget1)
            magic_index = 0
            for z, comp in enumerate(newcomps):
                if sourceclt.top + i in comp:
                    magic_index = z
                    break
            decos1 = [[0] + [0 for _ in newcomps[:magic_index]] + [1] + [0 for _ in newcomps[magic_index + 1:]] + [1],
                      [1] + [0 for _ in newcomps] + [-1]]
            decos2 = [[0] + [0 for _ in newcomps] + [1]]
            BottomLeft[row_cursor,     x] = Cob.mor(sourceclt, newTarget1, decos1, newcomps)
            BottomLeft[row_cursor + 1, x] = Cob.mor(sourceclt, newTarget2, decos2, newcomps)
            bl_nnz.append((row_cursor, x)); bl_nnz.append((row_cursor + 1, x))
            row_cursor += 2
        else:
            newTarget = AddCapToCLT(AddCupToCLT(targetclt, i)[0], i, "true")
            decos = [[0] + [0 for _ in Cob.components(sourceclt, newTarget)] + [1]]
            BottomLeft[row_cursor, x] = Cob.mor(sourceclt, newTarget, decos)
            bl_nnz.append((row_cursor, x))
            row_cursor += 1

    Newdiff = np.concatenate((np.concatenate((TopLeft, TopRight), axis=1),
                              np.concatenate((BottomLeft, BottomRight), axis=1)), axis=0)
    # Assemble nnz_hint from the four blocks:
    #   TopLeft block (N x N) at offsets (0, 0)       — nnz from Complex._nnz
    #   TopRight (N x M)                              — all zero
    #   BottomLeft (M x N) at offsets (N, 0)          — bl_nnz
    #   BottomRight (M x M) at offsets (N, N)         — from CapCup._nnz
    nnz_hint = list(Complex._nnz)
    nnz_hint.extend((N + r, c) for (r, c) in bl_nnz)
    nnz_hint.extend((N + t, N + s) for (t, s) in CapCup._nnz)
    return CobComplex(sourcegens + targetgens, Newdiff, nnz_hint=nnz_hint)

def grshiftclt(clt):
    return Cob.obj(clt.top, clt.bot, clt.arcs, clt.h+1, clt.q+1, clt.delta-0.5)

def grshiftcob(cob):
    if cob == 0:
        return 0
    newDecos = [deco[:-1] + [deco[-1]*-1] for deco in cob.decos]
    return Cob.mor(grshiftclt(cob.front), grshiftclt(cob.back), newDecos, cob.comps)
    
def AddNegCrossing(Complex, i):
    targetgens = [grshiftclt(clt) for clt in Complex.gens]
    CapCup = AddCap(AddCup(Complex, i), i)
    sourcegens = CapCup.gens
    N = len(Complex.gens)
    M = len(sourcegens)

    TopLeft = CapCup.diff

    # BottomRight = grshift of the original diff.  Build sparsely: iterate
    # only non-zero cells of Complex.diff rather than N*N Python cells.
    BottomRight = np.zeros((N, N), dtype=object)
    for t, s in Complex._nnz:
        BottomRight[t, s] = grshiftcob(Complex.diff[t, s])

    TopRight = np.zeros((M, N), dtype=object)

    # BottomLeft has the same diagonal-only structure as in AddPosCrossing:
    # non-zero only where the column index (source) matches the row's
    # originating generator.  Iterate columns directly.
    BottomLeft = np.zeros((N, M), dtype=object)
    bl_nnz = []
    col_cursor = 0
    for x, sourceclt in enumerate(Complex.gens):
        newTarget = targetgens[x]
        if sourceclt.arcs[sourceclt.top + i] == sourceclt.top + i + 1:
            newSource1 = AddCapToCLT(AddCupToCLT(sourceclt, i)[0], i)
            newSource2 = AddCapToCLT(AddCupToCLT(sourceclt, i)[1], i)
            newcomps = Cob.components(newSource2, newTarget)
            magic_index = 0
            for z, comp in enumerate(newcomps):
                if sourceclt.top + i in comp:
                    magic_index = z
                    break
            decos1 = [[0] + [0 for _ in newcomps] + [1]]
            decos2 = [[0] + [0 for _ in newcomps[:magic_index]] + [1] + [0 for _ in newcomps[magic_index + 1:]] + [1]]
            BottomLeft[x, col_cursor]     = Cob.mor(newSource1, newTarget, decos1, newcomps)
            BottomLeft[x, col_cursor + 1] = Cob.mor(newSource2, newTarget, decos2, newcomps)
            bl_nnz.append((x, col_cursor)); bl_nnz.append((x, col_cursor + 1))
            col_cursor += 2
        else:
            newSource = AddCapToCLT(AddCupToCLT(sourceclt, i)[0], i)
            decos = [[0] + [0 for _ in Cob.components(newSource, newTarget)] + [1]]
            BottomLeft[x, col_cursor] = Cob.mor(newSource, newTarget, decos)
            bl_nnz.append((x, col_cursor))
            col_cursor += 1

    Newdiff = np.concatenate((np.concatenate((TopLeft, TopRight), axis=1),
                              np.concatenate((BottomLeft, BottomRight), axis=1)), axis=0)
    # Assemble nnz_hint:
    #   TopLeft (M x M) at offset (0, 0)            — CapCup._nnz
    #   TopRight (M x N)                            — all zero
    #   BottomLeft (N x M) at offset (M, 0)         — bl_nnz (col_index offset 0 since BL is (N, M))
    #     ...but BottomLeft is placed to the LEFT of BottomRight, so in the
    #     concatenated matrix its columns start at 0. Rows offset M.
    #   BottomRight (N x N) at offset (M, M)        — from Complex._nnz shifted
    nnz_hint = list(CapCup._nnz)
    nnz_hint.extend((M + r, c) for (r, c) in bl_nnz)
    nnz_hint.extend((M + t, M + s) for (t, s) in Complex._nnz)
    return CobComplex(sourcegens + targetgens, Newdiff, nnz_hint=nnz_hint)

def BNbracket(string,pos=0,neg=0,start=1,options="unsafe",cleanup_field=None,signed_lift=False):
    """compute the Bar-Natan bracket for tangle specified by 'string', which is a concatenation of words <type>+<index>, separated by '.' read from right to left, for each elementary tangle slice, read from top to bottom, where:
    <type> is equal to:
        'pos': positive crossing
        'neg': negative crossing
        'cup': cup
        'cap': cap
    <index> is the index at which the crossing, cap or cup sits.
    'pos' and 'neg' are the numbers of positive and negative crossings.
    The first optional parameter 'start' is an integer which specifies the number of tangle ends at the top.
    The second optional parameter 'options' is either 'unsafe' (default) or 'safe'. The latter performs some sanity checks, but slows down the computation.
    The third optional parameter 'cleanup_field' enables intermediate (1,3)
    cleanup: whenever the running complex is in a (1,3) state, convert to
    BNAlgebra(cleanup_field), run the immersed-curve clean_up, and convert
    back via BNComplex.ToCob.  The round-trip is mathematically correct
    (final invariants match), but is EXPERIMENTAL and usually hurts
    performance: the Cob category has no field knowledge, so coefficients
    returned from ToCob are integer representatives in [0, p) and will
    inflate through subsequent slice operations until the final
    ToBNAlgebra(field) mods them out.  On a 3-braid test this made the
    total ~50x slower.  Only worth turning on for tangles with many (1,3)
    checkpoints AND a follow-up cob-level mod-p simplification (not yet
    implemented).  Pass None to disable (default).
    The fourth optional parameter 'signed_lift' controls the coefficient
    lift when converting back from BNComplex to Cob inside the
    intermediate cleanup: False (default) uses [0, p); True centers on
    (-p/2, p/2].  See OPEN_QUESTIONS.md item 1.
    E.g. 'BNbracket('cup0pos0',2)' is a (2,0)-tangle which is decomposed as a positive crossing followed by a cap.
    """
    stringlist=[[word[0:3],int(word[3:])] for word in string.split('.')]
    # Activate F_p arithmetic in the Cob layer for the duration of this
    # BNbracket call when intermediate cleanup is requested.  This keeps the
    # integer representatives returned by ToCob small (and simplify_decos
    # automatically collapses 1+1=2 -> 0 in F_2, etc.).  Reverted in the
    # finally-block at the end of this function.
    # F_p arithmetic in Cob is activated lazily the first time an actual
    # intermediate cleanup fires.  This avoids paying mod-p overhead on
    # tangles where cleanup_field was requested but no (1,3) checkpoint
    # ever produced simplification (e.g. 3braid_scaled, whose only
    # candidate is a 1-gen slice).
    _fp_state = {"prev": None, "active": False}

    def _activate_fp():
        if not _fp_state["active"]:
            _fp_state["prev"] = Cob.set_field(cleanup_field)
            _fp_state["active"] = True

    try:
        cx=CobComplex([Cob.obj(start,start,[start+i for i in range(start)]+[i for i in range(start)], 0,0,0)], [[0]])
        print("Computing the Bar-Natan bracket for the tangle\n\n"+string+"\n\n"+"with "+str(start)+" ends at the top, "+str(pos)+\
              " positive crossings, "+str(neg)+" negative crossings and "+str(len(stringlist))+" slices in total.")

        time0=time()
        time1=time0

        def _maybe_intermediate_cleanup(cx, is_last):
            # Safe only when all gens are (1,3)-CLTs (i.e., the running tangle has
            # shrunk to a 4-ended state).  The BN-algebra round-trip then gives
            # an aggressive simplification, because BNComplex.clean_up applies
            # immersed-curve arrow-shortening isotopies that Gaussian elimination
            # alone can't see.
            #
            # Skip conditions:
            #   - cleanup_field is None: feature not requested.
            #   - is_last: caller will do ToBNAlgebra + eliminateAll anyway,
            #     so intermediate cleanup here would be duplicate work.
            #   - len(cx.gens) <= 1: nothing for clean_up to do.
            if cleanup_field is None or not cx.gens or is_last or len(cx.gens) <= 1:
                return cx
            g0 = cx.gens[0]
            if not (g0.top == 1 and g0.bot == 3):
                return cx
            before = len(cx.gens)
            BN = cx.ToBNAlgebra(cleanup_field)
            BN.clean_up()
            cx2 = BN.ToCob(signed_lift=signed_lift)
            after = len(cx2.gens)
            # Now that ToCob has emitted integer representatives in [0, p),
            # flip Cob into F_p mode for the remaining slices so coefficients
            # stay small through subsequent AddCup/AddPos/etc.
            _activate_fp()
            if after < before:
                print("intermediate-cleanup: (1,3) reduction "+str(before)+" -> "+str(after)+" gens", end='\r')
            return cx2

        num_slices = len(stringlist)
        for i,word in enumerate(stringlist):

            time2=time()
            print("slice "+str(i)+"/"+str(len(stringlist))+": adding "+word[0]+" at index "+str(word[1])+" to tangle. ("+str(len(cx.gens))+" objects, "+str(round(time2-time1,1))+" sec)", end='\r')# monitor \n ->\r
            time1=time2

            if word[0]=="pos":
                cx=AddPosCrossing(cx, word[1])
                if options=="safe": cx.validate()
                cx.eliminateAll()

            elif word[0]=="neg":
                cx=AddNegCrossing(cx, word[1])
                if options=="safe": cx.validate()
                cx.eliminateAll()

            elif word[0]=="cup":
                cx=AddCup(cx, word[1])
                if options=="safe": cx.validate()
                cx.eliminateAll()

            elif word[0]=="cap":
                cx=AddCap(cx, word[1])
                if options=="safe": cx.validate()

            else:
                print("this should never execute")

            cx = _maybe_intermediate_cleanup(cx, is_last=(i == num_slices - 1))

        cx.shift_qhd(pos-2*neg,-neg,0.5*neg)

        print("Completed the computation successfully after "+str(round(time1-time0,1))+" second(s).                        ")
        return cx
    finally:
        if _fp_state["active"]:
            Cob.set_field(_fp_state["prev"])

def importCobcx(filename):
    with open("examples/data/CobComplexes/"+filename, "r") as text_file:
        data = text_file.read()
        return eval(data)

