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

from fractions import Fraction
from itertools import groupby

import Cob

# Turn an integer into superscript, and returns it as a string
def ToExponent(exponent):
    return str(exponent).translate(str.maketrans("-0123456789.", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹·"))

def inverse(num,field): #this only works over a field
    """Finds the multiplicative inverse of num over a field
       Note that num cannot be 0
       if field is 0, then the inversion is over Q
       if field is 1, then the 'field' is Z, and the only inverible elements are 1, -1
       if field is a prime p, then the inversion is over Z_p, using either fermats little theorem or Euclidean algorithm"""
    if field == 0:
        return Fraction(1)/num
    elif field == 1:
        if num in [1, -1]:
            return num
        else: 
            raise Exception("Can't invert over Z")
    else: 
        return pow(num, field -2, field) # num^{field-2} mod (field) #TODO: check that this actually works as intended
    #taken from https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
        # s = 0
        # S = 1
        # r = field
        # R = num
        # while r != 0:
            # q = R // r
            # R, r = r, R - q * r
            # S, s = s, S - q * s
        # #print((S*num)%self.field) # should be 1 if computed correctly
        # return S

class obj(object):
    """An object is a pair [idempotent,q,h,delta(optional)], where idempotent is either 0 (b=solid dot) or 1 (c=hollow dot). 
    """
    __slots__ = 'idem','q','h','delta'
    
    def __init__(self,idempotent,q,h,delta="default"):
        self.idem = idempotent
        self.q = q
        self.h = h
        if delta=="default":
            self.delta = q/2-h
        else:
            self.delta = delta
    
    def __repr__(self):
        return "BNAlgebra.obj({},{},{},{})".format(self.idem,self.q,self.h,self.delta)
    
    def idem2dot(self):
        if self.idem==0:
            return "●"#b (solid dot)
        else: 
            return "○"#c (hollow dot)
    
    def grading2TeX(self,switch="qhdelta"):
        if "q" in switch:
            q="q^{"+str(self.q)+"}"
        else:
            q=""
        if "h" in switch:
            h="h^{"+str(self.h)+"}"
        else:
            h=""
        if "delta" in switch:
            if 2*self.delta % 2 == 0: # delta is an integer
                delta="\\delta^{"+str(round(self.delta))+"}"
            else:
                #delta="δ"+ToExponent(round(2*self.delta))+"'²"
                delta="\\delta^{"+str(self.delta)+"}"
            
        else:
            delta=""
        
        return q+h+delta
    
    def obj2string(self,switch="idem",index=-1): 
        """ Returns a string version of the morphism including grading information if specified
            switch is a string containing idem, index, q, h, delta
            index specifies the index of the object when in a complex
            idem is either b or c
            q, h, delta are the grading informations presented as an exponent eg h^3
        """
            
        if "idem" in switch:
            idem=self.idem2dot()
        else:
            idem = ""
        
        if (index == -1) or ("index" not in switch):
            index = ""
        else:
            index = str(index)
        if "q" in switch:
            q="q"+ToExponent(self.q)
        else:
            q=""
        if "h" in switch:
            h="h"+ToExponent(self.h)
        else:
            h=""
        if "delta" in switch:
            if 2*self.delta % 2 == 0: # delta is an integer
                delta="δ"+ToExponent(round(self.delta))
            else:
                #delta="δ"+ToExponent(round(2*self.delta))+"'²"
                delta="δ"+ToExponent(self.delta)
            
        else:
            delta=""
        
        grading=q+h+delta
        
        if (grading == "") or (index == ""):
            return index+grading+idem
        else:
            return index+":"+grading+idem
    
    def shift_q(self,shift): #shift q, keep h fixed; create new object
        return obj(self.idem,self.q+shift,self.h,self.delta+shift/2)
        
    def shift_h(self,shift): #shift h, keep q fixed; create new object
        return obj(self.idem,self.q,self.h+shift,self.delta-shift)
    
    def ToCob(self): #turns a BN obj into a Cob obj (A CLT)
        arcs = []
        if self.idem == 0:
            arcs = [3,2,1,0]
        else:
            arcs = [1,0,3,2]
        return Cob.obj(1,3, arcs, self.h, self.q, self.delta)

def coeff_simplify(num,field): # Returns num mod p if field is a prime p, otherwise do nothing
    if field > 1:
        return num % field
    else: # This probably needs to be fixed for field=0
        return num

# OPEN_QUESTIONS.md item 3: cabled knots 5_1, 5_2, 6_1 fail d^2 = 0 over
# F_3, F_5, F_7.  The most likely culprit is a sign mismatch between the
# even-negative and odd-negative branches of mor.ToCob.  Setting
# _sign_convention = "corrected" flips the sign on the D-term in the
# even-negative branch.  Keep "legacy" as the default so historical
# outputs are unchanged; "corrected" is experimental and not yet
# mathematically verified.
_sign_convention = "legacy"

def set_sign_convention(convention):
    """Flip between legacy and experimental sign conventions in
    mor.ToCob (see OPEN_QUESTIONS.md item 3).  Returns the previous
    value so callers can restore it in a try/finally."""
    global _sign_convention
    if convention not in ("legacy", "corrected"):
        raise ValueError("sign_convention must be 'legacy' or 'corrected'")
    prev = _sign_convention
    _sign_convention = convention
    return prev

def signed_lift_coeff(num, field):
    """Centered representative of ``num`` modulo ``field``.

    Maps [0, p) to (-p/2, p/2]; pass-through when field <= 1.
    Used when converting BN-algebra coefficients back to Cob, where
    downstream Cob identities might prefer a signed lift.  See
    OPEN_QUESTIONS.md item 1 — the math question of which lift makes
    which identities hold is open, so this is opt-in via a flag.
    """
    if field is None or field <= 1:
        return num
    reduced = num % field
    half = field // 2
    if reduced > half:
        return reduced - field
    return reduced

class mor(object):
    """An element of Bar-Natan's algebra is a list of pairs [power,coeff]
    'power' is an integer, which determines the exponent of D (if positive) and the exponent of S (if negative)
    'coeff' is some non-zero integer (= coefficient in the base ring/field) # Alternatively, a Fraction object
    field == 0 means Q
    field == 1 means Z
    field == p means F_p for prime p
    """
    __slots__ = 'pairs','field'
    
    def __init__(self,pairs,field):
        self.pairs = pairs
        self.field = field
    
    def __repr__(self):
        return "BNAlgebra.mor({},{})".format(self.pairs,self.field)
    
    def simplify_mor(self,field):
        """simplify algebra elements by adding all coeffients of the same power of D or S, omitting those with coefficient 0. This is very similar to simplify_decos"""
        def droplast(l):
            return l[:-1]
        def add_coeffs(iterable):
            coeff=0
            for x in iterable:
                coeff+=x[-1]
            return coeff_simplify(coeff,field)
        self.pairs = [power+[add_coeffs(grouped)] for power,grouped in groupby(sorted(self.pairs),droplast)]
        self.pairs = [x for x in self.pairs if x[-1]!=0]
        if self.pairs == []:
            return 0
        return self
    
    def __add__(self, other):
        if other == 0:
            return self
        return mor(self.pairs+other.pairs,self.field).simplify_mor(self.field)
    
    def __radd__(self, other):
        if other == 0:
            return self
        return mor(self.pairs+other.pairs,self.field).simplify_mor(self.field)

    def __mul__(self, other):
        if other == 0:
            return 0
        if isinstance(other,mor):
            return mor([[a1[0]+a2[0],a1[1]*a2[1]] for a1 in self.pairs for a2 in other.pairs if a1[0]*a2[0]>=0],self.field).simplify_mor(self.field)
        # 'other' is assumed to be a non-zero integer
        return mor([[pair[0],other*pair[1]] for pair in self.pairs],self.field).simplify_mor(self.field)
        
    def __rmul__(self, other):
        if other == 0:
            return 0
        return mor([[a1[0]+a2[0],a1[1]*a2[1]] for a1 in self.pairs for a2 in other.pairs if a1[0]*a2[0]>=0],self.field).simplify_mor(self.field)
    
    def is_identity(self): #Returns true only if there is a single pair, which has no power of D, and the coefficient is +- 1
        if len(self.pairs)!=1:
            return False
        elif self.pairs[0][0]!=0:
            return False
        elif self.pairs[0][1] in [1,-1]: #TODO: check -1 in F_p
            return True
        else:
            return False
            
    def is_isomorphism(self): #Returns true only if there is a single pair, which has no power of D, and the coefficient is invertible over field
        if len(self.pairs)!=1: # if there is more than 1 pair
            return False
        elif self.pairs[0][0]!=0: # if the pair has a non-zero power of D or S
            return False
        elif self.pairs[0][1]==0: # if the coefficient is 0 it can't be inverted
            return False
        elif self.field == 1 and self.pairs[0][1] not in [-1, 1]: # if the coefficient is not +- 1 and field is Z, it can't be inverted
            return False 
        return True
    
    def contains_D(self): #Returns true if self contains a pair with a power of D
        return all([pair[0]<=0 for pair in self.pairs])==False
    
    def contains_S(self): #Returns true if self contains a pair with a power of S
        return all([pair[0]>=0 for pair in self.pairs])==False
    
    def __neg__(self): #returns a new mor that is the same except with the opposite sign on the coefficient
        return mor([[pair[0],(-1)*pair[1]] for pair in self.pairs],self.field)
    
    def label2TeX(self,switch="DS"):
        """The switch determines whether we only print powers of S or powers of D or both. The identity is included as D^0, not as S^0."""
        string=""
            
        for pair in self.pairs:
            if ((("D" in switch) and (pair[0]>=0)) or (("S" in switch) and (pair[0]<0))):
                coeff = pair[1]
                if (string != "") & (coeff > 0):# add plus sign if the next coefficient is positive, except for the first summand
                    string += "+"
                if coeff < 0: # add minus sign in any case
                    string += "-"
                    coeff = abs(coeff)
                
                if coeff != 0:# omit any summands with coefficient 0
                
                    exponent=abs(pair[0])
                    if exponent==1: # omit exponent 1 from the notation
                        exponent = ""
                    else:
                        exponent= "^{"+str(exponent)+"}"
                    
                    if coeff==1: # omit coefficients 1 and -1 from the notation
                        coeff = ""
                    else:
                        coeff = str(coeff) + "\cdot "
                    
                    if pair[0] > 0:# powers of D
                        string += coeff + "D" + exponent
                    if pair[0] < 0:
                        string += coeff + "S" + exponent
                    if pair[0] == 0:
                        string += coeff + "id"
        return string
    
    def __str__(self): #returns a string that represents self, expresed as a linear combination of powers of D and powers of H
        string=""
        for pair in self.pairs:
            coeff = pair[1]
            if (string != "") & (coeff > 0):# add plus sign if the next coefficient is positive, except for the first summand
                string += "+"
            if coeff < 0: # add minus sign in any case
                string += "-"
                coeff = abs(coeff)
            
            if coeff != 0:# omit any summands with coefficient 0
            
                exponent=abs(pair[0])
                if exponent==1: # omit exponent 1 from the notation
                    exponent = ""
                else:
                    exponent= ToExponent(exponent)
                
                if coeff==1: # omit coefficients 1 and -1 from the notation
                    coeff = ""
                else:
                    coeff = str(coeff) + "·"
                
                if pair[0] > 0:# powers of D
                    string += coeff + "D" + exponent
                if pair[0] < 0:
                    string += coeff + "S" + exponent
                if pair[0] == 0:
                    string += coeff + "id"
        return string
    
    def ToCob(self, sourceCLT, targetCLT, signed_lift=False, field=None):
        """Convert a BN-algebra morphism to a Cob morphism between the two
        given (1,3)-CLTs.

        When ``signed_lift`` is True and ``field`` is a prime > 1, each
        BN coefficient is lifted to its centered representative in
        (-p/2, p/2] before feeding into the Cob decos.  Default is the
        historical unsigned lift in [0, p); see OPEN_QUESTIONS.md item 1
        for the math question behind which lift is correct.

        Hardcoded to (1, 3)-tangles.  Callers passing CLTs of a different
        width will get a NotImplementedError — see OPEN_QUESTIONS.md item
        7 for the math deliverable needed to generalize to (1, 2k+1).
        """
        if not (sourceCLT.top == 1 and sourceCLT.bot == 3
                and targetCLT.top == 1 and targetCLT.bot == 3):
            raise NotImplementedError(
                "BNAlgebra.mor.ToCob is specific to (1, 3)-tangles "
                "(see OPEN_QUESTIONS.md item 7); got "
                "source=({}, {}), target=({}, {})".format(
                    sourceCLT.top, sourceCLT.bot,
                    targetCLT.top, targetCLT.bot))
        lift = (lambda c: signed_lift_coeff(c, field)) if signed_lift else (lambda c: c)
        decos = []

        for pair in self.pairs:
            c = lift(pair[1])
            if pair[0] > 0 :
                decos.append([pair[0]-1, 0, 1, c])
            elif pair[0] == 0:
                decos.append([0, 0, 0, c])
            elif pair[0] %2 == 0:
                power = int(Fraction(-1*pair[0], 2))
                # "corrected" flips the D-term sign (OPEN_QUESTIONS item 3).
                d_sign = ((-1) ** power) if _sign_convention == "corrected" else ((-1) ** (power - 1))
                decos.append([power, 0, 0, ((-1)** power) * c]) #(-H)^(n/2)
                decos.append([power - 1, 0, 1, d_sign * c]) #(-H)^(n/2-1) D
            elif pair[0] %2== 1:
                power = int(Fraction(-1*pair[0] -1, 2))
                decos.append([power, 0, ((-1)** power) * c ])#(-H)^((n-1)/2) S
            else:
                raise Exception("pair is not an integer?")
        mor = Cob.mor(sourceCLT, targetCLT, decos)
        mor.ReduceDecorations()
        return mor

class mor_alt(object):# work in progress
    """An element of Bar-Natan's algebra is a list of pairs [power,coeff]
    'power' is an integer, which determines the exponent of D (if positive) and the exponent of S (if negative)
    'coeff' is some non-zero integer (= coefficient in the base ring/field) # Alternatively, a Fraction object
    """
    import numpy as np # move this to the top of the file if used
    #__slots__ = 'pairs'
    
    def __init__(self,S,D,I):
        self.S = np.array(S) # list of coefficients
        self.D = np.array(D) # list of coefficients
        self.I = I #coefficient
    
    def simplify_mor(self,field):
        """simplify algebra elements by omitting superflous zeros."""
        self.S=[coeff_simplify(i,field) for i in self.S]
        self.D=[coeff_simplify(i,field) for i in self.D]
        self.I=coeff_simplify(self.I,field)
        while self.S[-1] ==0:
            del self.S[-1]
        while self.D[-1] ==0:
            del self.D[-1]
        return self
    
    def __add__(self, other):
        
        #newS = [a+b for a,b in zip_longest(self.S,other.S)]
        if len(self.S) < len(other.S):
            newS = other.S.copy()
            newS[:len(self.S)] += self.S
        else:
            newS = self.S.copy()
            newS[:len(other.S)] += other.S
        
        #newD = [a+b for a,b in zip_longest(self.D,other.D)]
        if len(self.D) < len(other.D):
            newD = other.D.copy()
            newD[:len(self.D)] += self.D
        else:
            newD = self.D.copy()
            newD[:len(other.D)] += other.D
        
        return mor(newS,newD,self.I+other.I).simplify_mor()

    def __mul__(self, other):
        newSmatrix = np.tensordot(self.S,other.S,axes=0)
        newS= [sum(np.diagonal(A[:, ::-1],len(other.S)-index)) for index in range(len(self.S)+len(other.S)-1)]
        
        newDmatrix = np.tensordot(self.D,other.D,axes=0)
        newD= [sum(np.diagonal(A[:, ::-1],len(other.D)-index)) for index in range(len(self.D)+len(other.D)-1)]
        
        return mor([[a1[0]+a2[0],a1[1]*a2[1]] for a1 in self.pairs for a2 in other.pairs if a1[0]*a2[0]>=0]).simplify_mor()
    
