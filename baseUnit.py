""""The BaseUnit class for managing primitive (non-derived) Units of Measure,
    like kg, m, s. For a complete list of BaseUnits refer to the UofM class in
    the file unitofMeasure.py. This file should not be modified.
"""

# Author: Jonathan Lenchner (lenchner@us.ibm.com)
#
# License: BSD 3-Clause

from sympy import *
class BaseUnit(Symbol):

    def __init__(self, sym):
        super().__init__()
        self._sym = Symbol(sym)

    def __str__(self):
        #return self._sym
        return self.name

def base_units(arg_string):
    """A batch constructor. Specify a single string of comma delimited names, e.g.,
       m,kg,s,mol,A,cd,K = base_units('m,kg,s,mol,A,cd,K'), which you can see done
       in the UofM class in the file unitOfMeasure.py.
    """

    syms = symbols(arg_string)
    base_units = []
    for sym in syms:
        bu = BaseUnit(sym.name)
        base_units.append(bu)
    return base_units
