from sympy import *
class BaseUnit(Symbol):

    def __init__(self, sym):
        super().__init__()
        self._sym = Symbol(sym)

    def __str__(self):
        #return self._sym
        return self.name

def base_units(arg_string):
    syms = symbols(arg_string)
    base_units = []
    for sym in syms:
        bu = BaseUnit(sym.name)
        base_units.append(bu)
    return base_units
