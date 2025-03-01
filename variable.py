from unitOfMeasure import *
from sympy import *


class Variable(Symbol):
    VAR_TO_UofM_DICT = dict()

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, args[0])

    def __init__(self, name, u_of_m = None):
        super().__init__()
        self._name = name
        if u_of_m is None:
            self._u_of_m = UofM.GuessUofM(measuredQtyDesc=self._name)
        else:
            self._u_of_m = u_of_m
        Variable.VAR_TO_UofM_DICT.update({name: self._u_of_m})


    @classmethod
    def SetAll(cls, var_to_UofM_dict):  #perhaps do this with parallel arrays
        Variable.VAR_TO_UofM_DICT = var_to_UofM_dict

    @classmethod
    def GetAll(cls):
        allVars = []
        keys = Variable.VAR_TO_UofM_DICT.keys()
        for key in keys:
            var = Variable(key, Variable.VAR_TO_UofM_DICT.get(key))
            allVars.append(var)

        return allVars

    def isDimensionless(self):  #this generally only applies to Constants (which are Variable),
                                # though in theory can apply to any Variable
        return self._u_of_m is not None and self._u_of_m._units == 1


    def __str__(self):
        return self._name


def variables(arg_string):
    syms =  symbols(arg_string)
    vars = []
    for sym in syms:
        var = Variable(sym.name)
        vars.append(var)
    return vars
