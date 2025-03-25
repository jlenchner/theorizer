"""The base Variable class. In the case of most physics equations, the system can correctly
    guess the associated unit of measure. Single variables can be initialized via the constructor
    just supplying the name of the variable as a string, e.g., v = Variable('v'), which will guess
    the unit of measure (class UofM) m/s for velocity. Any Variable constructed in like manner and
    beginning with 'v' will result in the same guess. See the UofM class for variables whose units
    of measure can be guessed in this manner. Of course, a fully constructed UofM can also be passed.
    Thus v = Variable('v') is equivalent to v = Variable('v', UofM(UofM.m/UofM.s)).

    If using guessed UofMs, Variables can be initialized in batch using the method variablees(), e.g.,
    d1, d2, m1, m2, Fg = variables('d1, d2, m1, m2, Fg'), where the first two variables would be guessed
    to be distances (UofM.m), the second two would be guessed to be masses (UofM.kg), and the last one
    would be guess to be a mass (UofM.kg*UofM.m/(UofM.s**2). Note that the variablees() method lives
    outside of the Variable class, in keeping with how symbols() is used in the base sympy package.
"""

# Author: Jonathan Lenchner (lenchner@us.ibm.com)
#
# License: BSD 3-Clause

from unitOfMeasure import *
from sympy import *


class Variable(Symbol):
    VAR_TO_UofM_DICT = dict()

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, args[0])

    def __init__(self, name, u_of_m = None):
        """
        :param name: Mandatory and should be specified as a string, typically
                     the string equivalent of the variable symbols
                     (e.g., v = Variable('v'))
        :param u_of_m: Fully constructed UofM. Omit if the UofM should be inferred
                    from the name. Example of use with fully constructed UofM:
                    v = Variable('v', UofM(UofM.m/UofM.s))
        """
        super().__init__()
        self._name = name
        if u_of_m is None:
            self._u_of_m = UofM.GuessUofM(measuredQtyDesc=self._name)
        else:
            self._u_of_m = u_of_m
        Variable.VAR_TO_UofM_DICT.update({name: self._u_of_m})


    @classmethod
    def SetAll(cls, var_to_UofM_dict):
        Variable.VAR_TO_UofM_DICT = var_to_UofM_dict

    @classmethod
    def GetAll(cls):
        allVars = []
        keys = Variable.VAR_TO_UofM_DICT.keys()
        for key in keys:
            var = Variable(key, Variable.VAR_TO_UofM_DICT.get(key))
            allVars.append(var)

        return allVars

    def isDimensionless(self):
        """This generally only applies to Constants (which are Variables),
            though in theory can apply to any Variable
        """

        return self._u_of_m is not None and self._u_of_m._units == 1


    def __str__(self):
        return self._name


def variables(arg_string):
    """The usual batch constructor when relying on guessed UofMs. Specify a single string
        of comma delimited names, e.g., d1, d2, m1, m2, Fg = variables('d1, d2, m1, m2, Fg'),
        which will create separate Variable objects for d1, d2, m1, m2 and Fg.
    """

    syms =  symbols(arg_string)
    vars = []
    for sym in syms:
        var = Variable(sym.name)
        vars.append(var)
    return vars
