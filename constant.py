"""The base Constant class. Constants are derived from the Variable class and are in all ways
    identical to Variables except that they take fixed values, stored in the field _value,
    and therefore are regarded as constant Variables. Additional common Constants can be added
    to the UofM class in unitOfMeasure.py: see the parallel lists (at class level) ALL_MEASURED_QUANTITIES
    and  ALL_UNITS. Note that mathematical constants, and even some physical constants are
    dimensionless, in which case the _u_of_m can be specified as UofM(1).
"""

# Author: Jonathan Lenchner (lenchner@us.ibm.com)
#
# License: BSD 3-Clause

from unitOfMeasure import *
from variable import *
from sympy import *


class Constant(Variable):  #should be able to automatically pick up units of measure for c, G, h (Plank's constant), h_bar
                           # and others.....
    CONSTANT_TO_UofM_DICT = dict()
    UNDEFINED = -999

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, args[0])

    def __init__(self, name, u_of_m=1, value=UNDEFINED):  #in other words, the default is that the constant is dimensionless
        """
            :param name: Mandatory and should be specified as a string, typically
                         the string equivalent of the Constant symbol
            :param u_of_m: A fully constructed UofM. For dimensionless UofMs, pass UofM(1).
            :param value: The value of the constant in the units specified in u_of_m. For
                         values of predefined constants refer to the guessValue() method
                         of this class.
         """

        self._name = name
        if u_of_m == 1:
            guess = UofM.GuessUofM(measuredQtyDesc=self._name)
            if guess is not None:
                self._u_of_m = guess
            else:
                self._u_of_m = 1
        else:
            self._u_of_m = u_of_m
        if value == Constant.UNDEFINED and self._u_of_m is not None:
            self._value = Constant.UNDEFINED  #just to seed it with something
            self.guessValue()
        else:
            self._value = value

        Constant.CONSTANT_TO_UofM_DICT.update({name: self._u_of_m})

    def guessValue(self):
        """
            Guesses the value associated with the constant based on the name. When a new constant is added
            the associated symbol as a strong should be added to UofM.ALL_MEASURED_QUANTITIES, the associated
            UofM should be added to UofM.ALL_UNITS and the value (in the supplied units) should be added below.
        """
        assumed_u_of_m = UofM.GuessUofM(measuredQtyDesc=self._name)
        if self._u_of_m == assumed_u_of_m:
            if self._name == 'c':
                self._value = 2.99792458e8
            elif self._name == 'G':
                self._value = 6.6743e-11
            elif self._name == 'h':
                self._value = 6.62607015e-34
            elif self._name == 'h-bar':
                self._value = 1.0545771817e-34
            elif self._name == 'k':
                self._value = 1.380649e-23
            elif self._name == "e":
                self._value = 1.602176634e-19
            elif self._name == "pi":
                self._value = 3.14159265359

        return self._value


    def __str__(self):
        return self._name


def constants(arg_string):
    """A batch constructor that can be use when relying on guessed UofMs and values. Specify a single
           string of comma delimited names, e.g., G, c, pi = constants('G,c,pi'),
           which will create separate Constant objects for G, c and pi.
       """
    syms =  symbols(arg_string)
    cons = []
    for sym in syms:
        con = Constant(sym.name)
        cons.append(con)
    return cons
