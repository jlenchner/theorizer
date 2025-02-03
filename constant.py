from unitOfMeasure import *
from variable import *
from sympy import *


class Constant(Variable):  #should be able to automatically pick up units of measure for c, G, h (Plank's constant), h_bar
                           # and others.....
    CONSTANT_TO_UofM_DICT = dict()
    UNDEFINED = -999

    def __init__(self, name, u_of_m=1, value=UNDEFINED):  #in other words, the default is that the constant is dimensionless
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
        assumed_u_of_m = UofM.GuessUofM(measuredQtyDesc=self._name)
        if self._u_of_m == assumed_u_of_m:
            if self._name == 'c':
                self._value = 2.998e8
            elif self._name == 'G':
                self._value = 6.6743e-11
            elif self._name == 'h':
                self._value = 6.62607015e-34
            elif self._name == 'h-bar':
                self._value = 1.0545771817e-34

        return self._value

    def __str__(self):
        return self._name


def constants(arg_string):
    syms =  symbols(arg_string)
    cons = []
    for sym in syms:
        con = Constant(sym.name)
        cons.append(con)
    return cons
