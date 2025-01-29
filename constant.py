from unitOfMeasure import *
from variable import *
from sympy import *


class Constant(Variable):  #should be able to automatically pick up units of measure for c, G, h (Plank's constant), h_bar
                           # and others.....
    CONSTANT_TO_UofM_DICT = dict()

    def __init__(self, name, u_of_m = 1):  #in other words, the default is that the constant is dimensionless
        self._name = name
        if u_of_m == 1:
            guess = UofM.GuessUofM(measuredQtyDesc=self._name)
            if guess is not None:
                self._u_of_m = guess
            else:
                self._u_of_m = 1
        else:
            self._u_of_m = u_of_m
        Constant.CONSTANT_TO_UofM_DICT.update({name: self._u_of_m})


    def __str__(self):
        return self._name


def constants(arg_string):
    syms =  symbols(arg_string)
    cons = []
    for sym in syms:
        con = Constant(sym.name)
        cons.append(con)
    return cons
