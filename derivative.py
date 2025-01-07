from sympy import *
from variable import *

class Derivative(Symbol):
    UNDEFINED = -999
    def __init__(self, name, u_of_m = None, derivativeOrder = UNDEFINED, depVariable = None, indepVariable = None):
        super().__init__()

        self._name = name
        self._u_of_m = u_of_m
        self._derivativeOrder = derivativeOrder
        self._depVariable = depVariable
        self._indepVariable = indepVariable

        if derivativeOrder == Derivative.UNDEFINED or depVariable is None or indepVariable is None:
            self.guessDerivativeInfo()

        if u_of_m is None and self._depVariable is not None and self._indepVariable is not None:
            #assign based on the u_of_ms of the dep and indep vbles!
            if self._derivativeOrder == 1:
                self._u_of_m = UofM(self._depVariable._u_of_m._units / self._indepVariable._u_of_m._units)
            else:
                self._u_of_m = UofM(self._depVariable._u_of_m._units / (self._indepVariable._u_of_m._units)**self._derivativeOrder)




    def guessDerivativeInfo(self):
        ableToGuess = False
        if self._name.startswith('d'):
            d_count = self._name.count('d')
            two_count = self._name.count('2')
            three_count = self._name.count('3')
            if d_count == 2:
                if two_count == 0 and three_count == 0:
                    l1 = self._name.find('d')
                    l2 = self._name.find('d', l1 + 1)
                    self._derivativeOrder = 1
                    self._depVariable = Variable(self._name[l1 + 1:l2])
                    self._indepVariable = Variable(self._name[l2 + 1:])
                    ableToGuess = True
                elif two_count == 2:
                    l1 = self._name.find('d2')
                    l2 = self._name.find('d', l1 + 1)
                    l3 = self._name.find('2', l2 + 1)
                    self._derivativeOrder = 2
                    self._depVariable = Variable(self._name[l1 + 2:l2])
                    self._indepVariable = Variable(self._name[l2 + 1:l3])
                    ableToGuess = True
                elif three_count == 2:
                    l1 = self._name.find('d3')
                    l2 = self._name.find('d', l1 + 1)
                    l3 = self._name.find('3', l2 + 1)
                    self._derivativeOrder = 3
                    self._depVariable = Variable(self._name[l1 + 2:l2])
                    self._indepVariable = Variable(self._name[l2 + 1:l3])
                    ableToGuess = True

        return ableToGuess

    def __str__(self):
        return self._name

def derivatives(arg_string):
    syms =  symbols(arg_string)
    derivatives = []
    for sym in syms:
        derivative = Derivative(sym.name)
        derivatives.append(derivative)
    return derivatives
