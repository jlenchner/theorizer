"""The base Derivatif class, used for creating and working with ordinary derivatives. (The
    funny class name is due to the fact that the class Derivative is used in sympy to handle symbolic
    differentiation, and the name Differential causes a namespace conflict. Derivatifs are similar to
    Variables but they contain three additional fields: _derivativeOrder (the order of the derivative),
    _depVariable (the dependent derivative) and _indepVariable (the independent variable).  In
    many cases, especially in physics, the system can correctly guess all values from the associated
    name. Single Derivatifs can be initialized via the constructor just supplying the name as a string,
    e.g., dxdt = Variable('dxdt'), which will guess the _u_of_m UofM(UofM.m/UofM.s), _derivativeOrder=1,
    _depVariable = x (=Variable('x')), _indepVariable = t (=Variable('t')).

    If using guessed values for the fields _u_of_m, _derivativeOrder, _depVariable and _ubdepVariable,
    Derivatifs can be initialized in batch using the method derivatives(), e.g.,
    dx1dt, d2x1dt2, dx2dt, d2x2dt2  = derivatives('dx1dt,d2x1dt2,dx2dt,d2x2dt2'), where the first two
    Derivatifs will be guessed to be first order derivatives, and the second two  will be guessed to be
    second order derivatives. Note that the derivatives() method lives outside of the Derivatif class,
    in keeping with how symbols() is used in the base sympy package.
"""

# Author: Jonathan Lenchner (lenchner@us.ibm.com)
#
# License: BSD 3-Clause

import sympy
from sympy import *
from variable import *

class Derivatif(Symbol):
    UNDEFINED = -999

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, args[0])

    def __init__(self, name, u_of_m = None, derivativeOrder = UNDEFINED, depVariable = None, indepVariable = None):
        """
        :param name: Mandatory and should be specified as a string, typically
                     the string equivalent of the Derivatif symbols
                     (e.g., for a Derivatif dxdt, the name would typically be 'dxdt')
        :param u_of_m: A fully constructed UofM, e.g., UofM(UofM.m/UofM.s) for a Derivatif
                     specifying velocity such as dxdt. Omit if parameters other than the name
                     are to be guess/inferred from the name.
        :param derivativeOrder: An integer specifying the order of the Derifatif, typically 1 or 2,
                     and occasionally 3. Omit if parameters other than the name
                     are to be guess/inferred from the name.
        :param depVariable:  The dependent Variable associated with the Derivatif, e.g., for the Derivatif
                     dxdt, one would likely specify depVariable = Variable('x'). Omit if parameters other
                     than the name are to be guess/inferred from the name.
        :param indepVariable: The independent Variable associated with the Derivatif, e.g., for the Derivatif
                     dxdt, one would likely specify indepVariable = Variable('t'). Omit if parameters other
                     than the name are to be guess/inferred from the name.
        """

        super().__init__()

        self._name = name
        self._u_of_m = u_of_m
        self._derivativeOrder = derivativeOrder
        self._depVariable = depVariable
        self._indepVariable = indepVariable


        if derivativeOrder == Derivatif.UNDEFINED or depVariable is None or indepVariable is None:
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
                if two_count <= 1 and three_count <= 1:
                    l1 = self._name.find('d')
                    l2 = self._name.find('d', l1 + 1)
                    self._derivativeOrder = 1
                    self._depVariable = Variable(self._name[l1 + 1:l2])
                    self._indepVariable = Variable(self._name[l2 + 1:])
                    ableToGuess = True
                elif two_count <= 3:
                    l1 = self._name.find('d2')
                    l2 = self._name.find('d', l1 + 1)
                    l3 = self._name.find('2', l2 + 1)
                    self._derivativeOrder = 2
                    self._depVariable = Variable(self._name[l1 + 2:l2])
                    self._indepVariable = Variable(self._name[l2 + 1:l3])
                    ableToGuess = True
                elif three_count <= 3:
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
    """The usual batch constructor when relying on guessed values for the parameters other than _name.
            Specify a single string of comma delimited names, dx1dt, d2x1dt2, dx2dt, d2x2dt2  =
            derivatives('dx1dt,d2x1dt2,dx2dt,d2x2dt2'), which will create separate Derivatif objects
            for dx1dt, d2x1dt2, dx2dt and d2x2dt2 with the expected orders, dependent and independent
            Derivatifs.
    """

    syms =  symbols(arg_string)
    derivatives = []
    for sym in syms:
        derivative = Derivatif(sym.name)
        derivatives.append(derivative)
    return derivatives
