import random

from sympy import *
from polynomial import Polynomial

class EquationSystem:

    def __init__(self, vars = [], equations = [], maxVarsPerEqn=Polynomial.NO_MAX):
        self._equations = equations
        self._vars = vars
        if len(vars) == 0 and len(equations) > 0:
            vars_as_set = set()
            for eqn in equations:
                vars_as_set.update(eqn._poly.free_symbols)
            for var in vars_as_set:   #Turn set into an array
                self._vars.append(var)
        self._maxVarsPerEqn = maxVarsPerEqn

    def getZeroes(self):   #not yet implemented
        return None

    def getVarNames(self):
        var_names = []
        for var in self._vars:
           var_names.append(str(var))
        var_names.sort()

        return var_names

    @classmethod
    def GenerateRandom(cls, vars, numEqns, maxVarsPerEqn):
        eqns = []

        while True:  #A crude way to make sure we use all the variables
            vars_used = set()
            for i in range(numEqns):
                eqn = Polynomial.GenerateRandom(vars=vars, max_vars=maxVarsPerEqn)
                vars_used = vars_used.union(eqn.getVarsUsed())
                eqns.append(eqn)

            if len(vars_used) < len(vars):
                eqns = []
            else:
                break

        return EquationSystem(vars=vars, equations=eqns, maxVarsPerEqn=maxVarsPerEqn)

    def replaceRandomEqnByIndex(self, eqnIndex):
        eqn = None
        vars_used = set()
        for i in range(len(self._equations)):
            if i != eqnIndex:
                vars_used = vars_used.union(self._equations[i].getVarsUsed())

        while True:
            eqn = Polynomial.GenerateRandom(vars=self._vars, max_vars=self._maxVarsPerEqn)
            all_vars = vars_used.union(eqn.getVarsUsed())
            if len(all_vars) == len(self._vars):
                self._equations[eqnIndex] = eqn
                break

        return eqn

    def replaceRamdomEqn(self):
        randIndex = random.randint(0, len(self._equations) - 1)
        return self.replaceRandomEqnByIndex(randIndex)

    def toString(self):
        varNames = self.getVarNames()

        exp = "Variables: \n"
        for i in range(len(varNames)):
            exp += varNames[i]
            if i < len(varNames) - 1:
                exp += ','
        exp += '\n'
        exp += "Equations: \n"
        for i in range(len(self._equations)):
            exp += self._equations[i].toString()
            if i < len(self._equations) - 1:
                exp += '\n'

        return exp









