import random

from sympy import *
from equation import Equation

class EquationSystem:

    def __init__(self, vars = [], measuredVars =[],  equations = [], maxVarsPerEqn=Equation.NO_MAX):
        self._equations = equations
        self._vars = vars
        if len(vars) == 0 and len(equations) > 0:
            vars_as_set = set()
            for eqn in equations:
                vars_as_set.update(eqn._poly.free_symbols)
            for var in vars_as_set:   #Turn set into an array
                self._vars.append(var)

        self._measuredVars = measuredVars
        if len(measuredVars) == 0 and len(vars) > 0:
            self._measuredVars = list(vars)[:len(vars)//2]

        self._maxVarsPerEqn = maxVarsPerEqn
        self._nonMeasuredVars = set(self._vars) - set(self._measuredVars)


    def getZeroes(self):   #not yet implemented
        return None

    def getVarNames(self):
        var_names = []
        for var in self._vars:
           var_names.append(str(var))

        return var_names
    
    def getMeasuredVars(self):
        measured_var_names = []
        for var in self._measuredVars:
            measured_var_names.append(str(var))
        
        return measured_var_names
    
    def getNonMeasuredVars(self):
        non_measured_var_names = []
        for var in self._nonMeasuredVars:
            non_measured_var_names.append(str(var))
        
        return non_measured_var_names
    
    def getEquations(self):
        equation_strings = []
        for eqn in self._equations:
            equation_strings.append(eqn.toString())
            
        return equation_strings

    @classmethod
    def GenerateRandom(cls, vars, measuredVars, numEqns, maxVarsPerEqn):
        eqns = []

        while True:  #A crude way to make sure we use all the variables
            vars_used = set()
            for i in range(numEqns):
                eqn = Equation.GenerateRandom(vars=vars, max_vars=maxVarsPerEqn)
                vars_used = vars_used.union(eqn.getVarsUsed())
                eqns.append(eqn)

            if len(vars_used) < len(vars):
                eqns = []
            else:
                break

        return EquationSystem(vars=vars, measuredVars= measuredVars, equations=eqns, maxVarsPerEqn=maxVarsPerEqn)

    def replaceRandomEqnByIndex(self, eqnIndex):
        eqn = None
        vars_used = set()
        for i in range(len(self._equations)):
            if i != eqnIndex:
                vars_used = vars_used.union(self._equations[i].getVarsUsed())

        while True:
            eqn = Equation.GenerateRandom(vars=self._vars, max_vars=self._maxVarsPerEqn)
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
        exp += '\n'

        measuredVarNames = self.getMeasuredVars()
        exp = "Measured Variables: \n"
        for i in range(len(measuredVarNames)):
            exp += measuredVarNames[i]
            if i < len(measuredVarNames) - 1:
                exp += ','

        exp += '\n'
        exp += '\n'
        exp += "Non Measured Variables: \n"
        nonMeasuredVarNames = self.getNonMeasuredVars()
        for i in range(len(nonMeasuredVarNames)):
            exp += nonMeasuredVarNames[i]
            if i < len(nonMeasuredVarNames) - 1:
                exp += ','

        exp += '\n'
        exp += '\n'
        exp += "Equations: \n"
        for i in range(len(self._equations)):
            exp += self._equations[i].toString()
            if i < len(self._equations) - 1:
                exp += '\n'

        return exp









