import random
import os
from variable import *

from sympy import *
from equation import Equation
from m2_functions import *

class EquationSystem:

    def __init__(self, vars = [],derivatives = [],  measuredVars =[],  equations = [], max_var_and_derivatives_per_eqn=Equation.NO_MAX):
        self._equations = equations
        self._vars = vars
        if len(vars) == 0 and len(equations) > 0:
            vars_as_set = set()
            for eqn in equations:
                vars_in_eqn = set()
                for sym in eqn._poly.free_symbols:
                    if isinstance(sym, Variable):
                        vars_in_eqn.add(sym)
                vars_as_set.update(vars_in_eqn)
            for var in vars_as_set:   #Turn set into an array
                self._vars.append(var)
        self._derivatives = derivatives
        if len(derivatives) == 0 and len(equations) > 0:
            derivs_as_set = set()
            for eqn in equations:
                derivs_in_eqn = set()
                for sym in eqn._poly.free_symbols:
                    if isinstance(sym, Derivative):
                        derivs_in_eqn.add(sym)
                derivs_as_set.update(derivs_in_eqn)
            for deriv in derivs_as_set:  # Turn set into an array
                self._derivatives.append(deriv)

        self._measuredVars = measuredVars
        if len(measuredVars) == 0 and len(vars) > 0:
            self._measuredVars = list(vars)[:len(vars)//2]

        self.max_var_and_derivatives_per_eqn = max_var_and_derivatives_per_eqn
        self._nonMeasuredVars = set(self._vars) - set(self._measuredVars)


    def getZeroes(self):   #not yet implemented
        return None

    def isConsistent(self):
        """
        Function that checks the consistency of an equation system.
        Inputs:
            eqnSystem: the system of equations to check
        Outputs:
            True if consistent, False otherwise
        """

        # Extracting the variables and equations from the system
        variables = self.getVarNames()
        # print("Variables: " + str(variables))
        equations = self.getEquations()
        # print("Equations: " + str(equations))

        temp_filename = "temp_results.txt"
        check_consistency(variables, equations, temp_filename)

        # Read the results from the file
        with open(temp_filename, "r") as file:
            result = file.read().strip()

        # Remove the temporary file
        os.remove(temp_filename)

        if "inconsistent" in result:
            return False
        else:
            return True

    def checkConsistencyReplacedSystems(self, index, num):
        """
        Checks the consistency of the system of equations after replacing a random equation
        num number of times.
        Inputs:
            eqnSystem: the system of equations to check
            index: the index of the equation that was replaced
            num: the number of times to replace the equation
        Outputs:
            numConsistent: the number of consistent systems found
        """

        numConsistent = 0
        for i in range(num):
            eqn = self.replaceRandomEqnByIndex(index)
            # print(f"Replaced equation {index} with: {str(eqn)}")
            # print("Current System: \n" + str(eqnSystem))
            if self.isConsistent():
                # print("Consistent \n")
                numConsistent += 1
        return numConsistent

    def project(self):
        """
        Projects the equations onto the measured variables.
        Inputs:
            eqnSystem: the system of equations to project
        Outputs:
            projectedEquations: the projected equations
        """
        # Extracting the variables and equations from the system
        variables = self.getVarNames()
        equations = self.getEquations()
        measured_vars = self.getMeasuredVars()
        non_measured_variables = self.getNonMeasuredVars()

        temp_filename = "temp.txt"

        projection(variables, equations, measured_vars, non_measured_variables, temp_filename)

        # Read the results from the file
        with open(temp_filename, "r") as file:
            result = file.read().strip()

        # Remove the temporary file
        os.remove(temp_filename)

        with open("projection_output.txt", "a") as file:
            file.write(str(self) + "\n")
            file.write(result + "\n")
            file.write("\n")

        return result

    def projectUnknownMeasuredVars(self):
        """
        Projects the equations onto the measured variables.
        Inputs:
            eqnSystem: the system of equations to project
        Outputs:
            projectedEquations: the projected equations
        """
        # Extracting the variables and equations from the system
        variables = self.getVarNames()
        equations = self.getEquations()
        measured_vars = self.getMeasuredVars()

        temp_filename = "temp.txt"
        result = ""

        for i in range(1, len(measured_vars)):
            print(f"Measuring {measured_vars[0:i]}")
            non_measured_variables = [var for var in variables if var not in measured_vars[0:i]]

            projection(variables, equations, measured_vars[0:i], non_measured_variables, temp_filename)

            # Read the results from the file
            with open(temp_filename, "r") as file:
                result = file.read().strip()
            basis = result.split("\n")[-1]

            # Remove the temporary file
            os.remove(temp_filename)

            if basis != "Polynomials of the Gröbner basis of the eliminated ideal:":
                # Add result to a file called projection_output.txt
                with open("unknown projection_output.txt", "a") as file:
                    file.write(str(self) + "\n")
                    file.write(result + "\n")
                    file.write("\n")
                break

        return result

    @classmethod
    def ProjectRandomSystems(cls, vars, derivatives, measured_vars, num):
        """
        Projects random systems of equations onto the measured variables.
        Inputs:
            num: the number of random systems to project
        Outputs:
            results: a list of results from the projections
        """
        results = []
        for i in range(num):
            eqnSystem = EquationSystem.GenerateRandom(vars, derivatives, measured_vars, 4, 6)
            print("System: \n" + str(eqnSystem))
            result = eqnSystem.projectUnknownMeasuredVars()
            results.append(result)
        return results

    @classmethod
    def CheckConsistencyRandomSystems(cls, vars, num):
        """
        Checks the consistency of random systems of equations.
        Inputs:
            num: the number of random systems to check
        Outputs:
            numConsistent: the number of consistent systems found
        """
        numConsistent = 0
        for i in range(num):
            eqnSystem = EquationSystem.GenerateRandom(vars, 4, 6)
            if eqnSystem.isConsistent():
                numConsistent += 1
        return numConsistent

    def getVarNames(self):
        var_names = []
        for var in self._vars:
           var_names.append(str(var))

        return var_names

    def getDerivNames(self):
        deriv_names = []
        for derirvative in self._derivatives:
           deriv_names.append(str(derirvative))

        return deriv_names
    
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
            equation_strings.append(str(eqn))
            
        return equation_strings

    @classmethod
    def GenerateRandom(cls, vars, derivatives, measuredVars, numEqns, max_var_and_derivatives_per_eqn):
        eqns = []

        while True:  #A crude way to make sure we use all the variables
            symbols_used = set()
            for i in range(numEqns):
                eqn = Equation.GenerateRandom(vars=vars, derivatives=derivatives, max_var_and_derivatives_per_eqn=max_var_and_derivatives_per_eqn)
                symbols_used = symbols_used.union(eqn.getSymbolsUsed())
                eqns.append(eqn)

            if len(symbols_used) < len(vars):
                eqns = []
            else:
                break

        return EquationSystem(vars=vars, measuredVars= measuredVars, equations=eqns, max_var_and_derivatives_per_eqn=max_var_and_derivatives_per_eqn)

    def replaceRandomEqnByIndex(self, eqnIndex):
        eqn = None
        symbols_used = set()
        for i in range(len(self._equations)):
            if i != eqnIndex:
                vars_used = symbols_used.union(self._equations[i].getSymbolsUsed())

        while True:
            eqn = Equation.GenerateRandom(vars=self._vars, max_vars=self._maxVarsPerEqn)
            all_symbols = symbols_used.union(eqn.getSymbolsUsed())
            if len(all_symbols) == len(self._vars):
                self._equations[eqnIndex] = eqn
                break

        return eqn

    def replaceRamdomEqn(self):
        randIndex = random.randint(0, len(self._equations) - 1)
        return self.replaceRandomEqnByIndex(randIndex)

    def __str__(self):
        varNames = self.getVarNames()
        exp = "Variables: \n"
        for i in range(len(varNames)):
            exp += varNames[i]
            if i < len(varNames) - 1:
                exp += ','
        exp += '\n'
        exp += '\n'

        derivNames = self.getDerivNames()
        exp += "Derivatives: \n"
        for i in range(len(derivNames)):
            exp += derivNames[i]
            if i < len(derivNames) - 1:
                exp += ','
        exp += '\n'
        exp += '\n'

        measuredVarNames = self.getMeasuredVars()
        exp += "Measured Variables: \n"
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
            exp += str(self._equations[i])
            if i < len(self._equations) - 1:
                exp += '\n'

        return exp





































