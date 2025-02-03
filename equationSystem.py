import random
import os
from variable import *
from constant import *
from derivative import *

from sympy import *
from equation import Equation
from m2_functions import *

class EquationSystem:
    _LastVarsDerivsAndConstants = []
    _LookupDict = None

    def __init__(self, vars = [],derivatives = [],  constants = [], measuredVars =[],  equations = [], max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX):
        self._equations = equations
        self._vars = vars
        if len(vars) == 0 and len(equations) > 0:
            vars_as_set = set()
            for eqn in equations:
                vars_in_eqn = set()
                for sym in eqn._poly.free_symbols:
                    if isinstance(sym, Variable) and not isinstance(sym, Constant):
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
                    if isinstance(sym, Derivatif):
                        derivs_in_eqn.add(sym)
                derivs_as_set.update(derivs_in_eqn)
            for deriv in derivs_as_set:  # Turn set into an array
                self._derivatives.append(deriv)
        self._constants = constants
        if len(constants) == 0 and len(equations) > 0:
            constants_as_set = set()
            for eqn in equations:
                constants_in_eqn = set()
                for sym in eqn._poly.free_symbols:
                    if isinstance(sym, Constant):
                        constants_in_eqn.add(sym)
                constants_as_set.update(constants_in_eqn)
            for const in constants_as_set:  # Turn set into an array
                self._constants.append(const)


        self._measuredVars = measuredVars
        if len(measuredVars) == 0 and len(vars) > 0:
            self._measuredVars = list(vars)[:len(vars)//2]

        self._max_vars_derivatives_and_constants_per_eqn = max_vars_derivatives_and_constants_per_eqn
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

    def getConstantNames(self):
        constant_names = []
        for constant in self._constants:
           constant_names.append(str(constant))

        return constant_names
    
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
    def GenerateRandom(cls, vars, derivatives, constants, measuredVars, numEqns, max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX):
        eqns = []

        while True:  #A crude way to make sure we use all the variables
            symbols_used = set()
            for i in range(numEqns):
                eqn = Equation.GenerateRandom(vars=vars, derivatives=derivatives, constants=constants,
                                              max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
                symbols_used = symbols_used.union(eqn.getSymbolsUsed())
                eqns.append(eqn)

            if len(symbols_used) < len(vars):
                eqns = []
            else:
                break

        return EquationSystem(vars=vars, derivatives=derivatives, constants=constants, measuredVars= measuredVars, equations=eqns,
                              max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)



    @classmethod
    def GenerateRandomDimensionallyConsistent(cls, vars, derivatives, constants, measuredVars, numEqns, max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX):
        eqns = []

        varsDerivsAndConstants = []
        varsDerivsAndConstants.extend(vars)
        varsDerivsAndConstants.extend(derivatives)
        varsDerivsAndConstants.extend(constants)

        varsDerivsAndConstantsRemainingToBeUsed = set(varsDerivsAndConstants)

        max_power = EquationSystem.DetermineMaxPower(vars, derivatives, constants)

        if varsDerivsAndConstants != EquationSystem._LastVarsDerivsAndConstants:
            Equation._logger.info("Using NEW lookup dictionary.")
            EquationSystem._LastVarsDerivsAndConstants = varsDerivsAndConstants
            EquationSystem._LookupDict = Equation.GetUofMToPrimitiveTermLookupTable(vars, derivatives, constants,
                                                    max_power=max_power,
                                                    max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
        else:
            Equation._logger.info("Using existing lookup dictionary.")


        symbols_used = set()
        u_of_ms = set()   #each equation should have a distinct u_of_m
        eqn = Equation.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivatives, constants=constants,
                                    u_of_mToTermLookupDict=EquationSystem._LookupDict, max_power=max_power,
                                    max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
        firstTerm = eqn.getTerms()[0]
        u_of_m = Equation.GetUofMForTerm(firstTerm)
        u_of_ms.add(u_of_m)
        symbols_used = eqn.getSymbolsUsed()
        eqns.append(eqn)

        Equation._logger.info("Eqn: " + str(eqn))

        #Note that additional equations should probably have different u_of_ms....

        varsDerivsAndConstantsRemainingToBeUsed = varsDerivsAndConstantsRemainingToBeUsed - symbols_used
        while len(varsDerivsAndConstantsRemainingToBeUsed) > 0 and len(eqns) < numEqns:
            nextVarDerivativeOrConstant = next(iter(varsDerivsAndConstantsRemainingToBeUsed))
            next_var = None
            next_derivative = None
            next_constant = None
            if isinstance(nextVarDerivativeOrConstant, Constant):
                next_constant = nextVarDerivativeOrConstant
            elif isinstance(nextVarDerivativeOrConstant, Variable):
                next_var = nextVarDerivativeOrConstant
            elif isinstance(nextVarDerivativeOrConstant, Derivatif):
                next_derivative = nextVarDerivativeOrConstant
            next_eqn = Equation.GenerateRandomDimensionallyConsistentEquationWithSpecifiedVarOrDerivative(vars=vars, derivatives=derivatives,
                                        constants=constants, u_of_mToTermLookupDict=EquationSystem._LookupDict,
                                        given_var=next_var, given_derivative=next_derivative, given_constant=next_constant,
                                        max_power = max_power,
                                        max_vars_derivatives_and_constants_per_eqn = max_vars_derivatives_and_constants_per_eqn)

            if next_eqn is not None:
                eqns.append(next_eqn)
                Equation._logger.info("Eqn: " + str(next_eqn))
                symbols_used = next_eqn.getSymbolsUsed()
                varsDerivsAndConstantsRemainingToBeUsed = varsDerivsAndConstantsRemainingToBeUsed - symbols_used
            else:
                varsDerivsAndConstantsRemainingToBeUsed.remove(nextVarDerivativeOrConstant)
                Equation._logger.info("Dropping var-deriv-or-constant: " + str(nextVarDerivativeOrConstant))

        while len(eqns) < numEqns:
            u_of_ms_in_use = set()
            for eqn in eqns:
                firstTerm = eqn.getTerms()[0]
                u_of_m_for_term = Equation.GetUofMForTerm(firstTerm)
                u_of_ms_in_use.add(u_of_m_for_term)
            while True:
                eqn = Equation.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivatives, constants=constants,
                                                            u_of_mToTermLookupDict=EquationSystem._LookupDict, max_power=max_power,
                                                            max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
                firstTerm = eqn.getTerms()[0]
                u_of_m_for_term = Equation.GetUofMForTerm(firstTerm)
                if u_of_m_for_term not in u_of_ms_in_use:
                    eqns.append(eqn)
                    Equation._logger.info("Eqn: " + str(eqn))
                    break

        return EquationSystem(vars=vars, derivatives=derivatives, constants=constants, measuredVars=measuredVars, equations=eqns,
                              max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
    @classmethod
    def DetermineMaxPower(cls, vars, derivatives, constants):
        varsDerivsAndConstants = []
        varsDerivsAndConstants.extend(vars)
        varsDerivsAndConstants.extend(derivatives)
        varsDerivsAndConstants.extend(constants)

        if len(varsDerivsAndConstants) > 10:
            return 2
        else:
            return 3

    def replaceRandomDimensionallyConsistentEqnByIndex(self, eqnIndex): #To be fully implemented
        #First make sure that all vars, derivatives and constants have u_of_ms!
        for var in self._vars:
            if var._u_of_m is None:
                Equation._logger.error("The variable " + str(var) + " has no _u_of_m!")
                return self._equations[eqnIndex]
        for deriv in self._derivatives:
            if deriv._u_of_m is None:
                Equation._logger.error("The derivative " + str(deriv) + " has no _u_of_m!")
                return self._equations[eqnIndex]
        for const in self._constants:
            if const._u_of_m is None:
                Equation._logger.error("The constant " + str(const) + " has no _u_of_m!")
                return self._equations[eqnIndex]

        varsDerivsAndConstants = []
        varsDerivsAndConstants.extend(self._vars)
        varsDerivsAndConstants.extend(self._derivatives)
        varsDerivsAndConstants.extend(self._constants)

        max_power = EquationSystem.DetermineMaxPower(self._vars, self._derivatives, self._constants)

        if varsDerivsAndConstants != EquationSystem._LastVarsDerivsAndConstants:
            EquationSystem._LastVarsDerivsAndConstants = varsDerivsAndConstants
            EquationSystem._LookupDict = Equation.GetUofMToPrimitiveTermLookupTable(vars, derivatives, constants,
                                max_power=max_power,
                                max_vars_derivatives_and_constants_per_eqn=self._max_vars_derivatives_and_constants_per_eqn)

        symbols_of_others = set()
        symbols_of_this_eqn = set()
        for i in range(len(self._equations)):
            if i != eqnIndex:
                symbols_of_others = symbols_of_others.union(self._equations[i].getSymbolsUsed())
            else:
                symbols_of_this_eqn = self._equations[i].getSymbolsUsed()

        symbols_needed = symbols_of_this_eqn - symbols_of_others  # these are the symbols needed in the new equation
        Equation._logger.info("EquationSystem.replaceRandomDimensionallyConsistentEqnByIndex(): Vars,derivs and constants needed = " + str(symbols_needed))
        replacement_eqn = None
        if len(symbols_needed) == 0:
            while True:
                replacement_eqn = Equation.GenerateRandomDimensionallyConsistent(vars=self._vars, derivatives=self._derivatives,
                                        constants=self._constants, u_of_mToTermLookupDict=EquationSystem._LookupDict)
                found_dup = False
                for eqn in self._equations:
                    if replacement_eqn.equalModUnnamedConstants(eqn):
                        Equation._logger.info("Generated equation: " + str(replacement_eqn)
                                + " is, modulo unnamed constants, equal to existing equation: " + str(eqn)
                                + ". Will generate a new equation!")
                        found_dup = True
                        continue
                if not found_dup:
                    break
        else:
            while True:
                random_additional_symbol = random.choice(list(symbols_needed))
                given_var = given_derivative = given_constant = None
                if isinstance(random_additional_symbol, Constant):
                    given_constant = random_additional_symbol
                elif isinstance(random_additional_symbol, Variable):
                    given_var = random_additional_symbol
                elif isinstance(random_additional_symbol, Derivatif):
                    given_derivative = random_additional_symbol
                replacement_eqn = Equation.GenerateRandomDimensionallyConsistentEquationWithSpecifiedVarOrDerivative(vars=self._vars,
                                derivatives=self._derivatives, constants=self._constants,
                                u_of_mToTermLookupDict=EquationSystem._LookupDict,
                                given_var=given_var, given_derivative=given_derivative, given_constant=given_constant,
                                max_power=max_power,
                                max_vars_derivatives_and_constants_per_eqn=self._max_vars_derivatives_and_constants_per_eqn)
                # should check that replacement_eqn is not pone of the existing equations
                if symbols_needed.issubset(replacement_eqn.getSymbolsUsed()):
                    if replacement_eqn.equalModUnnamedConstants(self._equations[eqnIndex]):
                        Equation._logger.info("Generated equation: " + str(replacement_eqn)
                                + " is, modulo unnamed constants, equal to existing equation: " + str(self._equations[eqnIndex])
                                + ". Will generate a new equation!")
                        continue
                    else:
                        break

        Equation._logger.info("EquationSystem.replaceRandomDimensionallyConsistentEqnByIndex(): Replacement eqn: " + str(replacement_eqn))
        self._equations[eqnIndex] = replacement_eqn
        return replacement_eqn


    def replaceRamdomDimensionallyConsistentEqn(self):
        randIndex = random.randint(0, len(self._equations) - 1)
        return self.replaceRandomDimensionallyConsistentEqnByIndex(randIndex)

    def replaceRandomEqnByIndex(self, eqnIndex):  #Need dimensionally consistent variant of this method
        eqn = None
        symbols_used = set()
        for i in range(len(self._equations)):
            if i != eqnIndex:
                symbols_used = symbols_used.union(self._equations[i].getSymbolsUsed())

        while True:
            eqn = Equation.GenerateRandom(vars=self._vars, derivatives=self._derivatives, constants=self._constants,
                                          max_vars_derivs_and_constants_per_eqn=self._max_vars_derivatives_and_constants_per_eqn)
            all_symbols = symbols_used.union(eqn.getSymbolsUsed())
            if len(all_symbols) == len(self._vars):
                self._equations[eqnIndex] = eqn
                break

        return eqn

    def replaceRamdomEqn(self):    #Need dimensionally consistent variant of this method
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

        constNames = self.getConstantNames()
        exp += "Constants: \n"
        for i in range(len(constNames)):
            exp += constNames[i]
            if i < len(constNames) - 1:
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



























