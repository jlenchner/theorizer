"""The EquationSystem class, including implementations for dimensionally consistent
   and not necessarily dimensionally consistent EquationSystems. EquatioSystems are also
   referred to as "theories" in our papers.
"""

# Author: Jonathan Lenchner (lenchner@us.ibm.com)
#
# License: BSD 3-Clause

import random
from variable import *
from constant import *
from derivative import *

from sympy import *
from equation import Equation
from m2_functions import *


class EquationSystem:
    _LastVarsDerivsAndConstants = []  # The last variables, derivatives and constants used to create an EquationSystem.
    # This list is saved so that when a number of EquationSystems are created in batch
    # the varius lookup dictionaries don't have to be recreated (since their creation
    # is quite time consuming).
    _LookupDict = None  # Unit of measure (class UofM) to term lookup dictionary for the list of variables,
    # derivatives and constants specified in EquationSystem._LastVarsDerivsAndConstants.
    _sigDict = None  # A dictionary whose keys are variable-derivative-constant "signatures"
    # and whose values are their relative frequencies in randomly generated terms.
    # For example if ([x,y,z] are the associated variables (class Variable), [dxdt, dydt]
    # the associated derivatives (class Derivatif), and [c,G] the associated constants
    # (class Constant), then the term c*x*y**2*z will have signature "112..1". The leading
    # '112' means that there are three variables appearing in the term, two of which
    # are raised to the 1st power and one of which is raised to the 2nd power. The '..'
    # indicates that there are no derivatives in the term, and the final 1 indicates that
    # there is one constant in the term, and it is raised to the 1st power. Another example
    # is G**2*dxdt**2dydt, which would have signature ".12.2". Just like EquationSystem._lookupDict,
    # this dictionary is associated with the variables, derivatives and constants specified in
    # EquationSystem._LastVarsDerivsAndConstants.

    def __init__(self, vars=[], derivatives=[], constants=[], measuredVars=[], equations=[],
                 max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX):
        """
            EquationSystem constructor.
            :param vars: A list of variables (class Variable) )used in the EuqationSystem. If not supplied,
                        inferred from the provided equations.
            :param derivatives: A list of derivatives (class Derivatif) )used in the EuqationSystem. If not
                        supplied, inferred from the provided equations.
            :param constants: A list of constants (class Constants) )used in the EuqationSystem.  If not
                        supplied, inferred from the provided equations.
            :param measuredVars: A list of measured variables for the equation system, in other words, the
                        variables that are directly measured in data collection, or, equivalently, the
                        variables that need to be generated in data generation.
            :param equations: A list of equations conmprising the EuqationSystem
            :param max_vars_derivatives_and_constants_per_eqn: The maximum number of distinct variables,
                        derivatives or constants that are allowed in any equation (only utilized if there is
                        a request to replace a given equation in the EuqationSystem with another equation.
        """
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
            for var in vars_as_set:  # Turn set into an array
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
            self._measuredVars = list(vars)[:len(vars) // 2]

        self._max_vars_derivatives_and_constants_per_eqn = max_vars_derivatives_and_constants_per_eqn
        self._nonMeasuredVars = set(self._vars) - set(self._measuredVars)

    def copy(self):
        """
            Performs a deep copy of the given EquationSystem. Changes to the original EquationSystem will
            not affect the copy and vice versa.
            :return: Copy of the EuqationSystem
        """
        return EquationSystem(vars=self._vars.copy(), derivatives=self._derivatives.copy(),
                              constants=self._constants.copy(),
                              measuredVars=self._measuredVars.copy(), equations=self._equations.copy(),
                              max_vars_derivatives_and_constants_per_eqn=self._max_vars_derivatives_and_constants_per_eqn)


    def isConsistent(self):
        """
            Checks the consistency of this EquationSystem.
            :return: True if consistent, False otherwise
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
            Checks the consistency of the current system of equations after replacing a random equation
            at a given index a given num number of times.
            :param index: the index of the equation that was replaced
            :param num: the number of times to replace the equation
            :return: the number of consistent systems found
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
            Projects the equations in the current EquationSystem onto the measured variables
            and writes results to the file projection_output.txt.
            :return: A string version of the file
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
            Projects the equations in the current EquationSystem onto the measured variables
            and writes results to the file projection_output.txt.
            :return: A string version of the file
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
            :param vars: List of variables
            :param derivatives: List of derivatives (class Derivatif)
            :param measured_vars: List of measured variables
            :param num: Number of random systems to project
            :return: a list of results from the projections
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
            Checks the consistency of a given number of randomly generated EquationSystems
            with a given set of variables.
            :param vars: Variables to use in generation of random systems
            :param num: the number of random systems to check
            :return: the number of consistent systems found
        """
        numConsistent = 0
        for i in range(num):
            eqnSystem = EquationSystem.GenerateRandom(vars, 4, 6)
            if eqnSystem.isConsistent():
                numConsistent += 1
        return numConsistent

    def getVarNames(self):
        """
            Returns a list of variable names for the variables used by the EquationSystem
            :return: List of variable names
        """
        var_names = []
        for var in self._vars:
            var_names.append(str(var))

        return var_names

    def getDerivNames(self):
        """
            Returns a list of derivative names for the derivatives used by the EquationSystem
            :return: List of derivative names
        """
        deriv_names = []
        for derirvative in self._derivatives:
            deriv_names.append(str(derirvative))

        return deriv_names

    def getConstantNames(self):
        """
            Returns a list of constant names for the constants used by the EquationSystem
            :return: List of constant names
        """
        constant_names = []
        for constant in self._constants:
            constant_names.append(str(constant))

        return constant_names

    def getMeasuredVars(self):
        """
            Returns a list of the names of the measured variables used by the EquationSystem
            :return: List of measured variable names
        """
        measured_var_names = []
        for var in self._measuredVars:
            measured_var_names.append(str(var))

        return measured_var_names

    def getNonMeasuredVars(self):
        """
            Returns a list of the names of the non-measured variables used by the EquationSystem
            :return: List of non-measured variable names
        """
        non_measured_var_names = []
        for var in self._nonMeasuredVars:
            non_measured_var_names.append(str(var))

        return non_measured_var_names

    def getEquations(self):
        """
            Returns a list of stringified versions of the equations in the EquationSystem
            :return: List of stringified equations
        """
        equation_strings = []
        for eqn in self._equations:
            equation_strings.append(str(eqn))

        return equation_strings

    @staticmethod
    def _get_var_by_name(seq, name):
        for v in seq:
            if str(v) == name:
                return v
        return None

    @staticmethod
    def _is_linear_in_var(eqn, var):
        """Return True if eqn’s polynomial is linear (deg ≤ 1) in `var`."""
        if var is None:
            return True
        try:
            p = Poly(eqn._poly, var)
            return p.degree(var) <= 1
        except Exception:
            # If Poly fails (rare), default to accepting to avoid false negatives
            return True
        
    @classmethod
    def GenerateRandom(cls, vars, derivatives, constants, measuredVars, numEqns,
                       max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX):
        """
            Creates a randomly generated, not-necessarily dimensionally consistent EquationSystem
            with the given variables, derivatives and constants abiding by the restriction on the
            maximum number of distinct variables, derivatives and constants in any one equation
            :param vars:List of variables (class Variable) to use
            :param derivatives: List of derivatives (class Derivatif) to use
            :param constants: List of constants (class Constant) to use
            :param measuredVars: List of measured variables (does not affect equation generation)
            :param numEqns: A positive integer indicating the number of equations to generate
            :param max_vars_derivatives_and_constants_per_eqn: The maximum number of distinct
                        variables, derivatives and constants in any one equation
            :return: The randomly generated EquationSystem
        """
        eqns = []

        while True:  # A crude way to make sure we use all the vars, derivatives and constatns
            symbols_used = set()
            for i in range(numEqns):
                eqn = Equation.GenerateRandom(vars=vars, derivatives=derivatives, constants=constants,
                                              max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
                symbols_used = symbols_used.union(eqn.getSymbolsUsed())
                eqns.append(eqn)

            if len(symbols_used) < len(vars) + len(derivatives) + len(constants):
                eqns = []
            elif len(eqns) != len(set(eqns)):  # means there are duplicate equations
                eqns = []
            else:
                break

        return EquationSystem(vars=vars, derivatives=derivatives, constants=constants, measuredVars=measuredVars,
                              equations=eqns,
                              max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)

    @classmethod
    def GenerateRandomDimensionallyConsistent(cls, vars, derivatives, constants, measuredVars, numEqns,
                                            max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX):
        """
            Creates a randomly generated, dimensionally consistent EquationSystem with the given
            variables, derivatives and constants. This version does NOT perform any global 'sanity'
            checks at the end, but it DOES enforce: degree(theta) <= 1 if a variable named 'theta'
            exists among `vars`.
        """
        while True:  # keep outer loop semantics identical to your original
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
                EquationSystem._LookupDict = Equation.GetUofMToPrimitiveTermLookupTable(
                    vars, derivatives, constants, max_power=max_power
                )
                EquationSystem._sigDict = Equation.GenerateVDCSigDistributionDict(
                    vars, derivatives, constants, max_power=max_power
                )
            else:
                Equation._logger.info("Using existing lookup & sig dictionaries.")

            # ---- helpers (local, self-contained) ----
            def _get_var_by_name(seq, name):
                for v in seq:
                    if str(v) == name:
                        return v
                return None

            theta_var = _get_var_by_name(vars, "theta")

            def _is_linear_in_theta(eqn):
                if theta_var is None:
                    return True
                try:
                    p = Poly(eqn._poly, theta_var)
                    return p.degree(theta_var) <= 1
                except Exception:
                    # If Poly construction fails (rare), be permissive
                    return True

            def _gen_dc_eqn():
                # resample until linear in theta
                while True:
                    e = Equation.GenerateRandomDimensionallyConsistent(
                        vars=vars, derivatives=derivatives, constants=constants,
                        u_of_mToTermLookupDict=EquationSystem._LookupDict,
                        sigDict=EquationSystem._sigDict,
                        max_power=max_power,
                        max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn
                    )
                    if _is_linear_in_theta(e):
                        return e

            def _gen_dc_eqn_given(given_var=None, given_derivative=None, given_constant=None):
                # resample until linear in theta (and not None)
                while True:
                    e = Equation.GenerateRandomDimensionallyConsistentEquationWithSpecifiedVarDerivOrConstant(
                        vars=vars, derivatives=derivatives, constants=constants,
                        u_of_mToTermLookupDict=EquationSystem._LookupDict,
                        sigDict=EquationSystem._sigDict,
                        given_var=given_var, given_derivative=given_derivative, given_constant=given_constant,
                        max_power=max_power,
                        max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn
                    )
                    if e is not None and _is_linear_in_theta(e):
                        return e

            # ---- first equation ----
            eqn = _gen_dc_eqn()
            symbols_used = eqn.getSymbolsUsed()
            eqns.append(eqn)
            Equation._logger.info("Eqn: " + str(eqn))

            # ---- cover remaining symbols at least once ----
            varsDerivsAndConstantsRemainingToBeUsed = varsDerivsAndConstantsRemainingToBeUsed - symbols_used
            while len(varsDerivsAndConstantsRemainingToBeUsed) > 0 and len(eqns) < numEqns:
                nextVarDerivativeOrConstant = next(iter(varsDerivsAndConstantsRemainingToBeUsed))
                next_var = next_derivative = next_constant = None
                if isinstance(nextVarDerivativeOrConstant, Constant):
                    next_constant = nextVarDerivativeOrConstant
                elif isinstance(nextVarDerivativeOrConstant, Variable):
                    next_var = nextVarDerivativeOrConstant
                elif isinstance(nextVarDerivativeOrConstant, Derivatif):
                    next_derivative = nextVarDerivativeOrConstant

                next_eqn = _gen_dc_eqn_given(
                    given_var=next_var, given_derivative=next_derivative, given_constant=next_constant
                )
                if next_eqn is not None:
                    eqns.append(next_eqn)
                    Equation._logger.info("Eqn: " + str(next_eqn))
                    symbols_used = next_eqn.getSymbolsUsed()
                    varsDerivsAndConstantsRemainingToBeUsed = varsDerivsAndConstantsRemainingToBeUsed - symbols_used
                else:
                    # fall back: drop the requirement if generator couldn't satisfy it
                    varsDerivsAndConstantsRemainingToBeUsed.remove(nextVarDerivativeOrConstant)
                    Equation._logger.info("Dropping var-deriv-or-constant: " + str(nextVarDerivativeOrConstant))

            # ---- fill up to numEqns with fresh UoMs not in use (same behavior as your code) ----
            while len(eqns) < numEqns:
                u_of_ms_in_use = set()
                for e in eqns:
                    firstTerm = e.getTerms()[0]
                    u_of_m_for_term = Equation.GetUofMForTerm(firstTerm)
                    u_of_ms_in_use.add(u_of_m_for_term)
                while True:
                    e = _gen_dc_eqn()
                    firstTerm = e.getTerms()[0]
                    u_of_m_for_term = Equation.GetUofMForTerm(firstTerm)
                    if u_of_m_for_term not in u_of_ms_in_use:
                        eqns.append(e)
                        Equation._logger.info("Eqn: " + str(e))
                        break

            # ---- No global sanity check here; return immediately ----
            return EquationSystem(
                vars=vars, derivatives=derivatives, constants=constants,
                measuredVars=measuredVars, equations=eqns,
                max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn
            )

    
    @classmethod
    def DetermineMaxPower(cls, vars, derivatives, constants):
        """
            A somewhat ad hoc method of determining the maximum power to raise variables,
            derivatives and (named) constants to within equations to keep the compute time
            manageable. Important to use this consistently across all methods.
            :param vars: List of variables (class Variable) to base determination upon
            :param derivatives: List of derivatives (class Derivatif) to base determination upon
            :param constants: List of constants to base determination upon
            :return: An integer giving the max power
        """
        varsDerivsAndConstants = []
        varsDerivsAndConstants.extend(vars)
        varsDerivsAndConstants.extend(derivatives)
        varsDerivsAndConstants.extend(constants)

        if len(varsDerivsAndConstants) > 12: #This number can be changed to anything up to about 20
                                             # (or even larger if running on powerful HW)
            return 2
        else:
            return 3

    def replaceRandomDimensionallyConsistentEqnByIndex(self, eqnIndex):
        """
            Replace an equation at the given index in this EquationSystem with a random dimensionally
            consistent equation (using the existing variables, derivatives and constants), with the
            additional constraint that the replacement is linear in the variable named 'theta'
            (degree <= 1) if that variable exists in this system.
            :param eqnIndex: index of the equation in ._equations to replace
            :return: the replacement equation (though the equation is replaced in place upon return)
        """
        # --- Existing safety checks ---
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

        # Assemble symbol universe
        varsDerivsAndConstants = []
        varsDerivsAndConstants.extend(self._vars)
        varsDerivsAndConstants.extend(self._derivatives)
        varsDerivsAndConstants.extend(self._constants)

        max_power = EquationSystem.DetermineMaxPower(self._vars, self._derivatives, self._constants)

        # Ensure lookup/sig dicts are up to date
        if varsDerivsAndConstants != EquationSystem._LastVarsDerivsAndConstants:
            EquationSystem._LastVarsDerivsAndConstants = varsDerivsAndConstants
            EquationSystem._LookupDict = Equation.GetUofMToPrimitiveTermLookupTable(
                vars=self._vars,
                derivatives=self._derivatives,
                constants=self._constants,
                max_power=max_power
            )
            EquationSystem._sigDict = Equation.GenerateVDCSigDistributionDict(
                vars=self._vars,
                derivatives=self._derivatives,
                constants=self._constants,
                max_power=max_power
            )

        # Figure out which symbols must appear in the replacement (to preserve coverage)
        symbols_of_others = set()
        symbols_of_this_eqn = set()
        for i in range(len(self._equations)):
            if i != eqnIndex:
                symbols_of_others = symbols_of_others.union(self._equations[i].getSymbolsUsed())
            else:
                symbols_of_this_eqn = self._equations[i].getSymbolsUsed()

        symbols_needed = symbols_of_this_eqn - symbols_of_others  # required in the new equation
        Equation._logger.info(
            "EquationSystem.replaceRandomDimensionallyConsistentEqnByIndex(): Vars,derivs and constants needed = " + str(
                symbols_needed))

        # Find the actual Variable object named "theta" (dimensionless var you added)
        theta_var = next((v for v in self._vars if str(v) == "theta"), None)

        def _is_linear_in_theta(eqn):
            """Return True if eqn has degree <= 1 in theta (or if theta not present)."""
            if theta_var is None:
                return True
            try:
                # eqn._poly is a SymPy polynomial in all symbols; project to theta and check degree
                p = Poly(eqn._poly, theta_var)
                return p.degree(theta_var) <= 1
            except Exception:
                # Be permissive on unexpected forms (rare); treat as OK
                return True

        replacement_eqn = None

        if len(symbols_needed) == 0:
            # No specific symbol required: sample until (1) linear in theta, (2) passes sanity checks, (3) not dup
            while True:
                replacement_eqn = Equation.GenerateRandomDimensionallyConsistent(
                    vars=self._vars,
                    derivatives=self._derivatives,
                    constants=self._constants,
                    u_of_mToTermLookupDict=EquationSystem._LookupDict,
                    sigDict=EquationSystem._sigDict
                )
                if not _is_linear_in_theta(replacement_eqn):
                    # Enforce theta-degree cap
                    continue

                if not self.sanityCheckReplacementEquation(eqnIndex, replacement_eqn):
                    Equation._logger.info(
                        "Equation " + str(replacement_eqn) + " did not pass EquationSystem sanityCheck!")
                    Equation._logger.info("Other equations: ")
                    for i in range(len(self._equations)):
                        if i != eqnIndex:
                            Equation._logger.info(str(self._equations[i]))
                    continue

                found_dup = False
                for eqn in self._equations:
                    if replacement_eqn.equalModUnnamedConstants(eqn):
                        Equation._logger.info("Generated equation: " + str(replacement_eqn)
                                            + " is, modulo unnamed constants, equal to existing equation: " + str(eqn)
                                            + ". Will generate a new equation!")
                        found_dup = True
                        break
                if not found_dup:
                    break

        else:
            # Must include some missing symbols from the old equation
            while True:
                random_additional_symbol = random.choice(list(symbols_needed))
                given_var = given_derivative = given_constant = None
                if isinstance(random_additional_symbol, Constant):
                    given_constant = random_additional_symbol
                elif isinstance(random_additional_symbol, Variable):
                    given_var = random_additional_symbol
                elif isinstance(random_additional_symbol, Derivatif):
                    given_derivative = random_additional_symbol

                replacement_eqn = Equation.GenerateRandomDimensionallyConsistentEquationWithSpecifiedVarDerivOrConstant(
                    vars=self._vars,
                    derivatives=self._derivatives,
                    constants=self._constants,
                    u_of_mToTermLookupDict=EquationSystem._LookupDict,
                    sigDict=EquationSystem._sigDict,
                    given_var=given_var,
                    given_derivative=given_derivative,
                    given_constant=given_constant,
                    max_power=max_power,
                    max_vars_derivatives_and_constants_per_eqn=self._max_vars_derivatives_and_constants_per_eqn
                )

                if replacement_eqn is None:
                    continue

                if not _is_linear_in_theta(replacement_eqn):
                    # Enforce theta-degree cap
                    continue

                # Check we actually covered the required symbols
                if not symbols_needed.issubset(replacement_eqn.getSymbolsUsed()):
                    continue

                if replacement_eqn.equalModUnnamedConstants(self._equations[eqnIndex]):
                    Equation._logger.info("Generated equation: " + str(replacement_eqn)
                                        + " is, modulo unnamed constants, equal to existing equation: "
                                        + str(self._equations[eqnIndex]) + ". Will generate a new equation!")
                    continue

                if not self.sanityCheckReplacementEquation(eqnIndex, replacement_eqn):
                    Equation._logger.info(
                        "Equation " + str(replacement_eqn) + " did not pass EquationSystem sanityCheck!")
                    Equation._logger.info("Other equations: ")
                    for i in range(len(self._equations)):
                        if i != eqnIndex:
                            Equation._logger.info(str(self._equations[i]))
                    continue

                break

        Equation._logger.info(
            "EquationSystem.replaceRandomDimensionallyConsistentEqnByIndex(): Replacement eqn: " + str(replacement_eqn))
        self._equations[eqnIndex] = replacement_eqn
        return replacement_eqn  


    def _is_trig_identity_eqn(self, idx_or_eqn):
        """Return True if the equation (by index or object) is sin^2+cos^2-1 = 0 (up to sign)."""
        eqn = self._equations[idx_or_eqn] if isinstance(idx_or_eqn, int) else idx_or_eqn

        sin_v = EquationSystem._get_var_by_name(self._vars, "sinTheta")
        cos_v = EquationSystem._get_var_by_name(self._vars, "cosTheta")
        if sin_v is None or cos_v is None:
            return False

        try:
            expr = eqn._poly.as_expr()
            ident = sin_v**2 + cos_v**2 - Integer(1)
            # accept either +ident==0 or -ident==0
            return simplify(expr - ident) == 0 or simplify(expr + ident) == 0
        except Exception:
            return False


    def replaceRandomDimensionallyConsistentEqn(self):
        """
        Replace a randomly selected equation with a random DC equation,
        but never replace sin_theta^2 + cos_theta^2 = 1 (just resample).
        """
        eligible = [i for i in range(len(self._equations)) if not self._is_trig_identity_eqn(i)]
        if not eligible:
            return -1, None  # nothing eligible to replace
        randIndex = random.choice(eligible)
        return randIndex, self.replaceRandomDimensionallyConsistentEqnByIndex(randIndex)


    def replaceRandomEqnByIndex(self, eqnIndex):
        """
            Replace an equation at the given index in this EquationSystem with a random not necessarily
            dimensionally consistent equation (using the exiting variables, derivatives and constants)
            :param eqnIndex: index of the equation in ._equations to replace
            :return: the replacement equation (though the equation is replaced in place upon return)
        """
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

    def replaceRamdomEqn(self):
        """
            Replace an equation at a randomly selected index in this EquationSystem with a random
            not necessarily dimensionally consistent equation (using the exiting variables,
            derivatives and constants)
            :return: A 2-tuple consisting of the index of the replaced equation and the replacement
                    equation itself (though the equation is replaced in place upon return)
        """
        randIndex = random.randint(0, len(self._equations) - 1)
        return randIndex, self.replaceRandomEqnByIndex(randIndex)

    def sanityCheck(self):
        """
            Performs various sanity checks this EquationSystem (for details see the method
            SanityCheckEquationList())
            :return: True of EquationSystem passes the sanity check(s), False otherwise
        """
        return EquationSystem.SanityCheckEquationList(self._equations)

    @classmethod
    def SanityCheckEquationList(cls, eqns):
        """
            Sanity checks a list of equations. Currently, does a single sanity check: verifying that
            all equations don't have the same unit of measure. There is an implicit assumption that
            the equations are all dimensionally consistent (though the method does not fail if not --
            the sanity check is just not meaningful in this case, since it is done based on the UofM
            of the first terms of all equations).
            :param eqns: List of equations to sanity check
            :return: True of the equations pass the sanity check, False otherwise.
        """
        uOofMsSoFar = []
        for i in range(len(eqns)):
            firstTermOfEqn = eqns[i].getTerms()[0]
            u_of_m = Equation.GetUofMForTerm(firstTermOfEqn)
            if u_of_m in uOofMsSoFar:
                j = uOofMsSoFar.index(u_of_m)
                i_vars = set(eqns[i]._variables)
                i_derivs = set(eqns[i]._derivatives)
                j_vars = set(eqns[j]._variables)
                j_derivs = set(eqns[j]._derivatives)
                if not i_vars.isdisjoint(j_vars) or not i_derivs.isdisjoint(j_derivs):
                    return False
                else:
                    uOofMsSoFar.append(u_of_m)
            else:
                uOofMsSoFar.append(u_of_m)

        return True

    def sanityCheckReplacementEquation(self, replacementIndex, eqn):
        """
            Performs the same sanity check as SanityCheckEquationList() but for this EquationSystem,
            with a specified replacement equation at a specified index
            :param replacementIndex: Index to slot in replacement equation for sanity check
            :param eqn: Replacement equation
            :return: True if EquationSystem with designated replacement passes the sanity check, False
                    otherwise.
        """
        u_of_m_for_replacement = Equation.GetUofMForTerm(eqn.getTerms()[0])
        for i in range(len(self._equations)):
            if i == replacementIndex:
                continue
            firstTermOfEqn = self._equations[i].getTerms()[0]
            u_of_m = Equation.GetUofMForTerm(firstTermOfEqn)
            if u_of_m_for_replacement == u_of_m:
                i_vars = set(self._equations[i]._variables)
                i_derivs = set(self._equations[i]._derivatives)
                j_vars = set(eqn._variables)
                j_derivs = set(eqn._derivatives)
                if not i_vars.isdisjoint(j_vars) or not i_derivs.isdisjoint(j_derivs):
                    return False

        return True
    
    def addEquation(self, eqn, keep_count=False):
        """
        Add a specified Equation to this system, ensuring dimensional consistency.
        If keep_count=True and the system already has equations, replace one
        (preferably same UoM) to keep the total count unchanged.

        Returns the index at which the equation ended up.
        """
        # Must be dimensionally consistent by itself
        if eqn.getUofM() is None:
            raise ValueError("Provided equation is not dimensionally consistent (UoM mismatch among its terms).")

        if not keep_count or len(self._equations) == 0:
            self._equations.append(eqn)
            return len(self._equations) - 1

        # keep_count=True: try to replace an equation with same UoM while preserving coverage
        target_uom = eqn.getUofM()
        desired_coverage = set(self._vars) | set(self._derivatives) | set(self._constants)

        def coverage_with(eqns_replacement_index=None):
            cov = set()
            for i, e in enumerate(self._equations):
                if i == eqns_replacement_index:
                    continue
                cov |= set(e.getSymbolsUsed())
            cov |= set(eqn.getSymbolsUsed())
            return cov

        # Prefer replacing an equation with the same UoM that does not reduce symbol coverage
        candidate_index = None
        for i, e in enumerate(self._equations):
            try:
                if e.getUofM() == target_uom and desired_coverage.issubset(coverage_with(eqns_replacement_index=i)):
                    candidate_index = i
                    break
            except Exception:
                pass

        if candidate_index is None:
            # Fallback: replace the last equation
            candidate_index = len(self._equations) - 1

        self._equations[candidate_index] = eqn
        return candidate_index


    def __str__(self):
        varNames = self.getVarNames()
        exp = "Variables: \n"
        for i in range(len(varNames)):
            exp += varNames[i]
            exp += ' (' + str(self._vars[i]._u_of_m) + ')'
            if i < len(varNames) - 1:
                exp += ', '
        exp += '\n'
        exp += '\n'

        derivNames = self.getDerivNames()
        exp += "Derivatives: \n"
        for i in range(len(derivNames)):
            exp += derivNames[i]
            exp += ' (' + str(self._derivatives[i]._u_of_m) + ')'
            if i < len(derivNames) - 1:
                exp += ' '
        exp += '\n'
        exp += '\n'

        constNames = self.getConstantNames()
        exp += "Constants: \n"
        for i in range(len(constNames)):
            exp += constNames[i]
            exp += ' = ' + str(self._constants[i]._value) + ' (' + str(self._constants[i]._u_of_m) + ')'
            if i < len(constNames) - 1:
                exp += ', '
        exp += '\n'
        exp += '\n'

        #Measured and non-measured variables are no longer reported
        # measuredVarNames = self.getMeasuredVars()
        # exp += "Measured Variables: \n"
        # for i in range(len(measuredVarNames)):
        #    exp += measuredVarNames[i]
        #    if i < len(measuredVarNames) - 1:
        #        exp += ', '
        #
        # exp += '\n'
        # exp += '\n'
        # exp += "Non Measured Variables: \n"
        # nonMeasuredVarNames = self.getNonMeasuredVars()
        # for i in range(len(nonMeasuredVarNames)):
        #    exp += nonMeasuredVarNames[i]
        #    if i < len(nonMeasuredVarNames) - 1:
        #        exp += ', '

        exp += '\n'
        exp += '\n'
        exp += "Equations: \n"
        for i in range(len(self._equations)):
            exp += str(self._equations[i])
            if i < len(self._equations) - 1:
                exp += '\n'

        return exp













