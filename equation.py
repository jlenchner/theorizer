from sympy import *
from variable import *
from derivative import *
from constant import *
import random
import math
import logging
import itertools


class Equation:
    PROB_OF_NUM_TERMS = [0, 42.0/70.0, 20.0/70.0, 6.0/70.0, 2/70.0]  #probs of a randomly generated poly having 1 term, 2 terms, etc.
    PROB_FACTORS_PER_TERM = [0.19, 0.47, 0.25, 0.09] #probs that a term in a randomly generated poly will have 1, 2, 3, 4 factors
    SAME_VARIABLE_FACTOR_BIAS = 2.0  #bias for repeating the same variable in a term (in other words, favoring x^2 over xy)
    SAME_DERIVATIVE_FACTOR_BIAS = 0.5 #bias for repeating the same derivative in a term (in other words, favoring dxdt^2 over dxdt*dydt)

    PROB_TERM_HAS_NON_UNITAL_INTEGER_CONSTANT = 0.2
    PROB_OF_SMALL_INTEGER_CONSTANTS = [0.0, 0.0, 15.0/34.0, 10.0/34.0, 9.0/34.0]  # prob that a small non-unital integer constant is [0,1,2,3,4,5]

    HOMOGENEOUS_EQN_PROB = 0.05

    NO_MAX = 999

    _logger = logging.getLogger()


    def __init__(self, exp):
        self._poly = Poly(exp)  #The equation is treated generically like a polynomial of symbols
        self._variables, self._derivatives, self._constants = Equation.InferVarsDerivativesAndConstantsFromExpression(exp)


    def isDimensionallyConsistent(self):
        for var in self._variables:
            if var._u_of_m is None:
                Equation._logger.warning("The variable " + str(var) + " has no associated units of measure!")
                return False
        for derivative in self._derivatives:
            if derivative._u_of_m is None:
                Equation._logger.warning("The derivative " + str(derivative) + " has no associated units of measure!")
                return False
        for constant in self._constants:
            if constant._u_of_m is None:
                Equation._logger.warning("The constant " + str(constant) + " has no associated units of measure!")
                return False

        UofMs = []
        terms = self.getTerms()
        for term in terms:  #may need to be a term class that contains the symbols along with the units of measure....
            UofM = Equation.GetUofMForTerm(term)
            UofMs.append(UofM)

        #Now check whether they are all equal
        firstUofM = None
        for UofM in UofMs:
            if firstUofM is None:
                firstUofM = UofM
            else:
                if firstUofM != UofM:
                    return False

        return True

    @classmethod
    def GetUofMForTerm(cls, term):
        return UofM(Equation.InnerGetUofMForTerm(term))

    @classmethod
    def InnerGetUofMForTerm(cls, term):
        if isinstance(term, Variable) or isinstance(term, Derivatif):
            return term._u_of_m._units
        elif isinstance(term, Mul):
            args = []
            for arg in term.args:
                raw_u_of_m = Equation.InnerGetUofMForTerm(arg)
                if raw_u_of_m is not None:
                    args.append(raw_u_of_m)
            return Mul(*args)
        elif isinstance(term, Pow):
            return Pow(term.args[0]._u_of_m._units, term.args[1])
        else:
            return None

    def getUofM(self):
        terms = self.getTerms()
        firstUofM = Equation.GetUofMForTerm(terms[0])
        for i in range(1, len(terms)):
            UofMforTerm = Equation.GetUofMForTerm(terms[i])
            if firstUofM != UofMforTerm:
                return None

        return firstUofM



    def add(self, exp, vars = [], derivatives = [], constants = []):
        self._poly = self._poly.add(Poly(exp))
        if len(vars) == 0 and len(derivatives) == 0:
            vars, derivatives, constants = Equation.InferVarsDerivativesAndConstantsFromExpression(exp)
        for var in vars:
            if var not in self._variables:
                self._variables.append(var)
        for deriv in derivatives:
            if deriv not in self._derivatives:
                self._derivatives.append(deriv)
        for constant in constants:
            if constant not in self._constants:
                self._constants.append(constant)



    @classmethod
    def InferVarsDerivativesAndConstantsFromExpression(cls, exp):
        vars = []
        derivtives = []
        constants = []
        for sym in exp.free_symbols:
            if isinstance(sym, Constant): #must go first since Constants are (constant) Variables!
                constants.append(sym)
            elif isinstance(sym, Variable):
                vars.append(sym)
            elif isinstance(sym, Derivatif):
                derivtives.append(sym)

        return vars, derivtives, constants



    @classmethod
    def GetRandomNumberOfTerms(cls):
        r = random.random()
        num_terms = len(Equation.PROB_OF_NUM_TERMS)
        prob_so_far = 0.0
        for i in range(len(Equation.PROB_OF_NUM_TERMS)):
            prob_so_far += Equation.PROB_OF_NUM_TERMS[i]
            if r < prob_so_far:
                num_terms = i + 1
                break

        return num_terms

    @classmethod
    def GenerateRandom(cls, vars, derivatives=[], constants=[], max_vars_derivatives_and_constants_per_eqn=NO_MAX):
        if max_vars_derivatives_and_constants_per_eqn == Equation.NO_MAX:
            max_vars_derivatives_and_constants_per_eqn = len(vars) + len(derivatives) + len(constants)

        num_terms = Equation.GetRandomNumberOfTerms()

        terms = None
        while True:
            firstTerm = Equation.GenerateRandomTerm(vars, derivatives, constants)
            terms = [firstTerm]
            vars_derivs_and_constants_in_use = firstTerm.free_symbols
            for i in range(1, num_terms):
                while True:
                    term = Equation.GenerateRandomTerm(vars, derivatives, constants)
                    vars_derivs_and_constants_in_term = term.free_symbols
                    new_vars_derivs_and_constants = vars_derivs_and_constants_in_term - vars_derivs_and_constants_in_use
                    if len(new_vars_derivs_and_constants) + len(vars_derivs_and_constants_in_use) <= max_vars_derivatives_and_constants_per_eqn \
                            and not Equation.TermAmongExistingTerms(terms, term):
                        vars_derivs_and_constants_in_use = vars_derivs_and_constants_in_use.union(new_vars_derivs_and_constants)
                        terms.append(term)
                        break
            if Equation.GetCommonFactors(terms) == set():
                if Equation.SanityCheckFromTerms(terms):
                    break
                else:
                    exp = terms[0]
                    for i in range(1, len(terms)):
                        exp = exp + terms[i]
                    Equation._logger.info("Equation: " + str(exp) + " has failed sanityCheck!")


        Equation.AssignRandomSignsToTerms(terms)

        exp = terms[0]
        for i in range(1, len(terms)):
            exp = exp + terms[i]

        eqn = Equation(exp)
        eqn.divideByCommonUnnamedConstants()

        return eqn

    @classmethod
    def GetAllDimensionallyConsistentTerms(cls, term_to_match, vars=[], derivatives=[], constants=[], max_power=3,
                                           max_varsDerivsAndConstants=NO_MAX):
        if max_varsDerivsAndConstants == Equation.NO_MAX:
            max_varsDerivsAndConstants = len(vars) + len(derivatives) + len(constants)

        allConsistentTerms = []
        u_of_m_to_match = Equation.GetUofMForTerm(term_to_match)
        for i in range(1, pow(max_power + 1, len(vars) + len(derivatives) + len(constants))):
            baseNum = Equation.ToBase(i, max_power + 1)
            s_baseNum = str(baseNum)
            if len(s_baseNum) < len(vars) + len(derivatives) + len(constants):
                s_baseNum = '0' * (len(vars) + len(derivatives) + len(constants) - len(s_baseNum)) + s_baseNum
            numZeroes = s_baseNum.count('0')
            numVarsDerivsAndConstants = len(vars) + len(derivatives) + len(constants) - numZeroes
            if numVarsDerivsAndConstants <= max_varsDerivsAndConstants:
                term = Equation.GenerateTermFromBaseNum(baseNum, vars, derivatives, constants)
                u_of_m = Equation.GetUofMForTerm(term)
                if u_of_m_to_match == u_of_m and not Equation.TermsEqualModUnnamedConstants(term, term_to_match):
                    allConsistentTerms.append(term)

        return allConsistentTerms

    @classmethod
    def ToBase(cls, num, base):
        res = 0
        if num == 0:
            return 0
        max_pow = ceiling(math.log(num, base))
        for i in range(max_pow+1):
            res += pow(10,i)*(num % base)
            num -= (num % base)
            num = int(num/base)

        return int(res)

    @classmethod
    def GenerateTermFromTupleAndBaseNum(cls, tuple, s_baseNum):
        term = None

        for i in range(len(s_baseNum)):
            val = int(s_baseNum[i])
            if val > 0:
                varDerivOrConstant = tuple[i]
                if term is not None:
                    if val == 1:
                        term = term*varDerivOrConstant
                    else:
                        term = term * (varDerivOrConstant**val)
                else:
                    if val == 1:
                        term = varDerivOrConstant
                    else:
                        term = varDerivOrConstant**val

        return term

    @classmethod
    def GenerateTermFromBaseNum(cls, baseNum, vars=[], derivatives=[], constants=[]):
        term = None

        varsDerivsAndConstants = []
        varsDerivsAndConstants.extend(vars)
        varsDerivsAndConstants.extend(derivatives)
        varsDerivsAndConstants.extend(constants)

        s_baseNum = str(baseNum)
        if  len(s_baseNum) < len(varsDerivsAndConstants):
            s_baseNum = '0' * (len(varsDerivsAndConstants) - len(s_baseNum)) + s_baseNum

        for i in range(len(s_baseNum)):
            val = int(s_baseNum[i])
            if val > 0:
                varDerivOrConstant = varsDerivsAndConstants[i]
                if term is not None:
                    if val == 1:
                        term = term*varDerivOrConstant
                    else:
                        term = term * (varDerivOrConstant**val)
                else:
                    if val == 1:
                        term = varDerivOrConstant
                    else:
                        term = varDerivOrConstant**val

        return term

    @classmethod
    def GenerateFixedNumberofDimensionallyConsistentTermsFromList(cls, num_terms, candidateTerms, existing_terms=[],
                                                                  max_vars_derivatives_and_constants_per_eqn=NO_MAX):
       #Note that num_terms is used for guidance. May not in principle be possible to obtain this many terms with no common factor with existing_terms
       MAX_TRIES = 20
       tries = 0
       while True:
            tries += 1
            if tries > MAX_TRIES:
                if num_terms < len(candidateTerms):
                    num_terms += 1
                    tries = 0
                else:
                    return []
            candidateTermsCopy = candidateTerms.copy()
            terms = []
            vars_derivs_and_constants_in_use = set()
            for term in existing_terms:
                vars_derivs_and_constants_in_use.update(term.free_symbols)

            for i in range(num_terms):
                while True:
                    rand_int = random.randint(0, len(candidateTermsCopy) - 1)
                    term = candidateTermsCopy[rand_int]
                    vars_derivs_and_constants_in_term = term.free_symbols
                    new_vars_derivs_and_constants = vars_derivs_and_constants_in_term - vars_derivs_and_constants_in_use
                    if len(new_vars_derivs_and_constants) + len(vars_derivs_and_constants_in_use) <= max_vars_derivatives_and_constants_per_eqn:
                        vars_derivs_and_constants_in_use = vars_derivs_and_constants_in_use.union(new_vars_derivs_and_constants)
                        c = Equation.GenerateRandomUnnamedConstant()
                        term = c * term
                        terms.append(term)
                        candidateTermsCopy.remove(Equation.GetUnnamedConstantStrippedTerm(term))
                        break

            all_terms = []
            all_terms.extend(existing_terms)
            all_terms.extend(terms)
            if Equation.GetCommonFactors(all_terms) == set():
                return terms


    @classmethod
    def GenerateRandomDimensionallyConsistent(cls, vars, derivatives, constants, u_of_mToTermLookupDict=None,
                                              max_power=3,
                                              max_vars_derivatives_and_constants_per_eqn=NO_MAX):  # Does not yet deal with picking coefficients
        if max_vars_derivatives_and_constants_per_eqn == Equation.NO_MAX:
            max_vars_derivatives_and_constants_per_eqn = len(vars) + len(derivatives) + len(constants)

        if u_of_mToTermLookupDict is None:
            u_of_mToTermLookupDict = Equation.GetUofMToPrimitiveTermLookupTable(vars=vars, derivatives=derivatives,
                                                constants=constants, max_power=max_power)

        terms = None
        while True:
            num_terms = Equation.GetRandomNumberOfTerms()

            firstTerm = Equation.GenerateRandomTerm(vars, derivatives, constants, max_power=max_power)
            terms = [firstTerm]
            u_of_m_to_match = Equation.GetUofMForTerm(firstTerm)
            candidateTerms = set()
            try:
                candidateTerms = u_of_mToTermLookupDict.get(u_of_m_to_match._units)[0].copy()
            except:  #should not happen; just for safety
                Equation._logger.error("There are no candidate terms for first term:" + str(firstTerm))
                continue

            if len(candidateTerms) < num_terms or Equation.GetCommonFactors(candidateTerms) != set():
                continue
            constantStrippedFirstTerm = Equation.GetUnnamedConstantStrippedTerm(firstTerm)
            try:
                candidateTerms.remove(constantStrippedFirstTerm)
            except: #I believe I have fixed the cause of this exception but just in case....
                Equation._logger.error("Could not remove constantStrippedFirstTerm=" + str(constantStrippedFirstTerm) + \
                                       ". First term was: " + str(firstTerm))
                continue

            additionalTerms = Equation.GenerateFixedNumberofDimensionallyConsistentTermsFromList(num_terms=num_terms - 1,
                                                                                                 candidateTerms=candidateTerms,
                                                                                                 existing_terms=terms)
            if additionalTerms == []: #could not get needed additional terms
                continue

            terms.extend(additionalTerms)
            if Equation.GetCommonFactors(terms) == set(): #this should always happen because it is
                # handled in Equation.GenerateFixedNumberofDimensionallyConsistentTermsFromList()
                if Equation.SanityCheckFromTerms(terms):
                    break
                else:
                    exp = terms[0]
                    for i in range(1, len(terms)):
                        exp = exp + terms[i]
                    Equation._logger.info("Equation: " + str(exp) + " has failed sanityCheck!")

        Equation.AssignRandomSignsToTerms(terms)

        exp = terms[0]
        for i in range(1, len(terms)):
            exp = exp + terms[i]

        eqn = Equation(exp)
        eqn.divideByCommonUnnamedConstants()

        return eqn

    #Will eventually delete this next method - once the new one is sufficently tested
    @classmethod
    def GetUofMToPrimitiveTermLookupTable_old(cls, vars, derivatives, constants, max_power=3,
                                          max_vars_derivatives_and_constants_per_term=4):
        # This impelementation is very slow. Can be sped up drastically by first choosing the up to
        # max_vars_derivatives_and_constants_per_eqn slots and filling them with integers between 0
        # and max_power
        uOfMToTermLookupDict = dict()

        for i in range(1, pow(max_power + 1, len(vars) + len(derivatives) + len(constants))):
            baseNum = Equation.ToBase(i, max_power + 1)
            s_baseNum = str(baseNum)
            if len(s_baseNum) < len(vars) + len(derivatives) + len(constants):
                s_baseNum = '0' * (len(vars) + len(derivatives) + len(constants) - len(s_baseNum)) + s_baseNum
            numZeroes = s_baseNum.count('0')
            numVarsDerivsAndConstants = len(vars) + len(derivatives) + len(constants) - numZeroes
            if numVarsDerivsAndConstants <= max_vars_derivatives_and_constants_per_term:
                primitiveTerm = Equation.GenerateTermFromBaseNum(baseNum, vars, derivatives, constants)
                u_of_m = Equation.GetUofMForTerm(primitiveTerm)
                lookupPair = uOfMToTermLookupDict.get(u_of_m._units)
                # (listOfPrimitiveTerms, vars_and_derivatives_in_use) = lookupDict.get(u_of_m._units)
                vars_derivs_and_constants_in_term = primitiveTerm.free_symbols
                if lookupPair is None:
                    uOfMToTermLookupDict[u_of_m._units] = ([primitiveTerm], vars_derivs_and_constants_in_term)
                else:
                    listOfPrimitiveTerms = lookupPair[0]
                    vars_derivatives_and_constants_in_use = lookupPair[1]
                    listOfPrimitiveTerms.append(primitiveTerm)
                    vars_derivatives_and_constants_in_use.update(vars_derivs_and_constants_in_term)

        # now sort by number of # of vars and derivs in use. This may not be necessary
        # sortedDict = sorted(uOfMToTermLookupDict.items(), key=lambda item: len(item[1][1]), reverse=True)
        # orphan_uofms = set()
        # for key, value in uOfMToTermLookupDict.items():
        #    if len(value[1]) == 1:
        #        orphan_uofms.add(key)
        # for orphan in orphan_uofms:
        #    del uOfMToTermLookupDict[orphan]

        return uOfMToTermLookupDict

    @classmethod
    def GetUofMToPrimitiveTermLookupTable(cls, vars, derivatives, constants, max_power=3,
                                          max_vars_derivatives_and_constants_per_term=4):
        uOfMToTermLookupDict = dict()
        vars_derivs_and_constants = set(vars).union(set(derivatives)).union(set(constants))
        for i in range(1, max_vars_derivatives_and_constants_per_term+1):
            sBaseNums = []
            for j in range(pow(max_power, i)):  #changed max_power+1 => max_power
                baseNum = Equation.ToBase(j, max_power) #changed max_power+1 => max_power
                for k in range(i):
                    baseNum += pow(10,k)
                s_baseNum = str(baseNum)
                sBaseNums.append(s_baseNum)
            tuples = itertools.combinations(vars_derivs_and_constants, i)
            for tuple in tuples:
                for s_baseNum in sBaseNums:
                    primitiveTerm = Equation.GenerateTermFromTupleAndBaseNum(tuple, s_baseNum)
                    u_of_m = Equation.GetUofMForTerm(primitiveTerm)
                    lookupPair = uOfMToTermLookupDict.get(u_of_m._units)
                    vars_derivs_and_constants_in_term = primitiveTerm.free_symbols
                    if lookupPair is None:
                        uOfMToTermLookupDict[u_of_m._units] = ([primitiveTerm], vars_derivs_and_constants_in_term)
                    else:
                        listOfPrimitiveTerms = lookupPair[0]
                        vars_derivatives_and_constants_in_use = lookupPair[1]
                        listOfPrimitiveTerms.append(primitiveTerm)
                        vars_derivatives_and_constants_in_use.update(vars_derivs_and_constants_in_term)

        return uOfMToTermLookupDict


    @classmethod
    def GenerateRandomDimensionallyConsistentEquationWithSpecifiedVarOrDerivative(cls, vars, derivatives, constants,
                                                                                  u_of_mToTermLookupDict,
                                                                                  given_var=None, given_derivative=None,
                                                                                  given_constant=None,
                                                                                  max_power=3,
                                                                                  max_vars_derivatives_and_constants_per_eqn=NO_MAX):
        num_terms = Equation.GetRandomNumberOfTerms()

        newFirstTerm = None
        allTermsForGivenVarDerivOrConstant = Equation.GetAllTermsWithGivenVarDerivativeOrConstant(vars=vars, derivatives=derivatives,
                                                                                  constants=constants,
                                                                                  given_var=given_var,
                                                                                  given_derivative=given_derivative,
                                                                                  given_constant=given_constant,
                                                                                  max_power=max_power, max_varsDerivsAndConstants=4)

        terms_to_keep = []
        while True:
            for term in allTermsForGivenVarDerivOrConstant:
                u_of_m_for_term = Equation.GetUofMForTerm(term)
                termsForUofM = u_of_mToTermLookupDict.get(u_of_m_for_term._units)[0]
                if len(termsForUofM) >= num_terms:
                    if Equation.GetCommonFactors(termsForUofM) == set():
                        terms_to_keep.append(term)

            if len(terms_to_keep) == 0 and num_terms > 2:
                num_terms -= 1
            elif len(terms_to_keep) == 0 and num_terms == 2:
                return None
            else:
                break

        newFirstTerm = None
        strippedNewFistTerm = None
        foundGoodFirstTerm = False
        eqn = None
        while True: #To get a viable equation
            for i in range(100):  #To get a viable first tierm
                newFirstTerm = Equation.GenerateRandomTermWithGivenVarDerivativeOrConstant(vars=vars,
                                                                derivatives=derivatives,
                                                                constants=constants,
                                                                givenVar=given_var,
                                                                givenDeriv=given_derivative,
                                                                givenConstant=given_constant,
                                                                max_power=max_power,
                                                                max_varsDerivsAndConstants=4)
                strippedNewFistTerm = Equation.GetUnnamedConstantStrippedTerm(newFirstTerm)
                if strippedNewFistTerm in terms_to_keep:
                    Equation._logger.info("Term is a good one!")
                    foundGoodFirstTerm = True
                    break
                else:
                    Equation._logger.info("Term is not good. Trying again....")

            if not foundGoodFirstTerm:
                randIndex = random.randint(0, len(terms_to_keep) - 1)
                strippedNewFistTerm = terms_to_keep[randIndex]
                c = Equation.GenerateRandomUnnamedConstant()
                newFirstTerm = c * strippedNewFistTerm

            u_of_m_for_term = Equation.GetUofMForTerm(strippedNewFistTerm)
            termsForUofM = u_of_mToTermLookupDict.get(u_of_m_for_term._units)[0].copy()
            termsForUofM.remove(strippedNewFistTerm)
            existing_terms = [newFirstTerm]
            addtional_terms = Equation.GenerateFixedNumberofDimensionallyConsistentTermsFromList(num_terms=num_terms-1,
                                                            candidateTerms=termsForUofM,
                                                            existing_terms=existing_terms,
                                                            max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
            if addtional_terms == []:
                continue

            existing_terms.extend(addtional_terms)
            Equation.AssignRandomSignsToTerms(existing_terms)
            exp = existing_terms[0]
            for i in range(1, len(existing_terms)):
                exp = exp + existing_terms[i]
            eqn = Equation(exp)
            if eqn.sanityCheck():
                break
            else:
                Equation._logger.info("Equation: " + str(exp) + " has failed sanityCheck!")

        eqn.divideByCommonUnnamedConstants()

        return eqn



    def getSymbolsUsed(self):  #returns a set of variables, derivs and constants
        return self._poly.free_symbols

    #def getUofM(self, sym):
    #    for var in self._variables:
    #        if var.name == str(sym):
    #            return var._u_of_m
    #    for deriv in self._derivatives:
    #        if deriv.name == str(sym):
    #            return deriv._u_of_m
    #    for const in self._constants:
    #        if const.name == str(sym):
    #            return const._u_of_m
    #
    #    return None

    @classmethod
    def GetAllTermsWithGivenVarDerivativeOrConstant(cls, vars=[], derivatives=[], constants=[], given_var=None, given_derivative=None,
                                            given_constant=None, max_power=3, max_varsDerivsAndConstants=NO_MAX):
        allTerms = []
        givenVarDerivOrConstantIndex = -1
        if given_var is not None:
            givenVarDerivOrConstantIndex = vars.index(given_var)
        elif given_derivative is not None:
            givenVarDerivOrConstantIndex = len(vars) + derivatives.index(given_derivative)
        elif given_constant is not None:
            givenVarDerivOrConstantIndex = len(vars) + len(derivatives) + constants.index(given_constant)

        for i in range(1, pow(max_power + 1, len(vars) + len(derivatives) + len(constants))):
            baseNum = Equation.ToBase(i, max_power + 1)
            s_baseNum = str(baseNum)
            if len(s_baseNum) < len(vars) + len(derivatives) + len(constants):
                s_baseNum = '0' * (len(vars) + len(derivatives) + len(constants) - len(s_baseNum)) + s_baseNum

            if s_baseNum[givenVarDerivOrConstantIndex] != '0':
                numZeroes = s_baseNum.count('0')
                numVarserivsAndConstants = len(vars) + len(derivatives) + len(constants) - numZeroes
                if numVarserivsAndConstants <= max_varsDerivsAndConstants:
                    term = Equation.GenerateTermFromBaseNum(baseNum, vars, derivatives, constants)
                    allTerms.append(term)

        return allTerms

    @classmethod
    def GenerateRandomTermWithGivenVarDerivativeOrConstant(cls, vars, derivatives=[], constants=[], givenVar=None, givenDeriv=None,
                                                           givenConstant=None, max_power=3,
                                                           max_varsDerivsAndConstants=4):
        #allTerms = Equation.GetAllTermsWithGivenVarOrDerivative(vars=vars, derivatives=derivatives, given_var=givenVar,
        #                                             given_derivative=givenDeriv, max_power=3, max_varsAndDerivs=3)
        #rand_index = random.randint(0, len(allTerms)-1)
        #return allTerms[rand_index]  #Need to do better than this! Need to verify there are enogh terms available with the same unit of measure!
        givenVarDerivOrConstant = None
        if givenVar is not None:
            givenVarDerivOrConstant = givenVar
        elif givenDeriv is not None:
            givenVarDerivOrConstant = givenDeriv
        elif givenConstant is not None:
            givenVarDerivOrConstant = givenConstant

        Equation._logger.info("Trying to generate a random term for var-deriv-or-constant: " + str(givenVarDerivOrConstant))
        if givenVarDerivOrConstant is None:
            Equation._logger.error("givenVarDerivOrConstant is None in Equation.GenerateRandomTermWithGivenVarDerivativeOrConstant!!!")

        while True:
            term = Equation.GenerateRandomTerm(vars=vars, derivatives=derivatives, constants=constants, max_power=max_power,
                                                           max_varsDerivsAndConstants=max_varsDerivsAndConstants)
            if givenVarDerivOrConstant in term.free_symbols:
                Equation._logger.info("Term found: " + str(term))
                return term

    @classmethod
    def GenerateRandomTerm(cls, vars, derivatives=[], constants=[], max_power=3,
                                                           max_varsDerivsAndConstants=4):
        vars_derivs_and_constants = []
        vars_derivs_and_constants.extend(vars)
        vars_derivs_and_constants.extend(derivatives)
        vars_derivs_and_constants.extend(constants)

        term = None
        while True: #Loop until max_ppwer and max_varsDerivsAndConstants are satisfied
            num_factors = len(Equation.PROB_FACTORS_PER_TERM)
            prob_so_far = 0.0
            r = random.random()
            for i in range(len(Equation.PROB_FACTORS_PER_TERM)):
                prob_so_far += Equation.PROB_FACTORS_PER_TERM[i]
                if r < prob_so_far:
                    num_factors = i+1
                    break

            vars_derivs_and_constants_so_far = []  #An array indicating counts of the vars, derivatives and constants
                                                   # that are already being used
            var_deriv_and_constant_probs = []
            for i in range(len(vars_derivs_and_constants)):
                vars_derivs_and_constants_so_far.append(0)
                var_deriv_and_constant_probs.append(1.0/len(vars_derivs_and_constants))

            for i in range(num_factors):
                r = random.random()
                var_deriv_or_constant = var_deriv_and_constant_probs[len(var_deriv_and_constant_probs) - 1]
                prob_so_far = 0.0
                for j in range(len(vars_derivs_and_constants)):
                    prob_so_far += var_deriv_and_constant_probs[j]
                    if r < prob_so_far:
                        if vars_derivs_and_constants_so_far[j] == 0:
                            old_prob = var_deriv_and_constant_probs[j]
                            if isinstance(var_deriv_and_constant_probs, Variable):  #Constants count here
                                var_deriv_and_constant_probs[j] *= Equation.SAME_FACTOR_VARIABLE_BIAS
                            elif isinstance(var_deriv_and_constant_probs, Derivatif): #should be the only other case
                                var_deriv_and_constant_probs[j] *= Equation.SAME_FACTOR_DERIVATIVE_BIAS
                            aggregate_weight = 1 + var_deriv_and_constant_probs[j] - old_prob
                            # must now rebalance the probs
                            for k in range(len(vars_derivs_and_constants)):
                                var_deriv_and_constant_probs[k] = var_deriv_and_constant_probs[k]/aggregate_weight
                        else:
                            if isinstance(var_deriv_and_constant_probs, Variable):  #prob bounces back down to original prob
                                var_deriv_and_constant_probs[j] /= Equation.SAME_FACTOR_VARIABLE_BIAS
                        vars_derivs_and_constants_so_far[j] += 1
                        break

            term = None
            highest_power = 1
            for i in range(len(vars_derivs_and_constants)):
                if vars_derivs_and_constants_so_far[i] == 1:
                    if term is not None:
                        term = term*vars_derivs_and_constants[i]
                    else:
                        term = vars_derivs_and_constants[i]
                elif vars_derivs_and_constants_so_far[i] > 1:
                    if term is not None:
                        term = term*vars_derivs_and_constants[i]**vars_derivs_and_constants_so_far[i]
                    else:
                        term = vars_derivs_and_constants[i]**vars_derivs_and_constants_so_far[i]
                    if vars_derivs_and_constants_so_far[i] > highest_power:
                        highest_power = vars_derivs_and_constants_so_far[i]

            if highest_power <= max_power and len(term.free_symbols) <= max_varsDerivsAndConstants:
                break


        c = Equation.GenerateRandomUnnamedConstant()
        term = c*term

        return term


    @classmethod
    def GenerateRandomUnnamedConstant(cls):  #should also test for common constant factors now
        c = 1
        r = random.random()
        if r < Equation.PROB_TERM_HAS_NON_UNITAL_INTEGER_CONSTANT:
            r = random.random()
            prob_so_far = 0.0
            for i in range(len(Equation.PROB_OF_SMALL_INTEGER_CONSTANTS)):
                prob_so_far += Equation.PROB_OF_SMALL_INTEGER_CONSTANTS[i]
                if r < prob_so_far:
                    c = i
                    break

        return c

    @classmethod
    def GetCommonFactors(cls, terms):
        common_factors = terms[0].free_symbols
        for i in range(1, len(terms)):
            factors = terms[i].free_symbols
            common_factors = common_factors.intersection(factors)

        return common_factors


    @classmethod
    def TermAmongExistingTerms(cls, existingTerms, term):
        for t in existingTerms:
            if t == term:
                return True

        return False

    @classmethod
    def AssignRandomSignsToTerms(cls, terms):  #guarantees at least one term has a positive and one term is negative
        #In rare cases just output an all positive (homogeneous) equation
        if random.random() < Equation.HOMOGENEOUS_EQN_PROB:
            return

        #First pick a term to be positive and a different one to be negative
        rands = []
        signs = []
        max_rand = -1.0
        min_rand = 2.0
        for i in range(len(terms)):
            rands.append(random.random())
            if rands[i] > max_rand:
                max_rand = rands[i]
            if rands[i] < min_rand:
                min_rand = rands[i]

        num_above_the_mean = 0
        mean = (max_rand + min_rand)/2.0
        for i in range(len(terms)):
            if rands[i] < mean:
                signs.append(-1)
            else:
                signs.append(1)
                num_above_the_mean += 1

        if num_above_the_mean < len(terms)/2.0: #flip the signs so the majority are positive
            for i in range(len(terms)):
                signs[i] = -1*signs[i]

        for i in range(len(terms)):
            if signs[i] == -1:
                terms[i] = -terms[i]

        return None

    def sanityCheck(self): #perform various different sanity checks to make sure an equation is of an appropriate form
        #First check is whether equation is just one var-deriv-or-constant is equal (mod unnamed constants) to another
        terms = self.getTerms()
        return Equation.SanityCheckFromTerms(terms)

    @classmethod
    def SanityCheckFromTerms(cls, terms):
        #check that there are not two terms, each with a single free symbol
        if len(terms) == 2 and len(terms[0].free_symbols) == 1 and len(terms[1].free_symbols) == 1:
            return False
        #check that no two terms are equal
        if len(terms) != len(set(terms)):
            return False
        else:
            for pair in itertools.combinations(set(terms), 2):
                if Equation.TermsEqualModUnnamedConstants(pair[0], pair[1]):
                    return False
        return True



    def getTerms(self):
        return self._poly.expr.args

    def __eq__(self, eqn):
        return self.getTerms() == eqn.getTerms()

    def __hash__(self):
        return hash(str(self))


    def equalModUnnamedConstants(self, eqn):
        ourTerms = self.getTerms()
        ourStrippedTerms = set()
        for term in ourTerms:
            ourStrippedTerms.add(Equation.GetUnnamedConstantStrippedTerm(term))
        otherTerms = eqn.getTerms()
        otherStrippedTerms = set()
        for term in otherTerms:
            otherStrippedTerms.add(Equation.GetUnnamedConstantStrippedTerm(term))

        return ourStrippedTerms == otherStrippedTerms

    @classmethod
    def TermsEqualModUnnamedConstants(cls, term1, term2):
        t1 = Equation.GetUnnamedConstantStrippedTerm(term1)
        t2 = Equation.GetUnnamedConstantStrippedTerm(term2)

        return t1 == t2

    @classmethod
    def TermsEqual(self, term1, term2):
        return term1 == term2

    @classmethod
    def GetUnnamedConstantStrippedTerm(cls, term):
        t = term
        if len(t.args) > 0 and isinstance(t.args[0], Integer):
            return t / t.args[0]
        else:
            return t

    @classmethod
    def GetUnnamedConstantForTerm(cls, term):
        if len(term.args) > 0 and isinstance(term.args[0], Integer):
            return term.args[0]
        else:
            return Integer(1)

    def divideByCommonUnnamedConstants(self):
        terms = self.getTerms()
        firstTerm = terms[0]
        firstConstant = Equation.GetUnnamedConstantForTerm(firstTerm)
        for i in range(1, len(terms)):
            term = terms[i]
            constant = Equation.GetUnnamedConstantForTerm(term)
            if constant != firstConstant and constant != -1*firstConstant:
                return
        if firstConstant < 0:
            firstConstant = -1*firstConstant


        exp = None
        for term in terms:
            if exp is None:
                exp = term/firstConstant
            else:
                exp += term/firstConstant

        self._poly = Poly(exp)




    @classmethod
    def SetLogging(cls, filename='theorizer.log', filemode='w', encoding='utf-8', level=logging.DEBUG):
        logging.basicConfig(filename=filename, filemode=filemode, encoding=encoding, level=level)

    def __str__(self):
        return str(self._poly.expr) + " (U of M: " + str(self.getUofM()) + ")"














