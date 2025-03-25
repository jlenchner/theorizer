"""The Equation class, including implementations for dimensionally consistent
   and not necessarily dimensionally consistent equations. Equations are also
   sometimes referred to as "axioms" in our papers.
"""

# Author: Jonathan Lenchner (lenchner@us.ibm.com)
#
# License: BSD 3-Clause

from sympy import *
from variable import *
from derivative import *
from constant import *
import random
import math
import logging
import itertools
from deprecated import deprecated  # In PyCharm debug mode this can cause the error:
                                   # ValueError: wrapper has not been initialized.
                                   # To remedy, goto Help > Find Action > Registry
                                   # and uncheck the box for python.debug.low.impact.monitoring.api.



class Equation:
    PROB_OF_NUM_TERMS = [0, 42.0/70.0, 20.0/70.0, 6.0/70.0, 2/70.0]  #Probs of a randomly generated poly having 1 term, 2 terms, etc.
                                                                     #Values have been extrapolated from the AI Feynman database.
    PROB_FACTORS_PER_TERM = [0.19, 0.47, 0.25, 0.09]  #Probs that a term in a randomly generated poly will have 1, 2, 3, etc. factors
                                                      #Values have been extrapolated from the AI Feynman database.
    SAME_VARIABLE_FACTOR_BIAS = 2.0  #Bias for repeating the same variable in a term (in other words, favoring x^2 over xy). Constants,
                                     # since they are variables, have the same bias.
    SAME_DERIVATIVE_FACTOR_BIAS = 0.5 #Bias for repeating the same derivative in a term (in other words, favoring dxdt^2 over dxdt*dydt)

    PROB_TERM_HAS_NON_UNITAL_INTEGER_CONSTANT = 0.2
    PROB_OF_SMALL_INTEGER_CONSTANTS = [0.0, 0.0, 15.0/34.0, 10.0/34.0, 9.0/34.0]  #Probs that a small non-unital integer constant is [0,1,2,3,4,5]
                                                                                  #Values have been extrapolated from the AI Feynman database.

    HOMOGENEOUS_EQN_PROB = 0.05

    NO_MAX = 999

    _logger = logging.getLogger()


    def __init__(self, exp):
        """
        :param exp: Is the Equation expressed in polynomial form. All component
                    Variables, Derivatifs and Constants must be pre-constructed.
                    Sample usage:
                    d1, d2, m1, m2, Fg = variables('d1, d2, m1, m2, Fg')
                    G = Constant('G')
                    eqn = Equation(Fg*d1**2 + 2*Fg*d1*d2 + Fg*d2**2 - G*m1*m2)
        """
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
        """
            Returns the unit of measure (class UofM) associated with a given term.
        """
        return UofM(Equation.InnerGetUofMForTerm(term))

    @classmethod
    def InnerGetUofMForTerm(cls, term):
        """
            This method should NOT be called directly. It is used internally
            by GetUofMForTerm().
        """
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
        """
            Returns the UofM for the Equation if the Equation is dimensionally
            consistent and None otherwise.
        """
        terms = self.getTerms()
        firstUofM = Equation.GetUofMForTerm(terms[0])
        for i in range(1, len(terms)):
            UofMforTerm = Equation.GetUofMForTerm(terms[i])
            if firstUofM != UofMforTerm:
                return None

        return firstUofM


    @classmethod
    def InferVarsDerivativesAndConstantsFromExpression(cls, exp):
        """
            Determines the Variablees, Derivatifs and Constants from
            a symbolic expression. This method should not be needed --
            its main job is as a helper method for Equation.__init__().
        """
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
        """
            Returns a random number of terms for randomly generated
            Equation using the values in Equation.PROB_OF_NUM_TERMS,
            which were in turn obtained from the AI Feynman database.
        """
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
    def GenerateRandom(cls, vars=[], derivatives=[], constants=[],
                       max_vars_derivatives_and_constants_per_eqn=NO_MAX):
        """
            Generates a random, not necessarily dimensionally consistent Equation,
            given provided variables (class Variable), deriviatives (class Derivatif)
            and constants (class Constant) as well as an integer value for
            max_vars_derivatives_and_constants_per_eqn

            :param vars: Provided set of variables (each of class Variable)
            :param derivatives: Provided set of deriviatives (each of class Derivatif)
            :param constants: Provided set of constants (each of class Constant)
            :param max_vars_derivatives_and_constants_per_eqn: the maximum number
                    variables, derivatives and constants that should appear in
                    the equation.
            :return: Equation
        """
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
    @deprecated(reason="This method is deprecated and will run very slowly for large numbers \
                    of varibales, derivatives and constants.  Instead, we recommend using \
                    EquationSystem._LookupDict to obtain terms dimensionally consistent \
                    with a given term.")
    def GetAllDimensionallyConsistentTerms(cls, term_to_match, vars=[], derivatives=[], constants=[], max_power=3,
                                           max_varsDerivsAndConstants=NO_MAX):
        """
            Get all terms dimensionally consistent with a provided term.
            :param term_to_match: The term you wish to match dimensionally
            :param vars: List of available variables (class Variable)
            :param derivatives: List of available derivatives (class Derivatif)
            :param constants: List of available constants (class Constant)
            :param max_power: An integer designated the maximum power that can appear
                            in any term associated with a variable, derivative or constant.
            :param max_varsDerivsAndConstants: The maximum number of distinct variables,
                            derivatives and constants that can appear in any one term.
            :return: A list of terms that are dimensionally consistent with the provided terms
        """
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
        """
            Given a non-negative integer, returns the analog of the integer in
            a given base. E.G.:
                Equation.ToBase(32, 3) = 1012, Equation.ToBase(41, 2) = 101001
            :param num: Integer to convert to the given base
            :param base: Positive integer indicating the give base
            :return: Integer equivalent of num in the provided base
        """
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
        """
            Used internally by the method Equation.GetUofMToPrimitiveTermLookupTable(). Takes a tuple
            of variables, derivatives and constants, i.e., something of the form (vdc_1,...,vdc_k),
            where each vdc_i is either a variable, derivative or constant, together with a given string
            of length k, each of whose charcters are in the set {'0', '1',...,'9'}, and returns a term
            of the form vdc_1**s_baseNum[0]*...*vdc_k**s_baseNum[k-1]
            :param tuple: A tuple of the form (vdc_1,...,vdc_k), where each vdc_i is either a variable,
                        derivative or constant
            :param s_baseNum: A string of the same length as the tuple size, each character of whic is
                        in the set {'0','1',...,'9'}
            :return: The resultant term
        """
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
        """
            Used internally by the deprecated method Equation.GenerateTermFromBaseNum
            :param baseNum: Positive integer indicating what powers to raise each of the
                        supplied variables, derivatives and constants to, in order to
                        produce the desired term.
            :param vars: List of supplied variables (class Variable)
            :param derivatives: List of supplied derivatives (class Derivatif)
            :param constants: List of supplied constants (class Constant)
            :return: Generated term
        """
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
    def GetImpliedSelectionProbabilitiesForCandidateTerms(cls, candidateTerms, sigDict):
        """
            Given a set of candidate terms (typically all dimensionally consistent), returns
            the implicit probability of randomly choosing each term, based on the signature of
            their exponents
            :param candidateTerms: A list of candidate terms to be chosen from
            :param sigDict: A pointer to the EquationSystem._sigDict dictionary, which for each
                        variable-derivative-constant signature contains a number between 1 and
                        10,000, representing their relative frequency.
            :return: A list of implied probabilities associated with each term. The probabilities
                        sill always sum to 1.
        """
        selectionCounts = []
        selectionProbabilities = []

        totalCount = 0
        for term in candidateTerms:
            sig = Equation.GetVDCSignatureForTerm(term)
            pre_count = sigDict.get(sig)
            if pre_count  is None:
                pre_count = 0
            count = max(pre_count, 1)
            totalCount += count
            selectionCounts.append(count)
        for i in range(len(selectionCounts)):
            selectionProbabilities.append(selectionCounts[i]/totalCount)

        return selectionProbabilities

    @classmethod
    def RemoveProbByIndexAndRenormalize(cls, termSeectionProbs, index):
        """
            Given a set of probabilities, removes the element at a given index
            and renormalizes the remaining elements
            :param termSeectionProbs: The given list of probabilities
            :param index: The index of the element to remove
            :return: The resulting array of probabilities. Will always sum to 1 and
                    contain one less element than the starting array.
        """
        cumProbWithoutIndex = 1.0 - termSeectionProbs[index]
        for i in range(len(termSeectionProbs)):
            if i != index:
                termSeectionProbs[i] = termSeectionProbs[i]/cumProbWithoutIndex

        del termSeectionProbs[index]


    @classmethod
    def PickRandomTermFromListGivenProbs(cls, candidateTerms, selectionProbs):
        """
            Picks a given term from a list given a parallel list of selection probabilities
            :param candidateTerms: List of terms to select from
            :param selectionProbs: List of selection probabilities
            :return: Selected term
        """
        r = random.random()
        probSoFar = 0
        for i in range(len(candidateTerms)):
            probSoFar += selectionProbs[i]
            if r < probSoFar:
                return candidateTerms[i]

        return candidateTerms[len(candidateTerms)-1]

    @classmethod
    def GenerateFixedNumberofDimensionallyConsistentTermsFromList(cls, num_terms, candidateTerms, existing_terms=[],
                                                                  max_vars_derivatives_and_constants_per_eqn=NO_MAX,
                                                                  sigDict=None):
       """
           Returns a specified number of dimensionally consistent terms from a candidate list
           of assumed to be dimensionally consistent candidates, beginning from a specified
           list of existing terms, which are also assumed to be dimensionally consistent. The
           maximum number of variables, derivatives and constants used across all terms is
           guaranteed not to exceed a specified maximum, if such a maximum is supplied.
           :param num_terms: Number of terms required, a positive integer
           :param candidateTerms: The candidate terms to be chosen from (should not include
                            the existing terms. Assumed to be dimensionally consistent with
                            the existing terms, or if there are no existing terms, they should
                            be dimensionally consistent with each other)
           :param existing_terms: The existing terms The existing terms, assumed to be dimen-
                            sionally consistent with each other.
           :param max_vars_derivatives_and_constants_per_eqn: Maximum number of distinct variables,
                            derivatives and constants to be used across all terms.
           :param sigDict: A pointer to EquationSystem._sigDict. An exception will be raised if
                            this is not passed.
           :return:Returns a List of randomly selected candidate terms.
       """
       if sigDict is  None:
           raise ValueError("A valid sigDict must be passed in call to Equation.GenerateFixedNumberofDimensionallyConsistentTermsFromList()!")

       termSelectionProbabilities = Equation.GetImpliedSelectionProbabilitiesForCandidateTerms(candidateTerms, sigDict)
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
            termSelectionProbabilitiesCopy = termSelectionProbabilities.copy()
            terms = []
            vars_derivs_and_constants_in_use = set()
            for term in existing_terms:
                vars_derivs_and_constants_in_use.update(term.free_symbols)

            for i in range(num_terms):
                while True:
                    # OLDER IMPLEMENTATION - PICKING TERMS UNIFORMLY AT RANDOM
                    # FROM THE LIST OF DIMENSIONALLY CONSISTENT CANDIDATES
                    #rand_int = random.randint(0, len(candidateTermsCopy) - 1)
                    #term = candidateTermsCopy[rand_int]
                    #vars_derivs_and_constants_in_term = term.free_symbols
                    #new_vars_derivs_and_constants = vars_derivs_and_constants_in_term - vars_derivs_and_constants_in_use
                    #if len(new_vars_derivs_and_constants) + len(vars_derivs_and_constants_in_use) <= max_vars_derivatives_and_constants_per_eqn:
                    #    vars_derivs_and_constants_in_use = vars_derivs_and_constants_in_use.union(new_vars_derivs_and_constants)
                    #    c = Equation.GenerateRandomUnnamedConstant()
                    #    term = c * term
                    #    terms.append(term)
                    #    candidateTermsCopy.remove(Equation.GetUnnamedConstantStrippedTerm(term))
                    #    break

                    term = Equation.PickRandomTermFromListGivenProbs(candidateTermsCopy, termSelectionProbabilitiesCopy)
                    term_index = candidateTermsCopy.index(term)
                    vars_derivs_and_constants_in_term = term.free_symbols
                    new_vars_derivs_and_constants = vars_derivs_and_constants_in_term - vars_derivs_and_constants_in_use
                    if len(new_vars_derivs_and_constants) + len(
                            vars_derivs_and_constants_in_use) <= max_vars_derivatives_and_constants_per_eqn:
                        vars_derivs_and_constants_in_use = vars_derivs_and_constants_in_use.union(
                            new_vars_derivs_and_constants)
                        c = Equation.GenerateRandomUnnamedConstant()
                        term = c * term
                        terms.append(term)
                        candidateTermsCopy.remove(Equation.GetUnnamedConstantStrippedTerm(term))
                        Equation.RemoveProbByIndexAndRenormalize(termSelectionProbabilitiesCopy, term_index)
                        break

            all_terms = []
            all_terms.extend(existing_terms)
            all_terms.extend(terms)
            if Equation.GetCommonFactors(all_terms) == set():
                return terms


    @classmethod
    def GenerateRandomDimensionallyConsistent(cls, vars, derivatives, constants, u_of_mToTermLookupDict=None,
                                              sigDict=None,
                                              max_power=3,
                                              max_vars_derivatives_and_constants_per_eqn=NO_MAX):
        """
            Generates a random dimensionally consistent Equation, given lists of variables
            (class Variable), derivatives (class Derivatif) and constants (class Constant),
            together with pointers to EquationSystem._lookupDict (a dictionary that given a
            unit of measure (class UofM) returns all terms with that unit of measure) and
            EquationSystem._sigDict (a dictionary that given a variable-derivative-constant
            exponent "signature", returns the frequency count for randomly picking a term with
            that signature), along with other parameters.
            :param vars: a list of variables (class Variable) to use in constructing the
                        Equation
            :param derivatives: a list of variables (class Variable) to use in constructing
                        the Equation
            :param constants: a list of constants (class Constant) to use in constructing
                        the Equation
            :param u_of_mToTermLookupDict: a pointer to EquationSystem._lookupDict (a dictionary
                        that given a unit of measure (class UofM) returns all terms with
                        that unit of measure). If None is passed, the dictionary will be created.
            :param sigDict: a pointer to EquationSystem._sigDict (a dictionary that given a
                        variable-derivative-constant exponent "signature", returns the frequency
                        count for randomly picking a term with that signature).  If None is passed,
                        the dictionary will be created.
            :param max_power: a non-negative integer specifying the maximum power to raise any
                        variable, derivative or constant to in the returned equation
            :param max_vars_derivatives_and_constants_per_eqn: a non-negative integer specifying
                        the maximum number of distinct variables, derivatives and constants allowed
                        in the Equation.
            :return: the randomly chosen Equation
        """
        if max_vars_derivatives_and_constants_per_eqn == Equation.NO_MAX:
            max_vars_derivatives_and_constants_per_eqn = len(vars) + len(derivatives) + len(constants)

        if u_of_mToTermLookupDict is None:
            u_of_mToTermLookupDict = Equation.GetUofMToPrimitiveTermLookupTable(vars=vars, derivatives=derivatives,
                                                constants=constants, max_power=max_power)
        if sigDict is None:
            sigDict = Equation.GenerateVDCSigDistributionDict(vars=vars, derivatives=derivatives, constants=constants,
                                                              max_power=max_power, max_varsDerivsAndConstants=4)

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
                                                                                                 existing_terms=terms,
                                                                                                 sigDict=sigDict)
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


    @classmethod
    def GetUofMToPrimitiveTermLookupTable(cls, vars, derivatives, constants, max_power=3,
                                          max_vars_derivatives_and_constants_per_term=4):
        """
            Generates the Unit of Measure to Primitive Term (meaning term without an unnamed
            constant) lookup dictionary, given a set of variables, derivatives and constants,
            and other parameters. The result is stored for the given set of variables, derivatives
            and constants in EquationSystem._lookupDict.
            :param vars: List of variables (class Variable) to use
            :param derivatives: List of derivatives (class Derivatif) to use
            :param constants: List of constants (class Constant) to use
            :param max_power: Non-negative integer specifying the maximum power to  raise any
                        variable, derivative or constant to in a term.
            :param max_vars_derivatives_and_constants_per_term: Non-negative integer specifying
                        the maximum number of distinct variables, derivatives and constants to
                        allow in any term.
            :return: The lookup dictionary
        """
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
    def GenerateRandomDimensionallyConsistentEquationWithSpecifiedVarDerivOrConstant(cls,
                                    vars, derivatives, constants,
                                    u_of_mToTermLookupDict,
                                    sigDict,
                                    given_var=None, given_derivative=None,
                                    given_constant=None,
                                    max_power=3,
                                    max_vars_derivatives_and_constants_per_eqn=NO_MAX):
        """
            Generate a random dimensionally consistent Equation starting from a given variable,
            derivative or constant. Note that precisely one of given_var, given_derivative and
            given_constant should be specified.
            :param vars: List of variables (class Variable) to use
            :param derivatives: List of derivatives (class Derivatif) to use
            :param constants: List of constants (class Constant) to use
            :param u_of_mToTermLookupDict: a pointer to EquationSystem._lookupDict (a dictionary
                        that given a unit of measure (class UofM) returns all terms with
                        that unit of measure). If None is passed, the dictionary will be created.
            :param sigDict: a pointer to EquationSystem._sigDict (a dictionary that given a
                        variable-derivative-constant exponent "signature", returns the frequency
                        count for randomly picking a term with that signature).  If None is passed,
                        the dictionary will be created.
            :param given_var: the given variable (class Variable) or None
            :param given_derivative: the given derivative (class Derivatif) or None
            :param given_constant: the given constant (class Constant) or None
            :param max_power: a non-negative integer specifying the maximum power to raise any
                        variable, derivative or constant to in the returned equation
            :param max_vars_derivatives_and_constants_per_eqn: a non-negative integer specifying
                        the maximum number of distinct variables, derivatives and constants allowed
                        in the Equation
            :return: the randomly chosen Equation
        """
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
                                                            max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn,
                                                            sigDict=sigDict)
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



    def getSymbolsUsed(self):
        """
            Returns a set of the variables, derivatives and constants used
        """
        return self._poly.free_symbols

    @classmethod
    def GetAllTermsWithGivenVarDerivativeOrConstant(cls, vars=[], derivatives=[], constants=[], given_var=None, given_derivative=None,
                                            given_constant=None, max_power=3, max_varsDerivsAndConstants=NO_MAX):
        """
            Returns all terms with the given variable (class Variable), derivative (class Derivatif) or
            constant (class Constant), subject to the constraints that the maximum power associated
            with any of these is max_power and the maximum number of distinct variables, derivatives or
            constants is max_varsDerivsAndConstants. Exactly one of the parameters given_var,
            given_derivative and given_constant should be specified, others should have the value None.
            :param vars: List of variables (class Variable) to use
            :param derivatives: List of derivatives (class Derivatif) to use
            :param constants: List of constants (class Constant) to use
            :param given_var: the given variable (class Variable) or None
            :param given_derivative: the given derivative (class Derivatif) or None
            :param given_constant: the given constant (class Constant) or None
            :param max_power: a non-negative integer specifying the maximum power to raise any
                        variable, derivative or constant to in any term
            :param max_varsDerivsAndConstants: a non-negative integer specifying
                        the maximum number of distinct variables, derivatives and constants allowed
                        in any term
            :return: A list of all terms satisfying to conditions specified
        """
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
    def GenerateRandomTermWithGivenVarDerivativeOrConstant(cls, vars=[], derivatives=[], constants=[], givenVar=None, givenDeriv=None,
                                                           givenConstant=None, max_power=3,
                                                           max_varsDerivsAndConstants=4):
        """
            Generates a random term with a supplied variable (class Variable), derivative (class Derivatif)
            or constant (class Constant), subject to the constraints that the maximum power associated
            with any of these is max_power and the maximum number of distinct variables, derivatives or
            constants is max_varsDerivsAndConstants. Exactly one of the parameters given_var,
            given_derivative and given_constant should be specified, others should have the value None.
            :param vars: List of variables (class Variable) to use
            :param derivatives: List of derivatives (class Derivatif) to use
            :param constants: List of constants (class Constant) to use
            :param givenVar: the given variable (class Variable) or None
            :param givenDeriv: the given derivative (class Derivatif) or None
            :param givenCnstant: the given constant (class Constant) or None
            :param max_power: a non-negative integer specifying the maximum power to raise any
                        variable, derivative or constant to in the term
            :param max_varsDerivsAndConstants: a non-negative integer specifying the maximum number
                        of distinct variables, derivatives and constants allowed in the term
            :return: the randomly generated term
        """

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
    def GenerateVDCSigDistributionDict(cls, vars, derivatives, constants, max_power=3, max_varsDerivsAndConstants=4):
        """
            Generate the variable-derivative-constant signature to frequency count dictionary. This dictionary
            is stored in EquationSystem._sigDict. It allows for selection of dimensionally consistent terms in a
            way that respects the probabilities established in this class' class-level variables. For a
            description of how the v-d-c signature is generated see the method Equation.GetVDCSignatureForTerm().
            :param vars: List of variables (class Variable) to use when creating the dictionary
            :param derivatives: List of derivatives (class Derivatif) to use when creating the dictionary
            :param constants: List of constants (class Constant) to use when creating the dictionary
            :param max_power: a non-negative integer specifying the maximum power to raise any
                            variable, derivative or constant to in any term
            :param max_varsDerivsAndConstants: a non-negative integer specifying the maximum number of distinct
                            variables, derivatives and constants allowed in any term
            :return: the dictionary
        """
        vdcSigDict = dict()
        for i in range(10000):
            term = Equation.GenerateRandomTerm(vars=vars, derivatives=derivatives,
                            constants=constants, max_power=max_power,
                            max_varsDerivsAndConstants=max_varsDerivsAndConstants) #seems like it is using 4 here but no_max later....
            sig = Equation.GetVDCSignatureForTerm(term)
            val = vdcSigDict.get(sig)
            if val is None:
                vdcSigDict[sig] = 1
            else:
                vdcSigDict[sig] = val + 1

        return vdcSigDict



    @classmethod
    def GetVDCSignatureForTerm(cls, term):
        """
            Obtains the variable-derivative-constant signature for a given term. It is easiest to
            understand how the signature works via a series of examples. Conisder the term c*x*y**2*z,
            where x, y, and z are variables and c is a constant.  No derivatives appear in the term.
            The signature will have segment designating the powers associated with the variables, then
            a segment designating the powers associated with the derivatives, followed by a segment
            designating the powers associated with the constants. The segments are separated by periods.
            In this case the signature is "112..1". The leading '112' means that there are three variables
            appearing in the term, two of which are raised to the 1st power and one of which is raised to
            the 2nd power. The '..' indicates that there are no derivatives in the term, and the final 1
            indicates that there is one constant in the term, and it is raised to the 1st power. The numbers
            each segment are sorted in increasing order, thus the initial '112' and not '121' or '211'.
            Another example is G**2*dxdt**2dydt, where G is a constant, and dxdt and dydt are derivatives.
            This term has signature ".12.2". J
            :param term: term to get the v-d-c signature for
            :return: the v-d-c signature as a string
        """
        var_powers = []
        deriv_powers = []
        const_powers = []
        if isinstance(term, Mul):
            for arg in term.args:
                if isinstance(arg, Constant):
                    const_powers.append(1)
                elif isinstance(arg, Variable):
                    var_powers.append(1)
                elif isinstance(arg, Derivatif):
                    deriv_powers.append(1)
                elif isinstance(arg, Pow):
                    if isinstance(arg.args[0], Constant):
                        const_powers.append(arg.args[1])
                    elif isinstance(arg.args[0], Variable):
                        var_powers.append(arg.args[1])
                    elif isinstance(arg.args[0], Derivatif):
                        deriv_powers.append(arg.args[1])
        elif isinstance(term, Pow):
            if isinstance(term.args[0], Constant):
                const_powers.append(term.args[1])
            elif isinstance(term.args[0], Variable):
                var_powers.append(term.args[1])
            elif isinstance(term.args[0], Derivatif):
                deriv_powers.append(term.args[1])
        elif isinstance(term, Constant):
            const_powers.append(1)
        elif isinstance(term, Variable):
            var_powers.append(1)
        elif isinstance(term, Derivatif):
            deriv_powers.append(1)

        var_powers.sort()
        deriv_powers.sort()
        const_powers.sort()
        var_sig = ""
        deriv_sig = ""
        const_sig = ""
        for p in var_powers:
            var_sig += str(p)
        for p in deriv_powers:
            deriv_sig += str(p)
        for p in const_powers:
            const_sig += str(p)

        return var_sig + '.' + deriv_sig + '.' + const_sig

    @classmethod
    def GenerateRandomTerm(cls, vars=[], derivatives=[], constants=[], max_power=3,
                                                           max_varsDerivsAndConstants=4):
        """
            Generate a random term given lists of variables, derivatives and constants, along with
            a maximum power any of these can be raised to, and a maximum number of distinct variables,
            derivatives or constants that can appear in the term
            :param vars: list of variables (class Variable) that can be used in the term
            :param derivatives: list of derivatives (class Derivatif) that can be used in the term
            :param constants: list of constants (class Constant) that can be used in the term
            :param max_power: a positive integer indicating the maximum power any variable, derivative
                            or constant can be raised to in the term
            :param max_varsDerivsAndConstants: a positive integer indicating the maximum number of
                            distinct variables, derivatives or constants that can appear in the term
            :return: the randomly generated term
        """
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
        """
            Generates a random unnamed constant (always a small integer > 0) for a term given the v
            alues specified in  Equation.PROB_TERM_HAS_NON_UNITAL_INTEGER_CONSTANT and
            quation.PROB_OF_SMALL_INTEGER_CONSTANTS
            :return: randomly generated small positive integer
        """
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
        """
            Returns a set consisting of variables, derivatives and constants that are common
            to all the given terms
            :param terms: list of terms
            :return: set of common factors
        """
        common_factors = terms[0].free_symbols
        for i in range(1, len(terms)):
            factors = terms[i].free_symbols
            common_factors = common_factors.intersection(factors)

        return common_factors


    @classmethod
    def TermAmongExistingTerms(cls, existingTerms, term):
        """
        Determines with a given term is contained in a list of existing terms
        :param existingTerms: list of existing terms
        :param term: term to test for
        :return: True if the term is contained in existingTerms, False otherwise
        """
        for t in existingTerms:
            if t == term:
                return True

        return False

    @classmethod
    def AssignRandomSignsToTerms(cls, terms):
        """
            Assigns random signs to a set of terms (assumed to be all the terms of an Equation).
            With a small probability (= Equation.HOMOGENEOUS_EQN_PROB) it will assign all positive signs
            to the terms. Otherwise, it will guarantee that there is at least one positive and one negative
            term.
            :param terms: List of terms to assign signs to
            :return: None. The changes are done in place to the list of terms
        """
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

    def sanityCheck(self):
        """
            Performs a variety of different sanity checks on an Equation. For details see
            Equation.SanityCheckFromTerms(). Returns True of the Equation passes the sanity checks
            and False otherwise.
            :return: A Boolean indicating whether or not the Equation has passed the sanity check(s)
        """
        terms = self.getTerms()
        return Equation.SanityCheckFromTerms(terms)

    @classmethod
    def SanityCheckFromTerms(cls, terms):
        """
            Performs a variety of different sanity checks on a set of terms (assumed to be the terms
            for a candidate Equation). Returns True if the terms pass the sanity checks and False
            otherwise. There are currently two sanity checks (though we expect others to be added
            over time): (i) a check that there are not just two terms, each containing a single
            variable, derivative or constant. Such an equation would imply that one of the two
            variables, derivatives or constants could be removed, and (ii) a check that there are
            not two identical terms modulo unnamed constants.
            :param terms: A list of terms to be sanity checked
            :return: A Boolean indicating whether the terms have passed the sanity check(s)
        """
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
        """
            Returns a tuple corresponding to the terms in the Equation
            :return: a tuple of terms
        """
        return self._poly.expr.args

    def __eq__(self, eqn):
        return self.getTerms() == eqn.getTerms()

    def __hash__(self):
        return hash(str(self))


    def equalModUnnamedConstants(self, eqn):
        """
            Determine whether the current equation is equal a given equation modulo unnamed
            constants. In other words, do the two equations have the same terms with their
            unnamed constants stripped off?
            :param eqn: Equation to compare this one to
            :return: Boolean indicating whether the two equations are equal or not
        """
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
        """
            Determine whether two given terms are equal modulo their unnamed constants
            :param term1: First term to compare
            :param term2: Second term to compare
            :return: True of the terms are equal modulo their unnamed constnats, False otherwise
        """
        t1 = Equation.GetUnnamedConstantStrippedTerm(term1)
        t2 = Equation.GetUnnamedConstantStrippedTerm(term2)

        return t1 == t2

    @classmethod
    def TermsEqual(self, term1, term2):
        """
            Determine whether two given terms are equal
            :param term1: First term to compare
            :param term2: Second term to compare
            :return: True if equal, False otherwise
        """
        return term1 == term2

    @classmethod
    def GetUnnamedConstantStrippedTerm(cls, term):
        """
            Given a term, return the term stripped of any unnamed constant it may have
            :param term: Term to be stripped
            :return: Result of stripping the unnamed constant from the given term
        """
        t = term
        if len(t.args) > 0 and isinstance(t.args[0], Integer):
            return t / t.args[0]
        else:
            return t

    @classmethod
    def GetUnnamedConstantForTerm(cls, term):
        """
            Get the unnamed constant associated with a provided term
            :param term: Provided term
            :return: Unnamed constant. Returns 1 if there is non unnamed constant.
        """
        if len(term.args) > 0 and isinstance(term.args[0], Integer):
            return term.args[0]
        else:
            return Integer(1)

    def divideByCommonUnnamedConstants(self):
        """
            Divide the current equation by any common unnamed constant all of its terms may have
            :return: None
        """
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
        """
            Setup runtime logging.
            :param filename: File to contain the output log
            :param filemode: Should generally be set to the default 'w' for 'write'
            :param encoding: Should generally be set to the default 'utf-8' though other encodings are
                        possible
            :param level: Possible levels are logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, and
                        logging.CRITICAL. Indicates the lowest level of debug message that will be output
                        the file specified by filename.
            :return: None
        """
        logging.basicConfig(filename=filename, filemode=filemode, encoding=encoding, level=level)

    def __str__(self):
        return str(self._poly.expr) + " (U of M: " + str(self.getUofM()) + ")"



















