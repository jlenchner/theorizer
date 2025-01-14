from sympy import *
from variable import *
from derivative import *
import random
import math


#This is a wrapper around the sympy Poly class with some extra functionality
class Equation:
    PROB_OF_NUM_TERMS = [0, 1.0/3.0, 1.0/2.0, 1.0/6.0]  #probs of a randomly generated poly having 1 term, 2 terms, etc.
    PROB_FACTORS_PER_TERM = [0.4, 0.4, 0.2]  #probs that a term in a randomly generated poly will have 1, 2, 3, etc. factors
    SAME_FACTOR_BIAS = 2.0  #bias for repeating the same factor in a term (in other words, favoring x^2 over xy)

    PROB_TERM_HAS_NON_UNITAL_CONSTANT = 0.2
    PROB_OF_SMALL_INTEGER_CONSTANTS = [0.0, 0.0, 0.8, 0.1, 0.05, 0.05]  # prob that a small non-unital integer constant is [0,1,2,3,4,5]

    NO_MAX = -1

    def __init__(self, exp):
        self._poly = Poly(exp)  #The equation is treated generically like a polynomial of symbols
        self._variables, self._derivatives = Equation.InferVarsAndDerivativesFromExpression(exp)


    def isDimensionallyConsistent(self):
        #should first check that all the variable, derivatives and constants have associated units of measure
        #Also should print warnings about entitites without units
        for var in self._variables:
            if var._u_of_m is None:
                print("The variable " + str(var) + " has no associated units of measure!")
                return False
        for derivative in self._derivatives:
            if derivative._u_of_m is None:
                print("The derivative " + str(derivative) + " has no associated units of measure!")
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



    def add(self, exp, vars = []):
        self._poly = self._poly.add(Poly(exp))
        if len(vars) == 0:
            vars = Equation.InferVarsFromExpression(exp)
        for var in vars:
            if var not in self._variables:
                self._variables.append(var)



    @classmethod
    def InferVarsAndDerivativesFromExpression(cls, exp):   #this will not retain the units of measure!
        vars = []
        derivtives = []
        for sym in exp.free_symbols:
            if isinstance(sym, Variable):
                vars.append(sym)
            elif isinstance(sym, Derivatif):
                derivtives.append(sym)

        return vars, derivtives





    @classmethod
    def GenerateRandom(cls, vars, derivatives=[], max_var_and_derivatives_per_eqn=NO_MAX):  #Does not yet deal with picking coefficients (constants)
                                                    #Does not yet deail with generation using derivatives
        if max_var_and_derivatives_per_eqn == Equation.NO_MAX:
            max_var_and_derivatives_per_eqn = len(vars) + len(derivatives)

        r = random.random()
        num_terms = len(Equation.PROB_OF_NUM_TERMS)
        prob_so_far = 0.0
        for i in range(len(Equation.PROB_OF_NUM_TERMS)):
            prob_so_far += Equation.PROB_OF_NUM_TERMS[i]
            if r < prob_so_far:
                num_terms = i+1
                break

        terms = None
        while True:
            firstTerm = Equation.GenerateRandomTerm(vars, derivatives)
            terms = [firstTerm]
            vars_and_derivs_in_use = firstTerm.free_symbols
            for i in range(1, num_terms):
                while True:
                    term = Equation.GenerateRandomTerm(vars, derivatives)
                    vars_and_derivs_in_term = term.free_symbols
                    new_vars_and_derivs = vars_and_derivs_in_term - vars_and_derivs_in_use
                    if len(new_vars_and_derivs) + len(vars_and_derivs_in_use) <= max_var_and_derivatives_per_eqn \
                            and not Equation.TermAmongExistingTerms(terms, term):
                        vars_and_derivs_in_use = vars_and_derivs_in_use.union(new_vars_and_derivs)
                        terms.append(term)
                        break
            if Equation.GetCommonFactors(terms) == set():
                break
            #else:
            #    print("Common factors found among all terms! Regenerating random poly!")

        Equation.AssignRandomSignsToTerms(terms)

        exp = terms[0]
        for i in range(1, len(terms)):
            exp = exp + terms[i]

        return Equation(exp)

    @classmethod
    def GetAllDimensionallyConsistentTerms(cls, term_to_match, vars=[], derivatives=[], max_power=3, max_varsAndDerivs=-1):
        if max_varsAndDerivs == -1:
            max_varsAndDerivs = len(vars) + len(derivatives)

        allConsistentTerms = []
        u_of_m_to_match = Equation.GetUofMForTerm(term_to_match)
        for i in range(1, pow(max_power + 1, len(vars) + len(derivatives))):
            baseNum = Equation.ToBase(i, max_power + 1)
            s_baseNum = str(baseNum)
            if len(s_baseNum) < len(vars) + len(derivatives):
                s_baseNum = '0' * (len(vars) + len(derivatives) - len(s_baseNum)) + s_baseNum
            numZeroes = s_baseNum.count('0')
            numVarsAndDerivs = len(vars) + len(derivatives) - numZeroes
            if numVarsAndDerivs <= max_varsAndDerivs:
                term = Equation.GenerateTermFromBaseNum(baseNum, vars, derivatives)
                u_of_m = Equation.GetUofMForTerm(term)
                if u_of_m_to_match == u_of_m and not Equation.TermsEqualModConstants(term, term_to_match):
                    allConsistentTerms.append(term)

        return allConsistentTerms

    @classmethod
    def ToBase(cls, num, base):
        res = 0
        max_pow = ceiling(math.log(num, base))
        for i in range(max_pow+1):
            res += pow(10,i)*(num % base)
            num -= (num % base)
            num = int(num/base)

        return int(res)

    @classmethod
    def GenerateTermFromBaseNum(cls, baseNum, vars=[], derivatives=[]):
        term = None

        varsAndDerivs = []
        varsAndDerivs.extend(vars)
        varsAndDerivs.extend(derivatives)

        s_baseNum = str(baseNum)
        if  len(s_baseNum) < len(varsAndDerivs):
            s_baseNum = '0' * (len(varsAndDerivs) - len(s_baseNum)) + s_baseNum

        for i in range(len(s_baseNum)):
            val = int(s_baseNum[i])
            if val > 0:
                varOrDeriv = varsAndDerivs[i]
                if term is not None:
                    if val == 1:
                        term = term*varOrDeriv
                    else:
                        term = term * (varOrDeriv**val)
                else:
                    if val == 1:
                        term = varOrDeriv
                    else:
                        term = varOrDeriv**val

        return term

    @classmethod
    def GenerateRandomDimensionallyConsistent(cls, vars, derivatives, max_var_and_derivatives_per_eqn=NO_MAX):  # Does not yet deal with picking coefficients
        if max_var_and_derivatives_per_eqn == Equation.NO_MAX:
            max_var_and_derivatives_per_eqn = len(vars) + len(derivatives)

        #Needs to be modified to generate the first random term and then use getAllDimensionallyConsistentTerms()

        r = random.random()
        num_terms = len(Equation.PROB_OF_NUM_TERMS)
        prob_so_far = 0.0
        for i in range(len(Equation.PROB_OF_NUM_TERMS)):
            prob_so_far += Equation.PROB_OF_NUM_TERMS[i]
            if r < prob_so_far:
                num_terms = i + 1
                break

        terms = None
        while True:
            firstTerm = Equation.GenerateRandomTerm(vars, derivatives)
            terms = [firstTerm]
            vars_and_derivs_in_use = firstTerm.free_symbols
            candidateTerms  = Equation.GetAllDimensionallyConsistentTerms(term_to_match=firstTerm, vars=vars, derivatives=derivatives)
            if len(candidateTerms) >= num_terms - 1:
                for i in range(num_terms - 1):
                    while True:
                        rand_int = random.randint(0, len(candidateTerms) - 1)
                        term = candidateTerms[rand_int]
                        vars_and_derivs_in_term = term.free_symbols
                        new_vars_and_derivs = vars_and_derivs_in_term - vars_and_derivs_in_use
                        if len(new_vars_and_derivs) + len(vars_and_derivs_in_use) <= max_var_and_derivatives_per_eqn:
                            vars_and_derivs_in_use = vars_and_derivs_in_use.union(new_vars_and_derivs)
                            terms.append(term)
                            candidateTerms.remove(term)
                            break

            if Equation.GetCommonFactors(terms) == set():
                break

        Equation.AssignRandomSignsToTerms(terms)

        exp = terms[0]
        for i in range(1, len(terms)):
            exp = exp + terms[i]

        return Equation(exp)

    def getSymbolsUsed(self):  #returns a set of variables
        return self._poly.free_symbols

    def getUofM(self, sym):
        for var in self._variables:
            if var.name == str(sym):
                return var._u_of_m

        return None


    @classmethod
    def GenerateRandomTerm(cls, vars, derivatives=[]):
        vars_and_derivs = []
        vars_and_derivs.extend(vars)
        vars_and_derivs.extend(derivatives)

        num_factors = len(Equation.PROB_FACTORS_PER_TERM)
        prob_so_far = 0.0
        r = random.random()
        for i in range(len(Equation.PROB_FACTORS_PER_TERM)):
            prob_so_far += Equation.PROB_FACTORS_PER_TERM[i]
            if r < prob_so_far:
                num_factors = i+1
                break

        vars_and_derivs_so_far = []  #An array indicating counts of the vars and derivatives that are already being used
        var_and_deriv_probs = []
        for i in range(len(vars_and_derivs)):
            vars_and_derivs_so_far.append(0)
            var_and_deriv_probs.append(1.0/len(vars_and_derivs))

        for i in range(num_factors):
            r = random.random()
            var_or_deriv = var_and_deriv_probs[len(var_and_deriv_probs) - 1]
            prob_so_far = 0.0
            for j in range(len(vars_and_derivs)):
                prob_so_far += var_and_deriv_probs[j]
                if r < prob_so_far:
                    var_or_deriv = vars_and_derivs[j]
                    bonus_prob_factpr = 1.0
                    if vars_and_derivs_so_far[j] == 0:
                        old_prob = var_and_deriv_probs[j]
                        var_and_deriv_probs[j] *= Equation.SAME_FACTOR_BIAS
                        aggregate_weight = 1 + var_and_deriv_probs[j] - old_prob
                        # must now rebalance the probs
                        for k in range(len(vars_and_derivs)):
                            var_and_deriv_probs[k] = var_and_deriv_probs[k]/aggregate_weight
                    vars_and_derivs_so_far[j] += 1
                    break

        term = None
        for i in range(len(vars_and_derivs)):
            if vars_and_derivs_so_far[i] == 1:
                if term is not None:
                    term = term*vars_and_derivs[i]
                else:
                    term = vars_and_derivs[i]
            elif vars_and_derivs_so_far[i] > 1:
                if term is not None:
                    term = term*vars_and_derivs[i]**vars_and_derivs_so_far[i]
                else:
                    term = vars_and_derivs[i]**vars_and_derivs_so_far[i]

        r = random.random()
        if r < Equation.PROB_TERM_HAS_NON_UNITAL_CONSTANT:
            c = Equation.GetConstant()
            term = c*term

        return term


    @classmethod
    def GetConstant(cls):  #should also test for common constant factors now
        c = 1
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

    def getTerms(self):
        return self._poly.expr.args

    @classmethod
    def TermsEqualModConstants(cls, term1, term2):
        t1 = term1
        t2 = term2
        if len(t1.args) > 0 and isinstance(t1.args[0], Integer):
            t1 = t1 / t1.args[0]
        if len(t2.args) > 0 and isinstance(t2.args[0], Integer):
            t2 = t2 / t2.args[0]

        return t1 == t2


    def __str__(self):
        return str(self._poly.expr)








