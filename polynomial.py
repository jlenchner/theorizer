from sympy import *
import random


#This is a wrapper around the sympy Poly class with some extra functionality
class Polynomial:
    PROB_OF_NUM_TERMS = [0, 1.0/3.0, 1.0/2.0, 1.0/6.0]  #probs of a randomly generated poly having 1 term, 2 terms, etc.
    PROB_FACTORS_PER_TERM = [0.4, 0.4, 0.2]  #probs that a term in a randomly generated poly will have 1, 2, 3, etc. factors
    SAME_FACTOR_BIAS = 2.0  #bias for repeating the same factor in a term (in other words, favoring x^2 over xy)

    PROB_TERM_HAS_NON_UNITAL_CONSTANT = 0.2
    PROB_OF_SMALL_INTEGER_CONSTANTS = [0.0, 0.0, 0.8, 0.1, 0.05, 0.05]  # prob that a small non-unital integer constant is [0,1,2,3,4,5]

    NO_MAX = -1

    def __init__(self, exp):
        self._poly = Poly(exp)

    def add(self, exp):
        self._poly = self._poly.add(Poly(exp))

    @classmethod
    def GenerateRandom(cls, vars, max_vars=NO_MAX):  #Does not yet deal with picking coefficients
        if max_vars == Polynomial.NO_MAX:
            max_vars = len(vars)

        r = random.random()
        num_terms = len(Polynomial.PROB_OF_NUM_TERMS)
        prob_so_far = 0.0
        for i in range(len(Polynomial.PROB_OF_NUM_TERMS)):
            prob_so_far += Polynomial.PROB_OF_NUM_TERMS[i]
            if r < prob_so_far:
                num_terms = i+1
                break

        terms = None
        while True:
            terms = []
            vars_in_use = set()
            for i in range(num_terms):
                while True:
                    term = Polynomial.GenerateRandomTerm(vars)
                    vars_in_term = term.free_symbols
                    new_vars = vars_in_term - vars_in_use
                    if len(new_vars) + len(vars_in_use) <= max_vars and not Polynomial.TermAmongExistingTerms(terms, term):
                        vars_in_use = vars_in_use.union(new_vars)
                        terms.append(term)
                        break
            if Polynomial.GetCommonFactors(terms) == set():
                break
            #else:
            #    print("Common factors found among all terms! Regenerating random poly!")

        Polynomial.AssignRandomSignsToTerms(terms)

        exp = terms[0]
        for i in range(1, len(terms)):
            exp = exp + terms[i]

        return Polynomial(exp)

    def getVarsUsed(self):  #returns a set of variables
        return self._poly.free_symbols


    @classmethod
    def GenerateRandomTerm(cls, vars):
        num_factors = len(Polynomial.PROB_FACTORS_PER_TERM)
        prob_so_far = 0.0
        r = random.random()
        for i in range(len(Polynomial.PROB_FACTORS_PER_TERM)):
            prob_so_far += Polynomial.PROB_FACTORS_PER_TERM[i]
            if r < prob_so_far:
                num_factors = i+1
                break

        vars_so_far = []  #An array indicating counts of the vars that are already being used
        var_probs = []
        for i in range(len(vars)):
            vars_so_far.append(0)
            var_probs.append(1.0/len(vars))

        for i in range(num_factors):
            r = random.random()
            var = vars[len(vars) - 1]
            prob_so_far = 0.0
            for j in range(len(vars)):
                prob_so_far += var_probs[j]
                if r < prob_so_far:
                    var = vars[j]
                    bonus_prob_factpr = 1.0
                    if vars_so_far[j] == 0:
                        old_prob = var_probs[j]
                        var_probs[j] *= Polynomial.SAME_FACTOR_BIAS
                        aggregate_weight = 1 + var_probs[j] - old_prob
                        # must now rebalance the probs
                        for k in range(len(vars)):
                            var_probs[k] = var_probs[k]/aggregate_weight
                    vars_so_far[j] += 1
                    break

        term = None
        for i in range(len(vars)):
            if vars_so_far[i] == 1:
                if term is not None:
                    term = term*vars[i]
                else:
                    term = vars[i]
            elif vars_so_far[i] > 1:
                if term is not None:
                    term = term*vars[i]**vars_so_far[i]
                else:
                    term = vars[i]**vars_so_far[i]

        r = random.random()
        if r < Polynomial.PROB_TERM_HAS_NON_UNITAL_CONSTANT:
            c = Polynomial.GetConstant()
            term = c*term

        return term


    @classmethod
    def GetConstant(cls):  #should also test for common constant factors now
        c = 1
        r = random.random()
        prob_so_far = 0.0
        for i in range(len(Polynomial.PROB_OF_SMALL_INTEGER_CONSTANTS)):
            prob_so_far += Polynomial.PROB_OF_SMALL_INTEGER_CONSTANTS[i]
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

    def toString(self):
        return str(self._poly.expr)
