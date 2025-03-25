"""
    This is the driver file for the various experiments. Use and freely modify this file,
    or, if preferred, create your own driver file and run that.
"""
import random

from equationSystem import *
from equation import *
from m2_functions import *
from unitOfMeasure import *
from derivative import *
from constant import *


def GeneratePseudoTheoriesWithReplacement(vars=[], derivatives=[], constants=[],
                                          measuredVars=[], minEqns=4, maxEqns=4,
                                          max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX,
                                          numTheories=10, numReplacements=5):
    """
        Helper method for generating random EquationSystems with specified variables, derivatives
        and constants, but also with a randomly selected number of equations, ranging from 4 to 8.
        :param vars: List of variables (class Variable) to use
        :param derivatives: List of derivatives (class Derivatif) to use
        :param constants: List of constants (class Constant) to use
        :param measuredVars: List of the variables that are measured (not currently used)
        :param minEqns: The minimum number of equations to generate
        :param maxEqns: The maximum number of equations to generate
        :param max_vars_derivatives_and_constants_per_eqn: The maximum number of distinct variables,
                    derivatives and constants that are allowed to appear in any one Equation.
        :param numTheories: Number of EuqationSystems to generate
        :param numReplacements: Number of random replacements to generate per generated EquationSystem
        :return: None
    """
    for i in range(numTheories):
        numEqns = random.randint(minEqns, maxEqns)
        eqnSystem = EquationSystem.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivatives,
                                                                         constants=constants,
                                                                         measuredVars=measuredVars, numEqns=numEqns,
                                                                         max_vars_derivatives_and_constants_per_eqn=max_vars_derivatives_and_constants_per_eqn)
        print("\n\nSystem " + str(i + 1) + ": " + str(eqnSystem) + "\n\n")
        print("\n\nReplacements: \n")
        if numReplacements > 0:
            for j in range(numReplacements):
                index, eqn = eqnSystem.copy().replaceRamdomDimensionallyConsistentEqn()
                print("Replacement eqn index: " + str(index))
                print("Replacement eqn: " + str(eqn))


def KeplerReplacement():
    """
        Starting with the Equations for discovering Kepler's 3rd Law, derive a new
        dimensionally consistent EquationSystem by replacing each of the equations
        in turn 10 times by a random dimensionally consistent equation. The replacement
        equations are printed out.
        :return: None
    """
    d1, d2, m1, m2, Fg = variables('d1, d2, m1, m2, Fg')
    w = Variable('w', UofM(1 / UofM.s))
    vars = [d1, d2, m1, m2, w, Fg]
    G = Constant('G')
    constants = [G]
    eqn0 = Equation(d1 * m1 - d2 * m2)
    eqn1 = Equation(Fg * d1 ** 2 + 2 * Fg * d1 * d2 + Fg * d2 ** 2 - G * m1 * m2)
    eqn2 = Equation(Fg - m2 * d2 * w ** 2)
    eqns = [eqn0, eqn1, eqn2]
    eqnSys = EquationSystem(vars=vars, derivatives=[], constants=constants, equations=eqns,
                            max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
    print("Kepler System Equation Replacements:\n")
    for i in range(10):
        newEqn = eqnSys.replaceRandomDimensionallyConsistentEqnByIndex(eqnIndex=0)
        print("New eqn for eqn0: " + str(newEqn))
    print("\n")
    for i in range(10):
        newEqn = eqnSys.replaceRandomDimensionallyConsistentEqnByIndex(eqnIndex=1)
        print("New eqn for eqn1: " + str(newEqn))
    print("\n")
    for i in range(10):
        newEqn = eqnSys.replaceRandomDimensionallyConsistentEqnByIndex(eqnIndex=2)
        print("New eqn for eqn2: " + str(newEqn))


def KeplerikeTheoriesWithReplacement():
    """
        Generate random dimensionally consistent EquationSystems using the Kepler variables
        and constants but none of the actual equations. Each of these "Kepler-like theories"
        should have a random number of equations between 4 and 8. The EquationSystems are
        printed out.
        :return: None
    """
    d1, d2, m1, m2, Fg = variables('d1, d2, m1, m2, Fg')
    w = Variable('w', UofM(1 / UofM.s))
    vars = [d1, d2, m1, m2, w, Fg]
    G = Constant('G')
    constants = [G]

    GeneratePseudoTheoriesWithReplacement(vars=vars, derivatives=[], constants=constants,
                                          measuredVars=vars, minEqns=4, maxEqns=8,
                                          max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX,
                                          numTheories=20, numReplacements=0)


def KeplerikeTheoriesWithDerivsAndhReplacement():
    """
        Generate random dimensionally consistent EquationSystems using the Kepler variables
        and constants plus an assortment of derivatives (but none of the actual equations).
        Each of these "Kepler-like theories with derivatives" again has a random number of
        equations between 4 and 8. The EquationSystems are printed out.
        :return: None
    """
    d1, d2, m1, m2, Fg = variables('d1, d2, m1, m2, Fg')
    w = Variable('w', UofM(1 / UofM.s))
    vars = [d1, d2, m1, m2, w, Fg]
    dx1dt, d2x1dt2, dx2dt, d2x2dt2 = derivatives('dx1dt,d2x1dt2,dx2dt,d2x2dt2')
    derivs = [dx1dt, d2x1dt2, dx2dt, d2x2dt2]
    G = Constant('G')
    constants = [G]

    GeneratePseudoTheoriesWithReplacement(vars=vars, derivatives=derivs, constants=constants,
                                          measuredVars=vars, minEqns=4, maxEqns=8,
                                          max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX,
                                          numTheories=200, numReplacements=10)


def TimeDilationReplacement():
    """
        Starting with the Equations for discovering the law of relativistic time dilation,
        derive a new dimensionally consistent EquationSystem by replacing each of the equations
        in turn 10 times by a random dimensionally consistent equation. The replacement
        equations are printed out.
        :return: None
    """
    dt = Variable('dt', UofM(UofM.s))
    dt0 = Variable('dt0', UofM(UofM.s))
    d, v = variables('d,v')
    L = Variable('L', UofM(UofM.m))
    f = Variable('f', UofM(1 / UofM.s))
    f0 = Variable('f0', UofM(1 / UofM.s))
    c = Constant('c')
    # vars = [dt, dt0, d, v, L, f, f0]
    vars = [dt, dt0, d, v, L]
    constants = [c]
    eqn0 = Equation(c * dt0 - 2 * d)
    eqn1 = Equation(c * dt - 2 * L)
    eqn2 = Equation(4 * L ** 2 - 4 * d ** 2 - v ** 2 * dt ** 2)
    eqns = [eqn0, eqn1, eqn2]
    eqnSys = EquationSystem(vars=vars, derivatives=[], constants=constants, equations=eqns,
                            max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
    print("Time Dilation System Equation Replacements:\n")
    for i in range(10):
        newEqn = eqnSys.replaceRandomDimensionallyConsistentEqnByIndex(eqnIndex=0)
        print("New eqn for eqn0: " + str(newEqn))
    print("\n")
    for i in range(10):
        newEqn = eqnSys.replaceRandomDimensionallyConsistentEqnByIndex(eqnIndex=1)
        print("New eqn for eqn1: " + str(newEqn))
    print("\n")
    for i in range(10):
        newEqn = eqnSys.replaceRandomDimensionallyConsistentEqnByIndex(eqnIndex=2)
        print("New eqn for eqn2: " + str(newEqn))


def TimeDilationLikeTheoriesWithReplacement():
    """
       Generate random dimensionally consistent EquationSystems using the variables and
       constants for discovering the law of relativistic time dilation but none of the
       actual equations. Each of these "Kepler-like theories" should have
       a random number of equations between 4 and 8. The EquationSystems are printed out.
       :return: None
   """
    dt = Variable('dt', UofM(UofM.s))
    dt0 = Variable('dt0', UofM(UofM.s))
    d, v = variables('d,v')
    L = Variable('L', UofM(UofM.m))
    f = Variable('f', UofM(1 / UofM.s))
    f0 = Variable('f0', UofM(1 / UofM.s))
    c = Constant('c')
    vars = [dt, dt0, d, v, L, f, f0]
    constants = [c]

    GeneratePseudoTheoriesWithReplacement(vars=vars, derivatives=[], constants=constants,
                                          measuredVars=vars, minEqns=4, maxEqns=4,
                                          max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)


if __name__ == "__main__":
    Equation.SetLogging()

    N = UofM(UofM.kg * UofM.m / UofM.s)
    F = Variable('F', N)

    one = UofM(1)
    pi = Constant('pi')

    # KeplerReplacement()
    # TimeDilationReplacement()
    # TimeDilationLikeTheoriesWithReplacement()
    # KeplerikeTheoriesWithReplacement()
    # KeplerikeTheoriesWithDerivsAndhReplacement()

    # sxit(0)
    x, y, z, w = variables('x,y,z,w')
    eqn = Equation(2 * x - 2 * y)
    terms = eqn.getTerms()
    for term in terms:
        c = Equation.GetUnnamedConstantForTerm(term)
    eqn.divideByCommonUnnamedConstants()

    terms = [x * y ** 3 * z, w * y ** 2, y * z, 3 * y * z * x]
    for term in terms:
        sig = Equation.GetVDCSignatureForTerm(term)
        print(str(term) + ': ' + sig)
    res = Equation.SanityCheckFromTerms(terms)
    terms = [x * y * z, w * y, y * z, x]
    res = Equation.SanityCheckFromTerms(terms)

    eqn1 = Equation(2 * x + 3 * y)
    eqn2 = Equation(3 * y + 2 * x)

    F = Variable('F')

    d1, d2, m1, m2, W, p, Fc, Fg, T, E = variables('d1,d2,m1,m2,W,p,Fc,Fg,T,E')
    dxdt, d2xdt2 = derivatives('dxdt,d2xdt2')
    dx1dt, d2x1dt2, dx2dt, d2x2dt2 = derivatives('dx1dt,d2x1dt2,dx2dt,d2x2dt2')
    G, c, pi = constants('G,c,pi')
    # vars = [Fc,Fg,W,d1,d2,m1,m2,p,T,E]
    vars = [Fc, Fg, W, d1, d2, m1, m2, p]
    # derivs = [dxdt, d2xdt2]
    derivs = [dx1dt, d2x1dt2, dx2dt, d2x2dt2]
    # constants = [G, c, pi]
    constants = [G, c]

    # sigDict = Equation.GenerateVDCSigDistributionDict(vars=vars, derivatives=derivs, constants=constants, max_power=3, max_varsDerivsAndConstants=4)
    # print("got here")

    v = Variable('v')
    v0 = Variable('v', UofM(UofM.m / UofM.s))

    res = pi.isDimensionless()
    term = Equation.GenerateTermFromBaseNum(31200113, vars, [])

    eqn = Equation(c * d2x1dt2 * m2 - c * Fg - d2x1dt2 * dx1dt * m1 + d2x1dt2 * dx2dt * m2)
    terms = eqn.getTerms()
    commonFactors = Equation.GetCommonFactors(terms)

    term1 = 3 * term
    term2 = term * 5
    term3 = d1
    terms = [term, term1, term2, term3]
    cf = Equation.GetCommonFactors(terms)
    term1 = 3 * W * p
    term2 = p * W * 3
    eq1 = Equation(2 * m1 * c ** 2 - d2xdt2 * G)
    eq2 = Equation(-d2xdt2 * G + m1 * 2 * c ** 2)
    res = eq1 == eq2
    eq2 = Equation(-d2xdt2 * 4 * G + m1 * 2 * c ** 2)
    res = eq1.equalModUnnamedConstants(eq2)
    res = Equation.TermsEqualModUnnamedConstants(term1, term2)
    res = Equation.TermsEqual(term1, term2)

    # vars.extend(constants)
    eqns = []
    eqns.append(Equation(d2xdt2 * d2 ** 2 * m2 - 2 * W * d1))
    eqns.append(Equation(G * Fc * d1 + G * W - 2 * d2xdt2 ** 2 * d1 ** 2 * d2))
    eqns.append(Equation(c ** 2 * d1 * m1 - Fg * d2 ** 2))
    eqns.append(Equation(-c * p + Fg * d2 + W))
    EquationSystem.SanityCheckEquationList(eqns)

    # print("Non-Dimensionally Consistent Set:\n\n")
    # for i in range(200):
    #    eqnSystem = EquationSystem.GenerateRandom(vars=vars, derivatives=derivs, measuredVars=vars, \
    #                                          constants=constants, numEqns=4, max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
    #    print(str(i + 1) + ": " + str(eqnSystem) + "\n\n")
    # res = eqnSystem._equations[0].isDimensionallyConsistent()

    # for i in range(100):
    #    eqn = Equation.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=[])
    #    print(str(i+1) + ": " + str(eqn))
    print("Dimensionally Consistent Set:\n\n")
    for i in range(200):
        eqnSystem = EquationSystem.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivs,
                                                                         constants=constants,
                                                                         measuredVars=vars, numEqns=4,
                                                                         max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
        print(str(i + 1) + ": " + str(eqnSystem) + "\n\n")
        # eqn = eqnSystem.replaceRandomDimensionallyConsistentEqnByIndex(1)

    # measured_vars = [d1,d2,m1,m2,p]
    eq2 = Equation(d1 ** 2 * Fg + 2 * d1 * d2 * Fg + d2 ** 2 * Fg - m1 * m2)
    res = eq2.isDimensionallyConsistent()
    eq2prime = Equation(d1 ** 2 * Fg + 2 * d1 * d2 * Fg + d2 ** 2 * Fg)
    res = eq2prime.isDimensionallyConsistent()
    vbles = eq2._variables
    terms = eq2.getTerms()
    m, s = base_units('m,s')  # see below
    units = UofM(
        m / s)  # should not be done in terms of symbols - perhaps a new class called units, which would be veneer over Symbols

    v = Variable('v')
    dxdt = Derivatif('dxdt')
    eq = Equation(v - dxdt)
    terms2 = eq.getTerms()
    res = eq.isDimensionallyConsistent()

    F, m, a = variables('F,m,a')
    dxdt, d2xdt2 = derivatives('dxdt,d2xdt2')
    eq3 = Equation(F - m * d2xdt2)
    print(eq3)
    res = eq3.isDimensionallyConsistent()
    vars = [F, m, a]
    derivatives = [dxdt, d2xdt2]

    eqn = Equation.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivatives, constants=[])
    print(str(eqn))
    res = eqn.isDimensionallyConsistent()

    # eq3 = Polynomial(Fc - m2*d2*w**2)
    # eq4 = Polynomial(Fg - Fc)
    # eq5 = Polynomial(w*p - 1)
    # eqns = [eq2, eq3, eq4, eq5]
    # eqnSystem = EquationSystem(vars=vars, measuredVars=measured_vars, equations=eqns)

    # print("System: \n" + eqnSystem.toString() + "\n")

    # isConsistent(eqnSystem)
    # print(checkConsistencyReplacedSystems(eqnSystem, 3, 10))
    # print(checkConsistencyRandomSystems(vars, 10))
    # project(eqnSystem)

    vars = [Fc, Fg, W, d1, d2, m1, m2, p]
    # Take the measured vars to be a random subset of the variables
    # measured_vars = random.sample(vars, len(vars)//2)

    measured_vars = vars.copy()
    random.shuffle(measured_vars)

    EquationSystem.ProjectRandomSystems(vars, derivatives, measured_vars, 10)




