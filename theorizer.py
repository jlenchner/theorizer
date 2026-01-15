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


def NewtonianGravity():
    d1, d2, d, m1, m2, Fg, Fc = variables('d1, d2, d, m1, m2, Fg, Fc')
    w = Variable('w', UofM(1 / UofM.s))
    p = Variable('p', UofM(UofM.s))
    vars = [d1, d2, d, m1, m2, w, p, Fg, Fc]
    G, pi = constants('G, pi')
    cons = [G, pi]
    eqn0 = Equation(d1 * m1 - d2 * m2)
    eqn1 = Equation(d - d1 - d2)
    eqn2 = Equation(Fg * d1**2 -  G * m1 * m2)
    eqn3 = Equation(Fc - m2 - d2 - w**2)
    eqn4 = Equation(Fc - Fg)
    eqn5 = Equation(p * w - 2 * pi)
    eqns = [eqn0, eqn1, eqn2, eqn3, eqn4, eqn5]
    eqnSys = EquationSystem(vars=vars, derivatives=[], constants=cons, measuredVars=vars, equations=eqns,
                                max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
    #for i in range(10):
    #    index, newEqn = eqnSys.replaceRamdomDimensionallyConsistentEqn()
    #    print("New eqn for index = " + str(index) + ": " + str(newEqn))
    #print("\n")
    return eqnSys

def SpecialRelativity():
    dt0, dt, d, v, f0, f = variables('dt0, dt, d, v, f0, f')
    L = Variable('L', UofM(UofM.m))
    vars = [dt0, dt, d, v, f0, f, L]
    c, one = constants('c, one')
    cons = [c, one]
    eqn0 = Equation(c * dt0 - 2 * d)
    eqn1 = Equation(c * dt - 2 * L)
    eqn2 = Equation(4 * L**2 + 4 * d**2 - v**2 * dt**2)
    eqn3 = Equation(f0 * dt0 - one)
    eqn4 = Equation(f * dt - one)
    eqns = [eqn0, eqn1, eqn2, eqn3, eqn4]
    eqnSys = EquationSystem(vars=vars, derivatives=[], constants=cons, measuredVars=vars, equations=eqns,
                                max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
    return eqnSys

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

NewtonianGravity()