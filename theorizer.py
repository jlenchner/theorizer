from equationSystem import *
from equation import *
from m2_functions import *
from unitOfMeasure import *
from derivative import *
from constant import *

    

if __name__ == "__main__":
    x, y,z ,w  = variables('x,y,z,w')
    eqn = Equation(2*x - 2*y)
    terms = eqn.getTerms()
    for term in terms:
        c = Equation.GetUnnamedConstantForTerm(term)
    eqn.divideByCommonUnnamedConstants()

    d1, d2, m1, m2, W, p, Fc, Fg, T = variables('d1,d2,m1,m2,W,p,Fc,Fg,T')
    d2xdt2 = Derivatif('d2xdt2')
    G, c = constants('G,c')
    #d1,d2,m1,m2,w,p,Fc,Fg = symbols('d1,d2,m1,m2,w,p,Fc,Fg')   #see if it is possible to do this using variables()
    vars = [Fc,Fg,W,d1,d2,m1,m2,p,T]
    derivs = [d2xdt2]
    constants = [G,c]
    term = Equation.GenerateTermFromBaseNum(31200113, vars, [])

    term1 = 3*term
    term2 = 5*term
    term3 = d1
    terms = [term, term1, term2, term3]
    cf = Equation.GetCommonFactors(terms)
    res = Equation.TermsEqualModUnnamedConstants(term1, term2)

    #vars.extend(constants)

    eqnSystem = EquationSystem.GenerateRandom(vars=vars, derivatives=derivs, measuredVars=vars, \
                                              constants=constants, numEqns=4, max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
    res = eqnSystem._equations[0].isDimensionallyConsistent()

    #for i in range(100):
    #    eqn = Equation.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=[])
    #    print(str(i+1) + ": " + str(eqn))
    for i in range(250):
        eqnSystem = EquationSystem.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivs, constants=constants,
                                              measuredVars=vars, numEqns=4, max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
        print(str(i+1) + ": " + str(eqnSystem) + "\n\n")
        res = eqnSystem._equations[0].isDimensionallyConsistent()
    #measured_vars = [d1,d2,m1,m2,p]
    eq2 = Equation(d1**2*Fg + 2*d1*d2*Fg + d2**2*Fg - m1*m2)
    res = eq2.isDimensionallyConsistent()
    eq2prime = Equation(d1 ** 2 * Fg + 2 * d1 * d2 * Fg + d2 ** 2 * Fg )
    res = eq2prime.isDimensionallyConsistent()
    vbles = eq2._variables
    terms = eq2.getTerms()
    m,s = base_units('m,s') #see below
    units = UofM(m/s)  #should not be done in terms of symbols - perhaps a new class called units, which would be veneer over Symbols

    v = Variable('v')
    dxdt = Derivatif('dxdt')
    eq = Equation(v - dxdt)
    terms2 = eq.getTerms()
    res = eq.isDimensionallyConsistent()


    F,m,a = variables('F,m,a')
    dxdt, d2xdt2 = derivatives('dxdt,d2xdt2')
    eq3 = Equation(F - m*d2xdt2)
    print(eq3)
    res = eq3.isDimensionallyConsistent()
    vars = [F,m,a]
    derivatives = [dxdt, d2xdt2]

    eqn = Equation.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivatives)
    print(str(eqn))
    res = eqn.isDimensionallyConsistent()

    #eq3 = Polynomial(Fc - m2*d2*w**2)
    #eq4 = Polynomial(Fg - Fc)
    #eq5 = Polynomial(w*p - 1)
    #eqns = [eq2, eq3, eq4, eq5]
    #eqnSystem = EquationSystem(vars=vars, measuredVars=measured_vars, equations=eqns)

    #print("System: \n" + eqnSystem.toString() + "\n")

    #isConsistent(eqnSystem)
    #print(checkConsistencyReplacedSystems(eqnSystem, 3, 10))
    #print(checkConsistencyRandomSystems(vars, 10))
    #project(eqnSystem)

    vars = [Fc,Fg,W,d1,d2,m1,m2,p]
    # Take the measured vars to be a random subset of the variables
    #measured_vars = random.sample(vars, len(vars)//2)
    
    measured_vars = vars.copy()
    random.shuffle(measured_vars)

    EquationSystem.ProjectRandomSystems(vars, derivatives, measured_vars, 10)




