from equationSystem import *
from equation import *
from m2_functions import *
from unitOfMeasure import *
from derivative import *

    

if __name__ == "__main__":
    d1, d2, m1, m2, W, p, Fc, Fg = variables('d1,d2,m1,m2,W,p,Fc,Fg')
    #d1,d2,m1,m2,w,p,Fc,Fg = symbols('d1,d2,m1,m2,w,p,Fc,Fg')   #see if it is possible to do this using variables()
    vars = [Fc,Fg,W,d1,d2,m1,m2,p]
    term = Equation.GenerateTermFromBaseNum(31200113, vars, [])

    term1 = 3*term
    term2 = 5*term
    term3 = d1
    terms = [term, term1, term2, term3]
    cf = Equation.GetCommonFactors(terms)
    res = Equation.TermsEqualModConstants(term1, term2)

    eqnSystem = EquationSystem.GenerateRandom(vars=vars, derivatives=[], measuredVars=vars, \
                                              numEqns=4, max_var_and_derivatives_per_eqn=4)
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



















