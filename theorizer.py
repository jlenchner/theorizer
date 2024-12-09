from equationSystem import *
from equation import *
from m2_functions import *
from unitOfMeasure import *

def isConsistent(eqnSystem):
    """
    Function that checks the consistency of an equation system.
    Inputs:
        eqnSystem: the system of equations to check
    Outputs:
        True if consistent, False otherwise
    """

    # Extracting the variables and equations from the system 
    variables = eqnSystem.getVarNames()
    #print("Variables: " + str(variables))
    equations = eqnSystem.getEquations()
    #print("Equations: " + str(equations))
    
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

def checkConsistencyReplacedSystems(eqnSystem, index, num):
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
        eqn = eqnSystem.replaceRandomEqnByIndex(index)
        #print(f"Replaced equation {index} with: {eqn.toString()}")
        #print("Current System: \n" + eqnSystem.toString())
        if isConsistent(eqnSystem):
            #print("Consistent \n")
            numConsistent += 1
    return numConsistent

def checkConsistencyRandomSystems(vars, num):
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
        if isConsistent(eqnSystem):
            numConsistent += 1
    return numConsistent

def project(eqnSystem):
    """
    Projects the equations onto the measured variables.
    Inputs:
        eqnSystem: the system of equations to project
    Outputs:
        projectedEquations: the projected equations
    """
    # Extracting the variables and equations from the system 
    variables = eqnSystem.getVarNames()
    equations = eqnSystem.getEquations()
    measured_vars = eqnSystem.getMeasuredVars()
    non_measured_variables = eqnSystem.getNonMeasuredVars()
    
    temp_filename = "temp.txt"

    projection(variables, equations, measured_vars, non_measured_variables, temp_filename)

    # Read the results from the file
    with open(temp_filename, "r") as file:
        result = file.read().strip()

    # Remove the temporary file
    os.remove(temp_filename)

    with open("projection_output.txt", "a") as file:
        file.write(eqnSystem.toString() + "\n")
        file.write(result + "\n")
        file.write("\n")

    return result

def projectUnknownMeasuredVars(eqnSystem):
    """
    Projects the equations onto the measured variables.
    Inputs:
        eqnSystem: the system of equations to project
    Outputs:
        projectedEquations: the projected equations
    """
    # Extracting the variables and equations from the system 
    variables = eqnSystem.getVarNames()
    equations = eqnSystem.getEquations()
    measured_vars = eqnSystem.getMeasuredVars()
    
    temp_filename = "temp.txt"
    result = ""

    for i in range(1,len(measured_vars)):
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
                file.write(eqnSystem.toString() + "\n")
                file.write(result + "\n")
                file.write("\n")
            break  

    return result


def projectRandomSystems(vars, measured_vars, num):
    """
    Projects random systems of equations onto the measured variables.
    Inputs:
        num: the number of random systems to project
    Outputs:
        results: a list of results from the projections
    """
    results = []
    for i in range(num):
        eqnSystem = EquationSystem.GenerateRandom(vars, measured_vars, 4, 6)
        print("System: \n" + eqnSystem.toString())
        result = projectUnknownMeasuredVars(eqnSystem)
        results.append(result)
    return results
    

if __name__ == "__main__":
    d1,d2,m1,m2,w,p,Fc,Fg = symbols('d1,d2,m1,m2,w,p,Fc,Fg')
    syms = [Fc,Fg,w,d1,d2,m1,m2,p]
    #measured_vars = [d1,d2,m1,m2,p]
    eq2 = Equation(d1**2*Fg + 2*d1*d2*Fg + d2**2*Fg - m1*m2)
    vbles = eq2._variables
    terms = eq2.getTerms()
    m,s = symbols('m,s')
    units = UofM(m/s)

    x = Symbol('x')
    var = Variable(x)

    F,m,d,dxdt,d2xdt2 = symbols("F,m,d, dxdt, d2xdt2")
    eq3 = Equation(F - m*d2xdt2)

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

    vars = [Fc,Fg,w,d1,d2,m1,m2,p]
    # Take the measured vars to be a random subset of the variables
    #measured_vars = random.sample(vars, len(vars)//2)
    
    measured_vars = vars.copy()
    random.shuffle(measured_vars)

    projectRandomSystems(vars, measured_vars, 10) 




