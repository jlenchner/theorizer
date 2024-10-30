import subprocess
import os

def abductive_inference(variables, measured_variables, non_measured_variables, axioms, q, filename, removed_axiom):
    # Macaulay2 executable path
    macaulay2_path = "/opt/homebrew/bin/M2"

    # Construct the Macaulay2 script
    macaulay2_script = f"""
needsPackage("PrimaryDecomposition", Reload => true)

-- Define the ring
R = QQ[{','.join(variables)}, MonomialOrder => Lex];
print ("Ring defined: " | toString R);

-- Locally defining the target polynomial
q = {q};
print ("Target polynomial defined: " | toString q);

-- Locally defining the axioms
axioms = toList([{','.join(axioms)}]);
print ("Axioms defined: " | toString axioms);

-- Locally defining the removed axiom
removedAxiom = {removed_axiom};
print ("Removed axiom defined: " | toString removedAxiom);

-- Defining the measured variables
measuredVariables = toList([{','.join(measured_variables)}]);
print("Measured variables defined: " | toString measuredVariables);

-- Defining the non-measured variables
nonMeasuredVariables = toList([{','.join(non_measured_variables)}]);
print("Non-measured variables defined: " | toString nonMeasuredVariables);

-- Define the ideal
I = ideal(append(axioms, q));
print ("Ideal defined: " | toString I);

-- Compute primary decomposition
PD = primaryDecomposition I;
print ("Primary decomposition computed");

isInIdeal = (q,J) -> (
    M = ideal(append(axioms,J));
    G = gens gb M;
    result = q % ideal(G) == 0;
    return result;
);

-- Function to check if q is in the Gröbner basis of an ideal
isInGB = (q,J) -> (
    M = ideal(append(axioms,J));
    eliminatedIdeal = eliminate(nonMeasuredVariables, M);
    GBproj = gens gb eliminatedIdeal;
    result = toList apply(flatten entries GBproj, g -> g == q);
    return member(true, result)
);

-- Open a file to write the results
f = openOut "decomposition_analysis.txt"

-- Write the removed axiom to the file
f << "Axiom removed: " << toString removedAxiom << endl << endl;

-- Write the primary decomposition to the file 
f << "Primary Decomposition:" << endl;
scan(PD, J -> f << toString J << endl << endl);

-- Write the primary decomposition to the file and check generators
f << "Primary Decomposition Analysis:" << endl;
savedPolys = {{}};
strongCandidates = {{}};
scan(PD, J -> (
    f << "Ideal: " << toString J << endl;
    idealGens = first entries gens J;
    scan(idealGens, g -> (
        if isInIdeal(q,g) then (
            f << "  Generator " << toString g << " implies q" << endl;
            savedPolys = append(savedPolys, g);
            if isInGB(q,g) then (
                f << "  Generator " << toString g << " is a strong candidate" << endl;
                strongCandidates = append(strongCandidates, g);
            );
        ) else (
            f << "  Generator " << toString g << " does not imply q" << endl;
        );
    ));
    f << endl;
));

-- Write saved polynomials
f << "Saved Polynomials:" << endl;
scan(savedPolys, g -> f << toString g << endl);

-- Write empty line
f << endl;

-- Write strong candidates
f << "Strong Candidates:" << endl;
scan(strongCandidates, g -> f << toString g << endl);

close f;

print "Results saved to decomposition_analysis.txt"
"""

    # Write the Macaulay2 script to a temporary file
    script_path = "temp_script.m2"
    with open(script_path, "w") as file:
        file.write(macaulay2_script)

    # Call Macaulay2 with the script
    result = subprocess.run([macaulay2_path, "--script", script_path], capture_output=True, text=True)

    # Remove the temporary script file
    os.remove(script_path)

    # Print the output from Macaulay2
    print("Output from Macaulay2:")
    print(result.stdout)

    # Print any errors
    if result.stderr:
        print("Errors:")
        print(result.stderr)

    # Read and return the contents of the saved file
    try:
        with open("decomposition_analysis.txt", "r") as file:
            os.rename("decomposition_analysis.txt", filename)
            return file.read()
    except FileNotFoundError:
        return "Error: The file 'decomposition_analysis.txt' was not found."
    except IOError as e:
        return f"Error reading the file: {e}"

def projection(variables, axioms, measured_variables, non_measured_variables, filename):
    # Macaulay2 executable path
    macaulay2_path = "/opt/homebrew/bin/M2"

    # Construct the Macaulay2 script
    macaulay2_script = f"""
-- Define the ring and ideal
R = QQ[{','.join(variables)}, MonomialOrder => Lex];
print("Ring defined: " | toString R);

-- Defining the axioms
axioms = toList([{','.join(axioms)}]);
print("Axioms defined: " | toString axioms);

-- Defining the measured variables
measuredVariables = toList([{','.join(measured_variables)}]);
print("Measured variables defined: " | toString measuredVariables);

-- Defining the non-measured variables
nonMeasuredVariables = toList([{','.join(non_measured_variables)}]);
print("Non-measured variables defined: " | toString nonMeasuredVariables);

-- Define the ideal
I = ideal(axioms);
print("Ideal defined: " | toString I);

-- Compute the Gröbner basis of the ideal
GB = gens gb I;
print("Gröbner basis: " | toString GB);

-- Define the Ring of Measured Variables
Rproj = QQ[measuredVariables, MonomialOrder => Lex];
print("Ring of Measured Variables defined: " | toString Rproj);

-- Define the projection map
eliminatedIdeal = eliminate(nonMeasuredVariables, I);
print("Eliminated Ideal: " | toString eliminatedIdeal);

-- Compute the Gröbner basis of the eliminated ideal
GBproj = gens gb eliminatedIdeal;
print("Gröbner basis of the eliminated ideal: " | toString GBproj);

-- Write the results to a file
f = openOut "projection_analysis.txt";
f << "Gröbner basis of the ideal:" << endl;
f << toString GB << endl << endl;
f << "Gröbner basis of the eliminated ideal:" << endl;
f << toString GBproj << endl;
f << endl;

-- Save the polynomials of the Gröbner basis of the eliminated ideal
f << "Polynomials of the Gröbner basis of the eliminated ideal:" << endl;
scan(flatten entries GBproj, g -> f << toString g << endl);

close f;

"""

    # Write the Macaulay2 script to a temporary file
    script_path = "temp_script.m2"
    with open(script_path, "w") as file:
        file.write(macaulay2_script)

    # Call Macaulay2 with the script
    result = subprocess.run([macaulay2_path, "--script", script_path], capture_output=True, text=True)

    # Remove the temporary script file
    os.remove(script_path)

    # Print the output from Macaulay2
    print("Output from Macaulay2:")
    print(result.stdout)

    # Print any errors
    if result.stderr:
        print("Errors:")
        print(result.stderr)

    # Read and return the contents of the saved file
    try:
        with open("projection_analysis.txt", "r") as file:
            os.rename("projection_analysis.txt", filename)
            return file.read()
    except FileNotFoundError:
        return "Error: The file 'projection_analysis.txt' was not found."
    except IOError as e:
        return f"Error reading the file: {e}"
    
def check_consistency(variables, axioms, filename):
    # Macaulay2 executable path
    macaulay2_path = "/opt/homebrew/bin/M2"

    # Construct the Macaulay2 script
    macaulay2_script = f"""
-- Define the ring and ideal
R = QQ[{','.join(variables)}, MonomialOrder => Lex];
-- print("Ring defined: " | toString R);

-- Defining the axioms
axioms = toList([{','.join(axioms)}]);
-- print("Axioms defined: " | toString axioms);

-- Define the ideal
I = ideal(axioms);
-- print("Ideal defined: " | toString I);

-- Check for consistency
isInConsistent = (I) -> (
    return dim I == -1;
);

-- Open a file to write the results
f = openOut "{filename}";

-- Check consistency and write the result
if isInConsistent(I) then (
    f << "inconsistent." << endl;
) else (
    f << "consistent." << endl;
);

close f;

"""
    
    # Write the Macaulay2 script to a temporary file
    script_path = "temp_script.m2"
    with open(script_path, "w") as file:
        file.write(macaulay2_script)

    # Call Macaulay2 with the script
    result = subprocess.run([macaulay2_path, "--script", script_path], capture_output=True, text=True)

    # Remove the temporary script file
    os.remove(script_path)

    # Print the output from Macaulay2
    #print("Output from Macaulay2:")
    #print(result.stdout)

    # Print any errors
    if result.stderr:
        print("Errors:")
        print(result.stderr)

    # Read and return the contents of the saved file
    try:
        with open(filename, "r") as file:
            return file.read()
    except FileNotFoundError:
        return f"Error: The file '{filename}' was not found."
    except IOError as e:
        return f"Error reading the file: {e}"


# Example Usage
vars=["d", "dt", "dt0", "L", "c", "f0", "f", "v"]
measured_vars=["c", "f0", "f", "v"]
axioms=["c*dt0-2*d", "4*L^2-4*d^2-v^2*dt^2", "f0*dt0-1", "f*dt-1", "c*dt-2*L"]
target="c^2*f0^2-c^2*f^2-f0^2*v^2"
target_var="f"
non_measured_variables = list(set(vars) - set(measured_vars))

#projection(vars, axioms, measured_vars, non_measured_variables, "projection_analysis.txt")
#abductive_inference(vars, measured_vars, non_measured_variables, axioms, target, "decomposition_analysis.txt", "c*dt0-2*d")
#check_consistency(vars, axioms, "consistency_analysis.txt")