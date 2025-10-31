import subprocess
import os

def projection(variables, axioms, measured_variables, non_measured_variables, filename):
    """
    Compute the algebraic projection (elimination) of an ideal using Macaulay2.

    Given:
      - variables: All ring variables (strings usable as M2 symbols).
      - axioms: Axioms as strings in M2 syntax (e.g., "x^2 + y - 1").
      - measured_variables: Subset of variables to keep after elimination.
      - non_measured_variables: Variables to eliminate.
      - filename: Destination file to write human-readable results to.

    Side effects:
      - Writes a result file (gröbner bases + eliminated basis polynomials).
      - Prints M2 stdout; prints M2 stderr if present.

    Returns:
      The textual contents of `filename` on success, or an error message string on failure.
    """
    
    macaulay2_path = "/opt/homebrew/bin/M2"

    # Constructing the Macaulay2 script
    macaulay2_script = f"""
-- Define the ring and ideal
print("Beginning projection process")
print("Variables: " | toString {','.join(variables)})

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

    # Use a temporary script and output file; move to `filename` on success.
    script_path = "temp_script.m2"
    with open(script_path, "w") as file:
        file.write(macaulay2_script)

    # Call Macaulay2 with the script
    result = subprocess.run([macaulay2_path, "--script", script_path], capture_output=True, text=True)

    # Remove the temporary script file
    os.remove(script_path)

    # (Optional) Uncomment to print Macaulay2 logs
    """
    print("Output from Macaulay2:")
    print(result.stdout)
    if result.stderr:
        print("Errors:")
        print(result.stderr)
    """

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
    """
    Check ideal consistency (non-empty variety) with Macaulay2.

    Considers the ideal I = <axioms> in QQ[variables] with Lex order.
    A common convention: dim(I) == -1 ↔ inconsistent (empty variety).

    Args:
      variables: Ring variables in M2 syntax.
      axioms: Axioms as M2 expressions (strings).
      filename: Destination text file written with either "consistent." or "inconsistent."
      macaulay2_path: Optional explicit path to M2 binary.

    Returns:
      The textual contents of `filename` on success, or an error message string on failure.
    """
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
    
    # Use a temporary script and output file; move to `filename` on success.
    script_path = "temp_script.m2"
    with open(script_path, "w") as file:
        file.write(macaulay2_script)

    # Call Macaulay2 with the script
    result = subprocess.run([macaulay2_path, "--script", script_path], capture_output=True, text=True)

    # Remove the temporary script file
    os.remove(script_path)

    # (Optional) Uncomment to print Macaulay2 logs
    """
    print("Output from Macaulay2:")
    print(result.stdout)
    if result.stderr:
        print("Errors:")
        print(result.stderr)
    """
    # Read and return the contents of the saved file
    try:
        with open(filename, "r") as file:
            return file.read()
    except FileNotFoundError:
        return f"Error: The file '{filename}' was not found."
    except IOError as e:
        return f"Error reading the file: {e}"
