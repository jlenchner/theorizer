import random
from theorizer import *
from consequence_generation import *
from consequence_data_generation import *
from system_data_generation import *
from equationSystem import *
from equation import *
from m2_functions import *
from unitOfMeasure import *
from derivative import *
from constant import *
import re
import os
from pathlib import Path
import shutil
import time
import psutil 

def save_equation_system_to_file(eqn_system_str, filename="temp.txt"):
    # Initialize dictionaries to store parsed data
    data = {
        'Variables': {'names': [], 'units': []},
        'Constants': {'names': [], 'units': []},
        'Derivatives': {'names': [], 'units': []},
        'Equations': {'expressions': [], 'units': []}
    }
    
    # Split the input string into sections
    sections = re.split(r'\n\n+', eqn_system_str.strip())
    
    # Helper function to parse units with proper handling of parentheses
    def parse_unit(unit_str):
        # Replace ** with ^ first
        unit_str = unit_str.replace('**', '^')
        return unit_str
    
    # Parse Variables section
    vars_section = sections[0]
    var_matches = re.finditer(r'(\w+) \((.*?)\)(?=\s*(?:,|$))', vars_section)
    for match in var_matches:
        data['Variables']['names'].append(match.group(1))
        data['Variables']['units'].append(parse_unit(match.group(2)))
    
    # Parse Derivatives section - improved regex to handle parentheses
    derivs_section = sections[1]
    deriv_matches = re.finditer(r'(\w+) \((.*?)\)(?=\s*(?:\w|$))', derivs_section)
    for match in deriv_matches:
        data['Derivatives']['names'].append(match.group(1))
        data['Derivatives']['units'].append(parse_unit(match.group(2)))
    
    # Parse Constants section - improved regex to handle parentheses
    consts_section = sections[2]
    const_matches = re.finditer(r'(\w+) = [\d.e+-]+ \((.*?)\)(?=\s*(?:,|$))', consts_section)
    for match in const_matches:
        data['Constants']['names'].append(match.group(1))
        data['Constants']['units'].append(parse_unit(match.group(2)))
    
    # Parse Equations section
    eqns_section = sections[3]
    # Split equations while preserving the U of M information
    eqn_blocks = re.split(r'\n(?=\S)', eqns_section.split('\n', 1)[1])
    for block in eqn_blocks:
        if '(U of M:' in block:
            # Split equation expression and unit
            parts = re.split(r'\(U of M:\s*(.*?)\)\s*$', block)
            expr = parts[0].strip()
            unit = parse_unit(parts[1])
            data['Equations']['expressions'].append(expr.replace('**', '^'))
            data['Equations']['units'].append(unit)
    
    # Write to file
    with open(filename, 'w') as f:
        # Write Variables
        f.write(f"Variables: {data['Variables']['names']}\n")
        
        # Write Constants
        f.write(f"Constants: {data['Constants']['names']}\n")
        
        # Write Derivatives
        f.write(f"Derivatives: {data['Derivatives']['names']}\n")
        
        # Write Equations
        f.write("Equations:\n")
        for expr in data['Equations']['expressions']:
            f.write(expr + '\n')
        
        # Write Units of Measure
        f.write(f"Units of Measure of Variables: {data['Variables']['units']}\n")
        f.write(f"Units of Measure of Constants: {data['Constants']['units']}\n")
        f.write(f"Units of Measure of Derivatives: {data['Derivatives']['units']}\n")
        
        # Write Equations Units
        f.write("Units of Measure of Equations:\n")
        for unit in data['Equations']['units']:
            f.write(unit + '\n')

def create_replacement_files(replacement_info, num_replacements, original_file="temp.txt"):
    # Read the original file
    with open(original_file, 'r') as f:
        original_content = f.read()
    
    # Parse the replacement info
    replacements = []
    current_replacement = {}
    lines = replacement_info.split('\n')
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.isdigit():  # This is an equation index
            if current_replacement:  # Save previous replacement if exists
                replacements.append(current_replacement)
            current_replacement = {'index': int(line) - 1}  # Convert to 0-based index
        elif '(U of M:' in line:  # Unit of measure
            # Extract the equation and unit
            parts = line.split('(U of M:')
            current_replacement['equation'] = parts[0].strip().replace('**', '^')
            unit = parts[1].split(')')[0].strip()
            current_replacement['unit'] = unit.replace('**', '^')
    
    if current_replacement:  # Add the last replacement
        replacements.append(current_replacement)
    
    # Ensure we have exactly num_replacements
    replacements = replacements[:num_replacements]
    
    # Parse the original file structure
    sections = {}
    current_section = None
    
    for line in original_content.split('\n'):
        if line.startswith('Variables:') or line.startswith('Constants:') or line.startswith('Derivatives:'):
            current_section = line.split(':')[0]
            sections[current_section] = line + '\n'
        elif line.startswith('Equations:'):
            current_section = 'Equations'
            sections[current_section] = [line + '\n']  # Start equations list
        elif line.startswith('Units of Measure of'):
            current_section = line
            sections[current_section] = line + '\n'
        elif current_section == 'Equations' and line.strip():
            sections['Equations'].append(line + '\n')
        elif current_section and current_section.startswith('Units of Measure of'):
            sections[current_section] += line + '\n'
    
    # For each replacement, create a new file
    for k, replacement in enumerate(replacements[:num_replacements]):
        # Create a copy of the original content
        new_content = {k: v for k, v in sections.items()}
        
        # Replace the specified equation
        eqn_index = replacement['index']
        if eqn_index < len(new_content['Equations']) - 1:  # -1 because first line is "Equations:"
            new_content['Equations'][eqn_index + 1] = replacement['equation'] + '\n'
        
        # Replace the corresponding unit of measure
        units_line = None
        for section in new_content:
            if section.startswith('Units of Measure of Equations'):
                units_line = section
                break
        
        if units_line:
            units = new_content[units_line].split('\n')[1:]  # Skip header line
            if eqn_index < len(units):
                units[eqn_index] = replacement['unit'] + '\n'
                new_content[units_line] = 'Units of Measure of Equations:\n' + ''.join(units)
        
        # Reconstruct the file content
        output_lines = []
        output_lines.append(new_content['Variables'])
        output_lines.append(new_content['Constants'])
        output_lines.append(new_content['Derivatives'])
        output_lines.extend(new_content['Equations'])
        
        # Add all units sections
        for section in new_content:
            if section.startswith('Units of Measure of') and section not in ['Variables', 'Constants', 'Derivatives', 'Equations']:
                output_lines.append(new_content[section])
        
        # Write to replacement file
        replacement_file = f"temp_replacement_{k+1}.txt"
        with open(replacement_file, 'w') as f:
            f.write(''.join(output_lines))
        
        print(f"Created replacement file: {replacement_file}")

def cleanup_temp_files(system_num):
    """Clean up all temporary files for a given system number"""
    patterns = [
        f"temp_system_{system_num}.txt",
        f"temp_system_{system_num}.dat",
        f"temp_system_{system_num}_*.dat",
        f"temp_consequence_{system_num}.txt",
        f"temp_consequence_{system_num}_reordered.txt",
        f"temp_data_{system_num}.dat",
        f"temp_data_{system_num}_*.dat",
        "temp_replacement_*.txt"
    ]
    
    for pattern in patterns:
        for f in Path(".").glob(pattern):
            try:
                f.unlink()
            except:
                pass

def get_next_system_num():
    system_num = 1
    while Path(f"Benchmarking/System_{system_num}").exists():
        system_num += 1
    return system_num

def main():
    # Start timing
    start_time = time.time()
    process = psutil.Process(os.getpid())
    Equation.SetLogging()

    # Define parameter ranges
    num_vars_options = [6, 7, 8, 9]        # Approximate range 6-10
    num_derivs_options = [2, 3, 4]       # Range 2-4
    num_eqns_options = [4, 5]         # Range 4-6
    numReplacements = 0

    # Create Benchmarking directory if it doesn't exist
    Path("Benchmarking").mkdir(exist_ok=True)

    # Track completed systems per configuration
    config_counts = {}

    # Generate all parameter combinations
    param_combinations = list(itertools.product(num_vars_options, num_derivs_options, num_eqns_options))

    for num_vars, num_derivs, num_eqns in param_combinations:
        config_key = f"vars_{num_vars}_derivs_{num_derivs}_eqns_{num_eqns}"
        config_dir = f"Benchmarking/{config_key}"
        Path(config_dir).mkdir(exist_ok=True)

        existing_systems = 0
        if Path(config_dir).exists():
            existing_systems = len([f for f in os.listdir(config_dir) if f.startswith("System_")])
        
        config_counts[config_key] = existing_systems
        
        print(f"\n=== Starting configuration: {num_vars} vars, {num_derivs} derivs, {num_eqns} eqns ===")
        attempts = 1

        while config_counts[config_key] < 3 and attempts < 10:
            print("Attempt ", attempts, "\n")
            attempts+=1
            # Generate variables and derivatives dynamically based on current config
            var_names = ['Fc', 'Fg', 'W', 'd1', 'd2', 'm1', 'm2', 'p', 'T', 'E'][:num_vars]
            deriv_names = ['dx1dt', 'd2x1dt2', 'dx2dt', 'd2x2dt2'][:num_derivs]
            
            # Create the symbols
            vars = variables(','.join(var_names))
            derivs = derivatives(','.join(deriv_names))
            G, c, pi = constants('G,c,pi')
            
            print(f"\nGenerating system {config_counts[config_key]+1}/3 for config {config_key}")

            # Track time and memory at key stages
            stage_times = {}
            mem_usage = {}


            # Find the next available system number within this configuration
            def get_next_system_num_with_config():
                system_num = 1
                while Path(f"{config_dir}/System_{system_num}").exists():
                    system_num += 1
                return system_num
            
            system_num = get_next_system_num_with_config()
            
            try:
                t = time.time()
                eqnSystem = EquationSystem.GenerateRandomDimensionallyConsistent(
                    vars=vars, 
                    derivatives=derivs,
                    constants=[G, c],
                    measuredVars=vars, 
                    numEqns=num_eqns,
                    max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX
                )
                eqn_system_str = str(eqnSystem)
                stage_times['generation'] = time.time() - t
                mem_usage['generation'] = process.memory_info().rss / (1024 * 1024)

                print(f"System {system_num} generated.")
                print(f"Variables: {var_names}")
                print(f"Derivatives: {deriv_names}")
                print(f"Equations: {num_eqns}")

                temp_system_file = f"temp_system_{system_num}.txt"
                save_equation_system_to_file(eqn_system_str, temp_system_file)

                print(" Now generating replacement axioms. \n")
                # Generate replacements first
                replacement_info = ""
                t = time.time()
                if numReplacements > 0:
                    for j in range(numReplacements):
                        index, eqn = eqnSystem.copy().replaceRandomDimensionallyConsistentEqn()
                        replacement_info += f"{index+1}\n{eqn}\n"  # +1 to make it 1-based index
                        print(f"Replacement {j+1}: Equation {index+1} with {eqn}")
                stage_times['replacement_generation'] = time.time() - t
                mem_usage['replacement_generation'] = process.memory_info().rss / (1024 * 1024)

                # Generate replacement files if needed
                if numReplacements > 0:
                    create_replacement_files(replacement_info, numReplacements, temp_system_file)

                print("Replacement Axioms Created. Now searching for roots. \n")
                
                # Step 1: Generate noiseless system data
                temp_system_dat = f"temp_system_{system_num}.dat"
                t = time.time()
                found_system_roots = run_noiseless_system_data_generation(temp_system_file, temp_system_dat)
                stage_times['system_root_search'] = time.time() - t
                mem_usage['system_root_search'] = process.memory_info().rss / (1024 * 1024)
                if not found_system_roots:
                    print("Could not find roots to the system. Skipping. \n")
                    #cleanup_temp_files(system_num)
                    continue
                print("System Data Generated. Now adding noise \n")
                    
                # Step 2: Add noise to system data
                t = time.time()
                added_system_noise = run_noisy_system_data_generation(temp_system_file, temp_system_dat)
                stage_times['system_noise_addition'] = time.time() - t
                mem_usage['system_noise_addition'] = process.memory_info().rss / (1024 * 1024)
                if not added_system_noise:
                    print("Failed to add noise")
                    #cleanup_temp_files(system_num)
                    continue
                print("Noisy System Data Generated. Now searching for consequence \n")
                    
                # Step 3: Generate consequence
                temp_consequence = f"temp_consequence_{system_num}.txt"
                t = time.time()
                found_consequence = run_consequence_generation(temp_system_file, temp_consequence)
                stage_times['consequence_generation'] = time.time() - t
                mem_usage['consequence_generation'] = process.memory_info().rss / (1024 * 1024)
                if not found_consequence:
                    print("Could not find a consequence. Skipping. \n")
                    #cleanup_temp_files(system_num)
                    continue
                print("Consequence generated. Now generating noiseless data for consequence \n")

                # Step 4: Generate noiseless consequence data
                temp_data_dat = f"temp_data_{system_num}.dat"
                t = time.time()
                found_roots = run_consequence_noiseless_data_generation(temp_consequence, temp_data_dat)
                stage_times['consequence_data_generation'] = time.time() - t
                mem_usage['consequence_data_generation'] = process.memory_info().rss / (1024 * 1024)
                if not found_roots:
                    print("Could not solve equation for roots. Skipping")
                    #cleanup_temp_files(system_num)
                    continue
                

                reordered_consequence = temp_consequence.replace('.txt', '_reordered.txt')
                if os.path.exists(reordered_consequence):
                    # Remove the original consequence file
                    if os.path.exists(temp_consequence):
                        os.remove(temp_consequence)
                    # Rename the reordered file to take its place
                    os.rename(reordered_consequence, temp_consequence)
                    print("Used reordered consequence file")
                
                print("Noiseless Data Generated. Now adding noise \n")

                # Step 5: Add noise to consequence data
                t = time.time()
                run_consequence_noisy_data_generation(temp_data_dat)
                stage_times['consequence_data_noise_addition'] = time.time() - t
                mem_usage['consequence_data_noise_addition'] = process.memory_info().rss / (1024 * 1024)
                print("Generated Noisy Data. System Complete. \n")

                total_time = time.time() - start_time
                peak_memory = max(mem_usage.values())
                print("\nPerformance Summary:")
                print(f"Total time: {total_time:.2f} seconds")
                print(f"Peak memory: {peak_memory:.2f} MB")
                print("Stage times:")
                for stage, t in stage_times.items():
                    print(f"  {stage}: {t:.2f}s")
                
                # All steps succeeded - create system directory within config directory
                system_dir = f"{config_dir}/System_{system_num}"
                Path(system_dir).mkdir(exist_ok=True)
                
                # Move system files
                shutil.move(temp_system_file, f"{system_dir}/system.txt")
                shutil.move(temp_system_dat, f"{system_dir}/system.dat")
                
                # Move noisy system data files
                for f in Path(".").glob(f"temp_system_{system_num}_*.dat"):
                    noise_level = f.stem.split('_')[-1]
                    shutil.move(str(f), f"{system_dir}/system_{noise_level}.dat")
                    
                # Move consequence files
                shutil.move(temp_consequence, f"{system_dir}/consequence.txt")
                shutil.move(temp_data_dat, f"{system_dir}/consequence.dat")
                
                # Move noisy consequence data files
                for f in Path(".").glob(f"temp_data_{system_num}_*.dat"):
                    noise_level = f.stem.split('_')[-1]
                    shutil.move(str(f), f"{system_dir}/consequence_{noise_level}.dat")
                
                # Move replacement files if they exist
                for f in Path(".").glob(f"temp_replacement_*.txt"):
                    replacement_num = f.stem.split('_')[-1]
                    shutil.move(str(f), f"{system_dir}/replacement_{replacement_num}.txt")
                
                # Save Performance to file
                with open(f"{system_dir}/performance.txt", "w") as f:
                    f.write(f"Total time: {total_time:.2f} seconds\n")
                    f.write(f"Peak memory: {peak_memory:.2f} MB\n")
                    for stage, t in stage_times.items():
                        f.write(f"{stage}_time: {t:.2f}\n")
                    for stage, m in mem_usage.items():
                        f.write(f"{stage}_memory: {m:.2f} MB\n")
                
                config_counts[config_key] += 1
                print(f"Successfully completed system {config_counts[config_key]}/3 for this configuration")
                attempts = 1

            except Exception as e:
                print(f"Error generating system: {str(e)}")
                cleanup_temp_files(system_num)
                continue

"""
def old_main_for_individual_systems():
    # Start timing
    start_time = time.time()
    process = psutil.Process(os.getpid())
    Equation.SetLogging()

    d1, d2, m1, m2, W, p, Fc, Fg, T, E = variables('d1,d2,m1,m2,W,p,Fc,Fg,T,E')
    dxdt, d2xdt2 = derivatives('dxdt,d2xdt2')
    dx1dt, d2x1dt2, dx2dt, d2x2dt2 = derivatives('dx1dt,d2x1dt2,dx2dt,d2x2dt2')
    G, c, pi = constants('G,c,pi')
    vars = [Fc, Fg, W, d1, d2, m1, m2, p]
    derivs = [dx1dt, d2x1dt2, dx2dt, d2x2dt2]
    consts = [G, c]
    numReplacements = 0

    print("Dimensionally Consistent Set:\n")
    Path("Benchmarking").mkdir(exist_ok=True)
    
    # Find the next available system number
    def get_next_system_num():
        system_num = 1
        while Path(f"Benchmarking/System_{system_num}").exists():
            system_num += 1
        return system_num
    
    for _ in range(10):  
        
        system_num = get_next_system_num()
        print(f"Generating system {system_num} \n")

        # Track time and memory at key stages
        stage_times = {}
        mem_usage = {}

        t = time.time()
        eqnSystem = EquationSystem.GenerateRandomDimensionallyConsistent(vars=vars, derivatives=derivs,
                                                                     constants=consts,
                                                                     measuredVars=vars, numEqns=4,
                                                                     max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX)
        eqn_system_str = str(eqnSystem)
        stage_times['generation'] = time.time() - t
        mem_usage['generation'] = process.memory_info().rss / (1024 * 1024)

        print(f"{system_num}: {eqn_system_str}\n")
        print(f"System {system_num} generated. Now generating replacement axioms. \n")

        temp_system_file = f"temp_system_{system_num}.txt"
        save_equation_system_to_file(eqn_system_str, temp_system_file)

        # Generate replacements first
        replacement_info = ""
        t = time.time()
        if numReplacements > 0:
            for j in range(numReplacements):
                index, eqn = eqnSystem.copy().replaceRandomDimensionallyConsistentEqn()
                replacement_info += f"{index+1}\n{eqn}\n"  # +1 to make it 1-based index
                print(f"Replacement {j+1}: Equation {index+1} with {eqn}")
        stage_times['replacement_generation'] = time.time() - t
        mem_usage['replacement_generation'] = process.memory_info().rss / (1024 * 1024)

        # Generate replacement files if needed
        if numReplacements > 0:
            create_replacement_files(replacement_info, numReplacements, temp_system_file)

        print("Replacement Axioms Created. Now searching for roots. \n")
        

        # Step 1: Generate noiseless system data
        temp_system_dat = f"temp_system_{system_num}.dat"
        t = time.time()
        found_system_roots = run_noiseless_system_data_generation(temp_system_file, temp_system_dat)
        stage_times['system_root_search'] = time.time() - t
        mem_usage['system_root_search'] = process.memory_info().rss / (1024 * 1024)
        if not found_system_roots:
            print("Could not find roots to the system. Skipping. \n")
            #cleanup_temp_files(system_num)
            continue
        print("System Data Generated. Now adding noise \n")

            
        # Step 2: Add noise to system data
        t = time.time()
        added_system_noise = run_noisy_system_data_generation(temp_system_file, temp_system_dat)
        stage_times['system_noise_addition'] = time.time() - t
        mem_usage['system_noise_addition'] = process.memory_info().rss / (1024 * 1024)
        if not added_system_noise:
            print("Failed to add noise")
            #cleanup_temp_files(system_num)
            continue
        print("Noisy System Data Generated. Now searching for consequence \n")
        
        # Step 3: Generate consequence
        temp_consequence = f"temp_consequence_{system_num}.txt"
        t = time.time()
        found_consequence = run_consequence_generation(temp_system_file, temp_consequence)
        stage_times['consequence_generation'] = time.time() - t
        mem_usage['consequence_generation'] = process.memory_info().rss / (1024 * 1024)
        if not found_consequence:
            print("Could not find a consequence. Skipping. \n")
            #cleanup_temp_files(system_num)
            continue
        print("Consequence generated. Now generating noiseless data for consequence \n")
        
        # Step 4: Generate noiseless consequence data
        temp_data_dat = f"temp_data_{system_num}.dat"
        t = time.time()
        found_roots = run_consequence_noiseless_data_generation(temp_consequence, temp_data_dat)
        stage_times['consequence_data_generation'] = time.time() - t
        mem_usage['consequence_data_generation'] = process.memory_info().rss / (1024 * 1024)
        if not found_roots:
            print("Could not solve equation for roots. Skipping")
            #cleanup_temp_files(system_num)
            continue
        
        reordered_consequence = temp_consequence.replace('.txt', '_reordered.txt')
        if os.path.exists(reordered_consequence):
            # Remove the original consequence file
            if os.path.exists(temp_consequence):
                os.remove(temp_consequence)
            # Rename the reordered file to take its place
            os.rename(reordered_consequence, temp_consequence)
            print("Used reordered consequence file")
        
        print("Noiseless Data Generated. Now adding noise \n")

        # Step 5: Add noise to consequence data
        t = time.time()
        run_consequence_noisy_data_generation(temp_data_dat)
        stage_times['consequence_data_noise_addition'] = time.time() - t
        mem_usage['consequence_data_noise_addition'] = process.memory_info().rss / (1024 * 1024)
        print("Generated Noisy Data. System Complete \n")

        total_time = time.time() - start_time
        peak_memory = max(mem_usage.values())
        print("\nPerformance Summary:")
        print(f"Total time: {total_time:.2f} seconds")
        print(f"Peak memory: {peak_memory:.2f} MB")
        print("Stage times:")
        for stage, t in stage_times.items():
            print(f"  {stage}: {t:.2f}s")
        
        # All steps succeeded - create directory and move files
        system_dir = f"Benchmarking/System_{system_num}"
        Path(system_dir).mkdir(exist_ok=True)
        
        # Move system files
        shutil.move(temp_system_file, f"{system_dir}/system.txt")
        shutil.move(temp_system_dat, f"{system_dir}/system.dat")
        
        # Move noisy system data files
        for f in Path(".").glob(f"temp_system_{system_num}_*.dat"):
            noise_level = f.stem.split('_')[-1]
            shutil.move(str(f), f"{system_dir}/system_{noise_level}.dat")
            
        # Move consequence files
        shutil.move(temp_consequence, f"{system_dir}/consequence.txt")
        shutil.move(temp_data_dat, f"{system_dir}/consequence.dat")
        
        # Move noisy consequence data files
        for f in Path(".").glob(f"temp_data_{system_num}_*.dat"):
            noise_level = f.stem.split('_')[-1]
            shutil.move(str(f), f"{system_dir}/consequence_{noise_level}.dat")
        
        # Move replacement files if they exist
        for f in Path(".").glob(f"temp_replacement_*.txt"):
            replacement_num = f.stem.split('_')[-1]
            shutil.move(str(f), f"{system_dir}/replacement_{replacement_num}.txt")
        
        # Save Performance to file
        with open(f"Benchmarking/System_{system_num}/performance.txt", "w") as f:
            f.write(f"Total time: {total_time:.2f} seconds\n")
            f.write(f"Peak memory: {peak_memory:.2f} MB\n")
            for stage, t in stage_times.items():
                f.write(f"{stage}_time: {t:.2f}\n")
            for stage, m in mem_usage.items():
                f.write(f"{stage}_memory: {m:.2f} MB\n")
"""

if __name__ == "__main__":
    main()