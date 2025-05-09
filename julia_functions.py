import subprocess
import tempfile
import os
import json

def ai_hilbert_axioms_and_data(env_path, variables, axioms, measured_variables, target_var, deg, deg_overall, deg_overall_q, deg_elementwise, deg_elementwise_q, target, filename, data_name):

    julia_code = """
    using DynamicPolynomials
    using JuMP
    using MosekTools
    using LinearAlgebra
    using CSV
    using DataFrames
    import JSON
    using Dates
    using DelimitedFiles
    using Plots
    using Statistics
    using Gurobi


    function all_monomials_up_to_max_deg(x, deg)
        if size(x,1) == 0
            [1]
        else
        [ x[1]^k * m for k=0:deg 
                for m=all_monomials_up_to_max_deg(x[2:end], deg)
        ]
        end
    end

    function mons_of_max_degree_and_unit(x, deg)
        [m
            for m=all_monomials_up_to_max_deg(x, deg)
        ]
    end

    function all_monomials_up_to_max_deg_overall(x, deg, deg_overall, deg_elem)
        if size(x,1) == 0
            [1]
        else
        [ x[1]^k * m for k=0:min(deg, deg_overall, deg_elem[1]) 
                    for m=all_monomials_up_to_max_deg_overall(x[2:end], deg, deg_overall-k, deg_elem[2:end])
        ]
        end
    end

    function mons_of_max_degree_and_unit_overall(x, deg, deg_overall, deg_elem)
        [m
            for m=all_monomials_up_to_max_deg_overall(x, deg, deg_overall, deg_elem)
        ]
    end
                    
    function degree_poly(p)
        maximum(degree.(monomials(p)))
    end

    function save_polynomial_to_file(poly, filename, description)
        open(filename*"_poly.txt", "a") do file
            write(file, string(description, ": ", poly) * "\n") 
        end
    end


    function create_variables(var_strings)
        variables = Symbol[]
        for var_string in var_strings
            # Remove any whitespace and non-alphanumeric characters (except underscore)
            cleaned_var = replace(strip(var_string), r"[^a-zA-Z0-9_]" => "")
            
            # Ensure the variable doesn't start with a number
            if !isempty(cleaned_var) && !isdigit(cleaned_var[1])
                push!(variables, Symbol(cleaned_var))
            end
        end
        return variables
    end

    

    
    variables = """ + json.dumps(variables) + """
    axiom_strings = """ + json.dumps(axioms) + """
    measured_variables = """ + json.dumps(measured_variables) + """
    target_variable = """ + json.dumps(target_var) + """
    target = """ + json.dumps(target) + """

    eval(:((@polyvar $(Symbol.(variables)...))))
    var_dict = Dict(var => eval(Symbol(var)) for var in variables)
    vars = [eval(Symbol(var)) for var in variables]

    function replace_vars(expr)
        for (var, poly) in var_dict
            expr = replace(expr, var => string(poly))
        end
        return expr
    end

    axioms = [eval(Meta.parse(replace_vars(axiom))) for axiom in axiom_strings]
    
    measured_vars = [eval(Symbol(var)) for var in measured_variables]

    target_poly = eval(Meta.parse(replace_vars(target)))

    deg_elementwise = """ + json.dumps(deg_elementwise) + """
    deg_elementwise_q = """ + json.dumps(deg_elementwise_q) + """
    deg = """ + str(deg) + """
    deg_overall = """ + str(deg_overall) + """
    deg_overall_q = """ + str(deg_overall_q) + """

    target_var = eval(Symbol(target_variable))

    println("Variables:", vars)
    println("Axioms:", axioms)
    println("Measured variables:", measured_vars)
    println("Degree:", deg)
    println("Overall degree:", deg_overall)
    println("Overall degree q:", deg_overall_q)
    println("Elementwise degree:", deg_elementwise)
    println("Elementwise degree q:", deg_elementwise_q)
    println("Target variable:", target_var)

    filename = """ + json.dumps(filename) + """
    data_name = """ + json.dumps(data_name) + """

    data_bin_star = DelimitedFiles.readdlm(data_name)

    function run_model(N_data_points)
        candidate_mons = [
            mons_of_max_degree_and_unit_overall(vars, deg, deg_overall, deg_elementwise)
            for ai=axioms
        ]
        @show size.(candidate_mons)

        model = Model(Gurobi.Optimizer)

        mons_q = mons_of_max_degree_and_unit_overall(measured_vars, deg, deg_overall_q, deg_elementwise_q)
        coeff_q =   @variable(model, [1:size(mons_q,1)], base_name="q")
        abs_coeff_q =   @variable(model, [1:size(mons_q,1)], base_name="abq")
        
        q = sum(ci .* mi for (ci, mi)=zip(coeff_q, mons_q)) 
        coeff_αs = [
            @variable(model, [1:size(X,1)], base_name="α$i")
            for (i,X)=enumerate(candidate_mons)
        ]
        
        @show size.(coeff_αs)
        αs = [sum(ci .* mi) for (ci, mi)=zip(coeff_αs, candidate_mons)]

        residual = q - sum(αᵢ * aᵢ for (αᵢ, aᵢ)=zip(αs,axioms));
        @show length(residual)
        eqs = coefficients(residual)
        @show size(eqs)

        # Constraints on q
        @constraint model abs_coeff_q.>=coeff_q
        @constraint model abs_coeff_q.>=-coeff_q

        # Constraints on the residual
        @variable(model, abseqs[i=1:size(eqs,1)])
        @constraint model abseqs.>=eqs
        @constraint model abseqs.>=-eqs

        # Constraints on fit to data
        @variable(model, t[i=1:N_data_points]>=0.0) 
        @constraint(model, imposeabs1[i=1:N_data_points], t[i]>=q(data_bin_star[i,:]))
        @constraint(model, imposeabs2[i=1:N_data_points], t[i]>=-q(data_bin_star[i,:]))
        
        # Constraints on the target variable
        @constraint model sum(coeff_q[degree.(mons_q, target_var).>0]) == 8.0

        # Constraints on sparsity of q
        # @variable(model, z[i=1:size(abs_coeff_q,1)], Bin)
        # @constraint(model, imposeLogical[i=1:size(abs_coeff_q,1)], !(z[i]) .=> {abs_coeff_q[i].==0.0}) #logical 
        # @constraint(model, sum(z)<=3)

        # Objective
        #weight_1=5
        #@objective model Min (weight_1/sqrt(N_data_points))*sum(t) + (1-weight_1/sqrt(N_data_points))*sum(abs_coeff_q)

        weight_1=0.9
        weight_2=0.01
        @objective model Min (weight_1/sqrt(N_data_points))*sum(t) + weight_2*sum(abs_coeff_q) + (1.0-weight_1-weight_2)*sum(abseqs)

        optimize!(model)
        @show termination_status(model)

        ######### Post-processing the model results and saving data  #########

        ############# RESIDUAL #####################################

        # Show the optimal residual
        value_residual = p -> sum(value.(coefficients(p)).* monomials(p))
        value_res = value_residual(residual)
        @show value_res

        # Extracting coefficients and monomials
        coeffs_res = value.(coefficients(value_res))
        mons_res = monomials(value_res)

        threshold = 0.05 # threshold for filtering monomials for coefficients of terms

        # Filter for coefficients
        filtered_coeffs_res = coeffs_res[abs.(coeffs_res) .> threshold]
        filtered_mons_res = mons_res[abs.(coeffs_res) .> threshold]

        value_res_filtered = sum(filtered_coeffs_res .* filtered_mons_res)
        @show value_res_filtered

        # Round the coefficients to n_digits-th decimal place
        n_digits = 2
        rounded_coeffs_res = round.(filtered_coeffs_res, digits=n_digits)
        value_res_rounded = sum(rounded_coeffs_res .* filtered_mons_res)
        @show value_res_rounded

        ############# Alpha POLYNOMIALS #####################################
        # Finding the output polynomial
        value_poly = p -> sum(value.(coefficients(p)).* monomials(p))
        value_αs = value_poly.(αs)
        @show value_αs

        # Extracting coefficients and monomials
        coeffs_αs = [value.(coefficients(v)) for v in value_αs]
        mons_αs = [monomials(v) for v in value_αs]

        threshold = 0.0005 # threshold for filtering monomials for coefficients of terms

        # Filter for coefficients
        filtered_coeffs_αs = [coeffs_αs[i][abs.(coeffs_αs[i]) .> threshold] for i in 1:size(coeffs_αs,1)]
        filtered_mons_αs = [mons_αs[i][abs.(coeffs_αs[i]) .> threshold] for i in 1:size(coeffs_αs,1)]

        value_αs_filtered = [sum(filtered_coeffs_αs[i] .* filtered_mons_αs[i]) for i in 1:size(filtered_coeffs_αs,1)]
        @show value_αs_filtered

        # Round the coefficients to n_digits-th decimal place
        n_digits = 2
        rounded_coeffs_αs = [round.(filtered_coeffs_αs[i], digits=n_digits) for i in 1:size(filtered_coeffs_αs,1)]
        value_αs_rounded = [sum(rounded_coeffs_αs[i] .* filtered_mons_αs[i]) for i in 1:size(filtered_coeffs_αs,1)]

        @show value_αs_rounded



        ############# Q POLYNOMIAL #####################################
        # Finding the output polynomial
        value_poly2 = p -> sum(value.(coefficients(p)).* monomials(p))
        value_q = value_poly2(q)
        @show value_q
        
        # Extracting coefficients and monomials
        coeffs = value.(coefficients(value_q))
        mons = monomials(value_q)

        threshold = 0.05 # threshold for filtering monomials for coefficients of terms

        # Fileter for coefficients
        filtered_coeffs = coeffs[abs.(coeffs) .> threshold]
        filtered_mons = mons[abs.(coeffs) .> threshold]

        value_q_filtered = sum(filtered_coeffs .* filtered_mons)
        @show value_q_filtered

        # Round the coefficients to n_digits-th decimal place
        n_digits = 2
        rounded_coeffs = round.(filtered_coeffs, digits=n_digits)
        value_q_rounded = sum(rounded_coeffs .* filtered_mons)
        @show value_q_rounded

        save_polynomial_to_file(value_q, filename, "q polynomial")
        save_polynomial_to_file(value_q_filtered, filename, "filtered polynomial")
        save_polynomial_to_file(value_q_rounded, filename, "rounded polynomial")

        distance = min(sum(abs.(coefficients(value_q-target_poly))), sum(abs.(coefficients(value_q+target_poly))))
        @show N_data_points
        @show distance

        

    end

    N_data_points = 100
    run_model(N_data_points)


    """


    # Create a temporary file to store the Julia code
    with tempfile.NamedTemporaryFile(mode='w', suffix='.jl', delete=False) as temp_file:
        temp_file.write(julia_code)
        temp_file_path = temp_file.name

    try:
        # Construct the Julia command
        julia_command = [
            "julia",
            "--project=" + env_path,
            "-e",
            f'using Pkg; Pkg.activate("{env_path}"); include("{temp_file_path}")'
        ]

        # Run the Julia command
        result = subprocess.run(julia_command, check=True, capture_output=True, text=True)
        print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error occurred:")
        print(e.stdout)
        print(e.stderr)
    finally:
        os.unlink(temp_file_path)






def ai_hilbert_axioms(env_path, variables, axioms, measured_variables, target_var, deg, deg_overall, deg_overall_q, deg_elementwise, deg_elementwise_q, target, filename):

    julia_code = """
    using DynamicPolynomials
    using JuMP
    using MosekTools
    using LinearAlgebra
    using CSV
    using DataFrames
    import JSON
    using Dates
    using DelimitedFiles
    using Plots
    using Statistics
    using Gurobi


    function all_monomials_up_to_max_deg(x, deg)
        if size(x,1) == 0
            [1]
        else
        [ x[1]^k * m for k=0:deg 
                for m=all_monomials_up_to_max_deg(x[2:end], deg)
        ]
        end
    end

    function mons_of_max_degree_and_unit(x, deg)
        [m
            for m=all_monomials_up_to_max_deg(x, deg)
        ]
    end

    function all_monomials_up_to_max_deg_overall(x, deg, deg_overall, deg_elem)
        if size(x,1) == 0
            [1]
        else
        [ x[1]^k * m for k=0:min(deg, deg_overall, deg_elem[1]) 
                    for m=all_monomials_up_to_max_deg_overall(x[2:end], deg, deg_overall-k, deg_elem[2:end])
        ]
        end
    end

    function mons_of_max_degree_and_unit_overall(x, deg, deg_overall, deg_elem)
        [m
            for m=all_monomials_up_to_max_deg_overall(x, deg, deg_overall, deg_elem)
        ]
    end
                    
    function degree_poly(p)
        maximum(degree.(monomials(p)))
    end

    function save_polynomial_to_file(poly, filename, description)
        open(filename*"_poly.txt", "a") do file
            write(file, string(description, ": ", poly) * "\n") 
        end
    end


    function create_variables(var_strings)
        variables = Symbol[]
        for var_string in var_strings
            # Remove any whitespace and non-alphanumeric characters (except underscore)
            cleaned_var = replace(strip(var_string), r"[^a-zA-Z0-9_]" => "")
            
            # Ensure the variable doesn't start with a number
            if !isempty(cleaned_var) && !isdigit(cleaned_var[1])
                push!(variables, Symbol(cleaned_var))
            end
        end
        return variables
    end
    
    variables = """ + json.dumps(variables) + """
    axiom_strings = """ + json.dumps(axioms) + """
    measured_variables = """ + json.dumps(measured_variables) + """
    target_variable = """ + json.dumps(target_var) + """
    target = """ + json.dumps(target) + """

    eval(:((@polyvar $(Symbol.(variables)...))))
    var_dict = Dict(var => eval(Symbol(var)) for var in variables)
    vars = [eval(Symbol(var)) for var in variables]

    function replace_vars(expr)
        for (var, poly) in var_dict
            expr = replace(expr, var => string(poly))
        end
        return expr
    end

    axioms = [eval(Meta.parse(replace_vars(axiom))) for axiom in axiom_strings]
    
    measured_vars = [eval(Symbol(var)) for var in measured_variables]

    target_poly = eval(Meta.parse(replace_vars(target)))

    deg_elementwise = """ + json.dumps(deg_elementwise) + """
    deg_elementwise_q = """ + json.dumps(deg_elementwise_q) + """
    deg = """ + str(deg) + """
    deg_overall = """ + str(deg_overall) + """
    deg_overall_q = """ + str(deg_overall_q) + """

    target_var = eval(Symbol(target_variable))

    println("Variables:", vars)
    println("Axioms:", axioms)
    println("Measured variables:", measured_vars)
    println("Degree:", deg)
    println("Overall degree:", deg_overall)
    println("Overall degree q:", deg_overall_q)
    println("Elementwise degree:", deg_elementwise)
    println("Elementwise degree q:", deg_elementwise_q)
    println("Target variable:", target_var)

    filename = """ + json.dumps(filename) + """

    candidate_mons = [
        mons_of_max_degree_and_unit_overall(vars, deg, deg_overall, deg_elementwise)
        for ai=axioms
    ]
    @show size.(candidate_mons)

    model = Model(Gurobi.Optimizer)

    mons_q = mons_of_max_degree_and_unit_overall(measured_vars, deg, deg_overall_q, deg_elementwise_q)
    coeff_q =   @variable(model, [1:size(mons_q,1)], base_name="q")

    q = sum(ci .* mi for (ci, mi)=zip(coeff_q, mons_q)) 
    coeff_αs = [
        @variable(model, [1:size(X,1)], base_name="α$i")
        for (i,X)=enumerate(candidate_mons)
    ]
    
    @show size.(coeff_αs)
    αs = [sum(ci .* mi) for (ci, mi)=zip(coeff_αs, candidate_mons)]

    residual = q - sum(αᵢ * aᵢ for (αᵢ, aᵢ)=zip(αs,axioms));
    @show length(residual)
    eqs = coefficients(residual)
    @show size(eqs)

    @constraint model eqs .== 0
    @constraint model sum(coeff_q[degree.(mons_q, target_var).>0]) == 1.0
    @objective model Max 0

    optimize!(model)
    @show termination_status(model)

    ######### Post-processing the model results and saving data  #########

    ############# RESIDUAL #####################################

    # Show the optimal residual
    value_residual = p -> sum(value.(coefficients(p)).* monomials(p))
    value_res = value_residual(residual)
    @show value_res

    # Extracting coefficients and monomials
    coeffs_res = value.(coefficients(value_res))
    mons_res = monomials(value_res)

    #@show eval_values_res

    threshold = 0.05 # threshold for filtering monomials for coefficients of terms

    # Filter for coefficients
    filtered_coeffs_res = coeffs_res[abs.(coeffs_res) .> threshold]
    filtered_mons_res = mons_res[abs.(coeffs_res) .> threshold]

    value_res_filtered = sum(filtered_coeffs_res .* filtered_mons_res)
    @show value_res_filtered

    # Round the coefficients to n_digits-th decimal place
    n_digits = 2
    rounded_coeffs_res = round.(filtered_coeffs_res, digits=n_digits)
    value_res_rounded = sum(rounded_coeffs_res .* filtered_mons_res)
    @show value_res_rounded

    ############# Alpha POLYNOMIALS #####################################
    # Finding the output polynomial
    value_poly = p -> sum(value.(coefficients(p)).* monomials(p))
    value_αs = value_poly.(αs)
    @show value_αs

    # Extracting coefficients and monomials
    coeffs_αs = [value.(coefficients(v)) for v in value_αs]
    mons_αs = [monomials(v) for v in value_αs]

    threshold = 0.05 # threshold for filtering monomials for coefficients of terms

    # Filter for coefficients
    filtered_coeffs_αs = [coeffs_αs[i][abs.(coeffs_αs[i]) .> threshold] for i in 1:size(coeffs_αs,1)]
    filtered_mons_αs = [mons_αs[i][abs.(coeffs_αs[i]) .> threshold] for i in 1:size(coeffs_αs,1)]

    value_αs_filtered = [sum(filtered_coeffs_αs[i] .* filtered_mons_αs[i]) for i in 1:size(filtered_coeffs_αs,1)]
    @show value_αs_filtered

    # Round the coefficients to n_digits-th decimal place
    n_digits = 2
    rounded_coeffs_αs = [round.(filtered_coeffs_αs[i], digits=n_digits) for i in 1:size(filtered_coeffs_αs,1)]
    value_αs_rounded = [sum(rounded_coeffs_αs[i] .* filtered_mons_αs[i]) for i in 1:size(filtered_coeffs_αs,1)]

    @show value_αs_rounded



    ############# Q POLYNOMIAL #####################################
    # Finding the output polynomial
    value_poly2 = p -> sum(value.(coefficients(p)).* monomials(p))
    value_q = value_poly2(q)
    @show value_q
    
    # Extracting coefficients and monomials
    coeffs = value.(coefficients(value_q))
    mons = monomials(value_q)

    threshold = 0.05 # threshold for filtering monomials for coefficients of terms

    # Fileter for coefficients
    filtered_coeffs = coeffs[abs.(coeffs) .> threshold]
    filtered_mons = mons[abs.(coeffs) .> threshold]

    value_q_filtered = sum(filtered_coeffs .* filtered_mons)
    @show value_q_filtered

    # Round the coefficients to n_digits-th decimal place
    n_digits = 2
    rounded_coeffs = round.(filtered_coeffs, digits=n_digits)
    value_q_rounded = sum(rounded_coeffs .* filtered_mons)
    @show value_q_rounded

    distance = min(sum(abs.(coefficients(value_q-target_poly))), sum(abs.(coefficients(value_q+target_poly))))
    @show distance

    save_polynomial_to_file(value_q, filename, "q polynomial")
    save_polynomial_to_file(value_q_filtered, filename, "filtered polynomial")
    save_polynomial_to_file(value_q_rounded, filename, "rounded polynomial")


    """


    # Create a temporary file to store the Julia code
    with tempfile.NamedTemporaryFile(mode='w', suffix='.jl', delete=False) as temp_file:
        temp_file.write(julia_code)
        temp_file_path = temp_file.name

    try:
        # Construct the Julia command
        julia_command = [
            "julia",
            "--project=" + env_path,
            "-e",
            f'using Pkg; Pkg.activate("{env_path}"); include("{temp_file_path}")'
        ]

        # Run the Julia command
        result = subprocess.run(julia_command, check=True, capture_output=True, text=True)
        print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error occurred:")
        print(e.stdout)
        print(e.stderr)
    finally:
        os.unlink(temp_file_path)