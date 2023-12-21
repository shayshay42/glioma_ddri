using Pkg, ArgParse

function main(drug="rg", num_patients=1000, seed=123)
    # vp_filename = "$(drug)_$(num_patients)_vp.jls"
    # Activate and install project dependencies
    Pkg.activate(".")
    Pkg.instantiate()

    include("../../model/$(drug)_pkpd2.jl")
    include("../../model/$(drug)_dosing2.jl")
    include("../../model/$(drug)_params.jl")
    include("../../scripts/setup/init_integrate.jl")

    # Step 0: Verify the existence of model, dosing, and params files
    verify_required_files([
        "model/$(drug)_dosing2.jl",
        "model/$(drug)_params.jl",
        "model/$(drug)_pkpd2.jl"
    ])

    # # Step 1: Create bounds on the varying parameters
    # include("scripts/setup/bounds4.jl")




# USE THE STRUCT FROM THE BEGINNING




    # Step 2: Create virtual population
    # run(`julia scripts/setup/generate_vp4.jl --drug $drug --nbr $num_patients --filename $vp_filename --seed $seed`)

    # Step 3: Compute the lower and upper bound of tumor AUC for each patient
    # run(`julia scripts/setup/scaling_v4.jl --drug $drug --filename $vp_filename`)

    include("./setup.jl")
    patients = generate_patients_struct(num_patients, seed, drug)

    # Step 4: Compute the optimal doses
    

    # # Steps 4-12: Additional pipeline steps
    # # Replace `--args` with actual arguments as needed
    # run(`julia scripts/analysis/population_tumor_trajectories.jl --args`)
    # run(`julia scripts/gdc_optimization.jl --args`)
    # run(`julia scripts/analysis/simulate_save_powerbi.jl --args`)
    # run(`julia scripts/analysis/efast_histogram_alldose.jl --args`)
    # run(`julia scripts/analysis/outputs_violin_drugconditions.jl --args`)
    # run(`julia scripts/analysis/correlation_heatmap.jl --args`)
    # run(`julia scripts/analysis/dose_viz_vp4.jl --args`)
    # run(`julia scripts/analysis/kaplan_meier.jl --args`)

    # # Additional steps can be added in the same way
end

function verify_required_files(file_paths::Array{String,1})
    for path in file_paths
        if !isfile(path)
            error("Required file not found: $path")
        end
    end
    println("All required files verified.")
end

main()

