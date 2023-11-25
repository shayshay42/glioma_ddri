using Plots, DifferentialEquations, MAT, Statistics, Serialization, ModelingToolkit
using HDF5
pyplot()

include("../../utilities/utils.jl")

for drug in ["gdc"] #["rg", "gdc"]

    include("../../models/$(drug)_pkpd2.jl")
    include("../../models/$(drug)_dosing2.jl")

    population = deserialize("./assets/$(drug)_vp5.jls")
    # min_tumor = deserialize("./assets/$(drug)_min_tumor_max_dose.jls")

    #deserialize optimum .jls from results
    # optimum = deserialize("./results/optimum/optima_$(drug)_lfbgs_auc_oct27_dms.jls")
    # #mkae it into a Float64 array
    # opt = Array{Float64, 2}(optimum[1:end-2,:])
    #save it in the same directory as a HDF5
    # h5write("./results/optimum/optima_$(drug)_lfbgs_auc_oct24_dms.h5", "optimum", opt)
    dosetimes = eval(Symbol("$(drug)_dosetimes"))
    for (cond, mask) in [("TMZ_treat"=>[1,0,0]),("TMZ_treat&adjuv"=>[1,1,0]),("TMZ_DDRi"=>[1,1,1]),("DDRi_only"=>[0,0,1])]

        # if cond != "DDRi_only"
        #     continue
        # end

        include("../../models/$(drug)_params.jl")
        println("This is the new conditon: ",cond)
        tmz_treat_dose = mask[1]*75.0*avg_human_surface_area
        tmz_adjuv_dose = mask[2]*150.0*avg_human_surface_area
        doses = ones(length(dosetimes)).*(mask[3]*dose_amount)

        u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        tspan = (0,end_time+7).*hours
        ode_params = [ode_params; 1.0; 0.1]
        p = [ode_params;doses]
        prob = ODEProblem(pk_pd!,u0,tspan,p)
        sys = modelingtoolkitize(prob)
        sys = structural_simplify(sys)
        prob_jac = ODEProblem(sys, u0, tspan, p, jac=true)
        sol = solve(prob_jac, callback=hit, saveat=0.1)
        # display(plotter(sol))

        # nb_patients = size(vp, 1)
        nb_patients = size(population, 2)

        c_solutions = Array{Float64, 2}(undef,  nb_patients, size(sol.t)[1])
        d1_solutions = Array{Float64, 2}(undef,  nb_patients, size(sol.t)[1])
        d2_solutions = Array{Float64, 2}(undef,  nb_patients, size(sol.t)[1])

        Threads.@threads for i in 1:nb_patients
            println("\rPatient: $i")
            p = [population[:, i];1.0;0.1;doses]

            p_prob = remake(prob_jac, p=p)
            p_sol = solve(p_prob, callback=hit, saveat=0.1)#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
            
            sols = Array(p_sol)

            c_solutions[i,:]  = sols[1,:]
            d1_solutions[i,:] = sols[5,:]
            d2_solutions[i,:] = sols[7,:]
        end

        if cond == "DDRi_only"
            data_matrix = d2_solutions
            timepoints = sol.t ./24

            title = "$drug Plasma Presence for Virtual Population $cond"
            means = vec(mean(data_matrix, dims=1))
        
            # Define pastel colors for the mean ribbon
            mean_line_color = :pink#RGB(0.4, 0.76, 0.64)  # A pastel green for the mean line
            ribbon_fill_color = :lightblue#RGB(0.85, 0.9, 0.85)  # A lighter pastel green for the ribbon fill
        
            # Define a neutral color for the population trajectories
            population_trajectory_color = :lightgray#RGB(0.8, 0.8, 0.8)  # Light grey
                    # Plot the population trajectories in light grey
            plot(timepoints, data_matrix'[:,1], alpha=0.3, linecolor=population_trajectory_color, grid=false, label="Drug Trajectories", legend=:topright)
            plot!(timepoints, data_matrix'[:,2:end], alpha=0.3, linecolor=population_trajectory_color, grid=false, label=false)
        
            # Calculate the standard deviation and ensure it does not go below 0
            std_dev = vec(std(data_matrix, dims=1))
            # lower_ribbon = max.(means - std_dev, 0)  # Ensure ribbon does not go negative
        
            # Plot the mean and standard deviation with the pastel colors
            plot!(timepoints, means, 
                ribbon=(std_dev), 
                fillalpha=0.5, linewidth=3, linecolor=mean_line_color, fillcolor=ribbon_fill_color, 
                label="Mean ± SD")
        
            # Set the title and axis labels
            title!(title)
            xlabel!("Time (hours)")
            ylabel!("Drug Amount (mg)")
        
            savefig("./results/pop_plots/$drug/$(drug)_model_alldose_$cond.png")
            savefig("./results/pop_plots/$drug/$(drug)_model_alldose_$cond.svg")
        end

        data_matrix = c_solutions
        timepoints = sol.t ./24

        title = "Tumor Trajectories of Virtual Population \n Treated with $cond in $(drug) Model"
        means = vec(mean(data_matrix, dims=1))

                # Define pastel colors for the mean ribbon
        mean_line_color = :pink#RGB(0.4, 0.76, 0.64)  # A pastel green for the mean line
        ribbon_fill_color = :lightblue#RGB(0.85, 0.9, 0.85)  # A lighter pastel green for the ribbon fill

        # Define a neutral color for the population trajectories
        population_trajectory_color = :lightgray#RGB(0.8, 0.8, 0.8)  # Light grey
                # Plot the population trajectories in light grey
        plot(timepoints, data_matrix'[:,1], alpha=0.3, linecolor=population_trajectory_color, grid=false, label="Tumor Trajectories", legend=:topleft, ylim=(0,40))
        plot!(timepoints, data_matrix'[:,2:end], alpha=0.3, linecolor=population_trajectory_color, grid=false, label=false, ylim=(0,40))

        # Calculate the standard deviation and ensure it does not go below 0
        std_dev = vec(std(data_matrix, dims=1))
        # lower_ribbon = max.(means - std_dev, 0)  # Ensure ribbon does not go negative

        # Plot the mean and standard deviation with the pastel colors
        plot!(timepoints, means, 
            ribbon=(std_dev), 
            fillalpha=0.5, linewidth=3, linecolor=mean_line_color, fillcolor=ribbon_fill_color, 
            label="Mean ± SD", ylim=(0,40))

        # Set the title and axis labels
        title!(title)
        xlabel!("Time (days)")
        ylabel!("Tumor Volume (ml)")

        savefig("./results/pop_plots/$(drug)/$(drug)_model_$(cond)_tumorvolume.png")
        savefig("./results/pop_plots/$(drug)/$(drug)_model_$(cond)_tumorvolume.svg")
    end
end

include("../../models/rg_dosing2.jl")
for (cond, mask) in [("TMZ_treat"=>tmz_treat_dose),("TMZ_adjuv"=>tmz_adjuv_dose),("rg"=>1800.0*avg_human_surface_area),("gdc"=>600.0)]

    if cond != "TMZ_treat" || cond != "TMZ_adjuv"
        drug="rg"
        state = 3
        trajectory_var = "d1_solutions"
    else
        drug=cond
        state = 6
        trajectory_var = "d2_solutions"
    end

    include("../../models/$(drug)_pkpd2.jl")
    include("../../models/$(drug)_dosing2.jl")
    include("../../models/$(drug)_params.jl")

    population = deserialize("./assets/$(drug)_vp4.jls")
    # min_tumor = deserialize("./assets/$(drug)_min_tumor_max_dose.jls")

    # #deserialize optimum .jls from results
    # optimum = deserialize("./results/optimum/optima_$(drug)_lfbgs_auc_oct27_dms.jls")
    # #mkae it into a Float64 array
    # opt = Array{Float64, 2}(optimum[1:end-2,:])


    u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    tspan = (0.92,2.1).*hours
    ode_params = [ode_params; 1.0; 0.1]
    p = ode_params
    prob = ODEProblem(pk_pd!,u0,tspan,p)
    sys = modelingtoolkitize(prob)
    sys = structural_simplify(sys)
    prob_jac = ODEProblem(sys, u0, tspan, p, jac=true)

    condition(u, t, integrator) = t == 24.0
    affect!(integrator) = integrator.u[state] += mask
    cb = DiscreteCallback(condition, affect!)
    sol = solve(prob_jac, callback=cb, saveat=0.1,tstops=[24.0])
    # display(plotter(sol))

    # nb_patients = size(vp, 1)
    nb_patients = size(population, 2)

    c_solutions = Array{Float64, 2}(undef,  nb_patients, size(sol.t)[1])
    d1_solutions = Array{Float64, 2}(undef,  nb_patients, size(sol.t)[1])
    d2_solutions = Array{Float64, 2}(undef,  nb_patients, size(sol.t)[1])

    Threads.@threads for i in 1:nb_patients
        println("\rPatient: $i")
        p = [population[:, i];1.0;0.1]

        p_prob = remake(prob_jac, p=p)
        p_sol = solve(p_prob, callback=cb, saveat=0.1,tstops=[24.0])#, abstol=1e-10, reltol=1e-10,dtmax=1)#, alg_hints=[:stiff])
        
        sols = Array(p_sol)

        c_solutions[i,:]  = sols[1,:]
        d1_solutions[i,:] = sols[5,:]
        d2_solutions[i,:] = sols[7,:]
    end

    if cond != "TMZ_treat" || cond != "TMZ_adjuv"
        data_matrix = d1_solutions
    else
        data_matrix = d2_solutions
    end
    timepoints = sol.t

    title = "$cond Plasma Presence for Virtual Population"
    means = vec(mean(data_matrix, dims=1))

    # Define pastel colors for the mean ribbon
    mean_line_color = :pink#RGB(0.4, 0.76, 0.64)  # A pastel green for the mean line
    ribbon_fill_color = :lightblue#RGB(0.85, 0.9, 0.85)  # A lighter pastel green for the ribbon fill

    # Define a neutral color for the population trajectories
    population_trajectory_color = :lightgray#RGB(0.8, 0.8, 0.8)  # Light grey
            # Plot the population trajectories in light grey
    plot(timepoints, data_matrix'[:,1], alpha=0.3, linecolor=population_trajectory_color, grid=false, label="Drug Trajectories", legend=:topright)
    plot!(timepoints, data_matrix'[:,2:end], alpha=0.3, linecolor=population_trajectory_color, grid=false, label=false)

    # Calculate the standard deviation and ensure it does not go below 0
    std_dev = vec(std(data_matrix, dims=1))
    # lower_ribbon = max.(means - std_dev, 0)  # Ensure ribbon does not go negative

    # Plot the mean and standard deviation with the pastel colors
    plot!(timepoints, means, 
        ribbon=(std_dev), 
        fillalpha=0.5, linewidth=3, linecolor=mean_line_color, fillcolor=ribbon_fill_color, 
        label="Mean ± SD")

    # Set the title and axis labels
    title!(title)
    xlabel!("Time (hours)")
    ylabel!("Drug Amount (mg)")

    savefig("./results/pop_plots/$(cond)_$(drug)_model_singledose.png")
    savefig("./results/pop_plots/$(cond)_$(drug)_model_singledose.svg")
end

