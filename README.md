# glioma_ddri
.
├── LICENSE.md
├── Manifest.toml
├── Project.toml
├── README.md
├── assets
├── model
│   ├── gdc_dosing2.jl
│   ├── gdc_params.jl
│   ├── gdc_pkpd2.jl
│   ├── rg_dosing2.jl
│   ├── rg_params.jl
│   └── rg_pkpd2.jl
├── results
│   ├── correlation_heatmap
│   ├── dose_matrix_heatmap
│   ├── outputs_violin
│   ├── population_spaghetti
│   └── sensitivity_heatmap
├── scripts
│   ├── analysis
│   │   ├── correlation_heatmap.jl
│   │   ├── dose_viz_vp4.jl
│   │   ├── efast_histogram_alldose.jl
│   │   ├── kaplan_meier.jl
│   │   ├── outputs_violin_drugconditions.jl
│   │   ├── population_tumor_trajectories.jl
│   │   └── simulate_save_powerbi.jl
│   ├── cluster_commands.txt
│   ├── gdc_optimization.jl
│   ├── rg_optimization.jl
│   └── setup
│       ├── bounds4.jl
│       ├── generate_vp4.jl
│       └── scaling_v4.jl
└── utilities
    └── utils.jl

] activate .
add Plots Random
add DifferentialEquations ModelingToolkit
add Optim OptimizationOptimJL OptimizationOptimisers
add Zygote Enzyme Flux ForwardDiff GlobalSensitivity SciMLSensitivity

add Plots Random DifferentialEquations ModelingToolkit Optim OptimizationOptimJL OptimizationOptimisers Zygote Enzyme Flux ForwardDiff GlobalSensitivity SciMLSensitivity

Assumes that the order of parameters in the model equations is the same in pkpd!.jl and params.jl and that the pk_params.jl assets has ordered keys with the same order in params.jl and in pk!.jl. These are required in both compute_dose.jl and generate_vp.jl