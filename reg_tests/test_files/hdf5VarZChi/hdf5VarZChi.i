Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

realms:

  - name: realm_1
    mesh: 2cm_ped_35K_mks.g
    use_edges: yes

    time_step_control:
     target_courant: 2.0
     time_step_change_factor: 1.2
   
    equation_systems:
      name: theEqSys
      max_iterations: 1 

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        mixture_fraction: solve_scalar

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-3

        - MixtureFraction:
            name: myZ
            max_iterations: 1
            convergence_tolerance: 1.e-3
        - MixtureFraction:
            name: myZ
            max_iterations: 1
            convergence_tolerance: 1.e-3

    material_properties:

      target_name: block_1

      table_file_name: SLFM_CGauss_C2H4_ZMean_ZScaledVarianceMean_logChiMean.h5

      specifications:
 
        - name: density
          type: hdf5table
          independent_variable_set: [mixture_fraction, scalar_variance, scalar_dissipation]
          table_name_for_property: density
          table_name_for_independent_variable_set: [ZMean, ZScaledVarianceMean, ChiMean]
          aux_variables: temperature
          table_name_for_aux_variables: temperature

        - name: viscosity
          type: hdf5table
          independent_variable_set: [mixture_fraction, scalar_variance, scalar_dissipation]
          table_name_for_property: mu
          table_name_for_independent_variable_set: [ZMean, ZScaledVarianceMean, ChiMean]

    initial_conditions:
      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0
          velocity: [0,0]  
          mixture_fraction: 0.0
  
    boundary_conditions:

    - wall_boundary_condition: bc_floor
      target_name: surface_1
      wall_user_data:
        velocity: [0,0,0]

    - inflow_boundary_condition: bc_inflow
      target_name: surface_2
      inflow_user_data:
        velocity: [0,0,15.0]
        mixture_fraction: 1.0

    - wall_boundary_condition: bc_pipe
      target_name: surface_3
      wall_user_data:
        velocity: [0,0,0]

    - open_boundary_condition: bc_side
      target_name: surface_4
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        mixture_fraction: 0.0

    - open_boundary_condition: bc_top
      target_name: surface_5
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        mixture_fraction: 0.0

    solution_options:
      name: myOptions
      turbulence_model: wale  

      options:
        - hybrid_factor:
            velocity: 1.0
            mixture_fraction: 1.0

        - upw_factor:
            velocity: 1.0
            mixture_fraction: 1.0

        - laminar_schmidt:
            mixture_fraction: 0.90 

        - turbulent_schmidt:
            mixture_fraction: 0.9

        - source_terms:
           continuity: density_time_derivative

        - limiter:
            pressure: no
            velocity: yes 
            mixture_fraction: yes 

        - projected_nodal_gradient:
            pressure: element
            velocity: edge
            mixture_fraction: edge

    output:
      output_data_base_name: hdf5VarZChi.e
      output_frequency: 10
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - mixture_fraction
       - scalar_dissipation
       - density
       - viscosity
       - temperature

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 50
      time_step: 1.0e-5
      time_stepping_type: automatic
      time_step_count: 0
      second_order_accuracy: no 

      realms:
        - realm_1
