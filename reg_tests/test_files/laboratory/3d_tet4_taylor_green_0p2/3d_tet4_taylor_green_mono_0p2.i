Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1
    error_estimator: errest_1

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
    recompute_preconditioner: no    
    muelu_xml_file_name: ../../../xml/matches_ml_default.xml

realms:

  - name: realm_1
    mesh: mesh/3d_tet4_taylor_green_0p2.exo
    use_edges: no
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 2
   
      solver_system_specification:
        uvwp: solve_scalar

      systems:
        - LowMachMonoEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-2

    initial_conditions:

      - user_function: icUser
        target_name: block_1
        user_function_name:
         velocity: TaylorGreen
         pressure: TaylorGreen 
         
    material_properties:
      target_name: block_1

      specifications:

        - name: density
          type: constant
          value: 1.0

        - name: viscosity
          type: constant
          value: 625.0e-6

    boundary_conditions:

    - periodic_boundary_condition: bc_left_right
      target_name: [surface_1, surface_2]
      periodic_user_data:
        search_tolerance: 0.0001 

    - periodic_boundary_condition: front_back
      target_name: [surface_3, surface_4]
      periodic_user_data:
        search_tolerance: 0.0001 

    - periodic_boundary_condition: top_bot
      target_name: [surface_5, surface_6]
      periodic_user_data:
        search_tolerance: 0.0001 

    solution_options:
      name: myOptions
      turbulence_model: laminar

      use_consolidated_solver_algorithm: yes
  
      options:

        - element_source_terms:
            uvwp: [uvwp_lumped_time_advection_diffusion]

    turbulence_averaging:
      time_filter_interval: 10.0
      specifications:
        - name: one
          target_name: block_1
          reynolds_averaged_variables:
            - velocity

          compute_q_criterion: yes
          compute_vorticity: yes
          compute_lambda_ci: yes
          compute_mean_resolved_ke: yes

    output:
      output_data_base_name: output/3d_tet4_taylor_green_mono_0p2.e
      output_frequency: 20
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - vorticity
       - q_criterion
       - dual_nodal_volume

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 1.0
      time_step: 0.04 
      time_stepping_type: fixed 
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
