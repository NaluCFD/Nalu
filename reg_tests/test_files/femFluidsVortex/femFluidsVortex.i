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
    method: cg
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0
    recompute_preconditioner: no
    muelu_xml_file_name: ../../xml/milestone.xml

realms:

  - name: realm_1
    mesh: ../../mesh/hex8_20_10_1.g
    activate_fem: yes
    use_edges: no
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 3
   
      solver_system_specification:
        pressure: solve_cont
        velocity: solve_scalar
        dpdx: solve_scalar
       
      systems:
        - LowMachFemEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-2

    initial_conditions:

      - user_function: icUser
        target_name: [block_1]
        user_function_name:
         velocity: wind_energy_taylor_vortex
        user_function_parameters:
         velocity: [0.0,0.0,0.5,1.0,1.0] 
         
    material_properties:
      target_name: block_1

      specifications:

        - name: density
          type: constant
          value: 1.0

        - name: viscosity
          type: constant
          value: 1.8e-5

    boundary_conditions:

    - periodic_boundary_condition: bc_front_back
      target_name: [surface_1, surface_2]
      periodic_user_data:
        search_tolerance: 0.0001 

    - periodic_boundary_condition: left_right
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
            momentum: [momentum_time_derivative, advection, diffusion]
            continuity: [advection]
            dpdx: [interior_png]

        - consistent_mass_matrix_png:
            pressure: yes

    output:
      output_data_base_name: output/fem_hex8_cmm.e
      output_frequency: 5
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - vorticity
       - q_criterion

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 25
      time_step: 0.05
      time_stepping_type: adaptive 
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
