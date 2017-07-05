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
    tolerance: 1e-6
    max_iterations: 300
    kspace: 300
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres 
    preconditioner: muelu
    tolerance: 1e-6
    max_iterations: 300
    kspace: 300
    output_level: 0
    recompute_preconditioner: no
    muelu_xml_file_name: ../../xml/muelu_STV_HO.xml

realms:

  - name: realm_1
    mesh: ../../mesh/tquad4_16.g
    use_edges: no
    polynomial_order: 4
    automatic_decomposition_type: rcb     

    equation_systems:
      name: theEqSys
      max_iterations: 1
   
      solver_system_specification:
        pressure: solve_cont
        velocity: solve_scalar
        dpdx: solve_scalar

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-2

    initial_conditions:

      - user_function: icUser
        target_name: block_1
        user_function_name:
         velocity: SteadyTaylorVortex
         pressure: SteadyTaylorVortex 
         
    material_properties:
      target_name: block_1

      specifications:

        - name: density
          type: constant
          value: 1.0

        - name: viscosity
          type: constant
          value: 0.001

    boundary_conditions:

    - inflow_boundary_condition: bc_left
      target_name: surface_1
      inflow_user_data:
        user_function_name:
          velocity: SteadyTaylorVortex

    - inflow_boundary_condition: bc_right
      target_name: surface_2
      inflow_user_data:
        user_function_name:
          velocity: SteadyTaylorVortex

    - inflow_boundary_condition: bc_top
      target_name: surface_3
      inflow_user_data:
        user_function_name:
          velocity: SteadyTaylorVortex

    - inflow_boundary_condition: bc_bottom
      target_name: surface_4
      inflow_user_data:
        user_function_name:
          velocity: SteadyTaylorVortex

    solution_options:
      name: myOptions
      turbulence_model: laminar
  
      options:

        - hybrid_factor:
            velocity: 0.0

        - limiter:
            pressure: no
            velocity: no

        - element_source_terms:
            momentum: [SteadyTaylorVortex, momentum_time_derivative]

        - consistent_mass_matrix_png:
            pressure: yes
            velocity: no

    solution_norm:
      output_frequency: 75
      file_name: steadyTaylorVortex_P4_R1.dat
      spacing: 12
      percision: 6
      target_name: block_1
      dof_user_function_pair:
       - [velocity, SteadyTaylorVortexVelocity]
       - [dpdx, SteadyTaylorVortexGradPressure]

    output:
      output_data_base_name: steadyTaylorVortex_P4_R1.e
      output_frequency: 75
      output_node_set: no 
      output_variables:
       - dpdx
       - velocity
       - pressure
       - velocity_exact
       - dpdx_exact

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 75
      time_step: 12.5e-4
      time_stepping_type: fixed 
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
