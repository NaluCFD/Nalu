Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1
    error_estimator: errest_1

linear_solvers:

  - name: solve_adv_diff
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-4
    max_iterations: 200
    kspace: 200
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres
    preconditioner: muelu
    tolerance: 1e-4
    max_iterations: 200
    kspace: 200
    output_level: 0
    recompute_preconditioner: false
    muelu_xml_file_name: ../../xml/milestone.xml

realms:

  - name: realm_1
    mesh: ../../mesh/hex8_4_small.g
    use_edges: no
    polynomial_order: 5
    automatic_decomposition_type: rib

    equation_systems:
      name: theEqSys
      max_iterations: 2
   
      solver_system_specification:
        velocity: solve_adv_diff
        mixture_fraction: solve_adv_diff
        pressure: solve_cont
        dpdx: solve_adv_diff

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-2

        - MixtureFraction:
            name: myZ
            max_iterations: 1
            convergence_tolerance: 1.e-2
            clipping_delta: 0.01

    initial_conditions:
      - user_function: ic_1
        target_name: block_1
        user_function_name:
         velocity: VariableDensity
         pressure: VariableDensity
         mixture_fraction: VariableDensity

    material_properties:

      target_name: block_1

      specifications:
 
        - name: density
          type: mixture_fraction
          primary_value: 0.1
          secondary_value: 1.0

        - name: viscosity
          type: constant
          value: 0.001

    boundary_conditions:

    - inflow_boundary_condition: bc_left
      target_name: surface_1
      inflow_user_data:
        user_function_name:
          velocity: VariableDensity
          mixture_fraction: VariableDensity

    - inflow_boundary_condition: bc_right
      target_name: surface_2
      inflow_user_data:
        user_function_name:
          velocity: VariableDensity
          mixture_fraction: VariableDensity

    - inflow_boundary_condition: bc_top
      target_name: surface_3
      inflow_user_data:
        user_function_name:
          velocity: VariableDensity
          mixture_fraction: VariableDensity

    - inflow_boundary_condition: bc_bottom
      target_name: surface_4
      inflow_user_data:
        user_function_name:
         velocity: VariableDensity
         mixture_fraction: VariableDensity

    - inflow_boundary_condition: bc_front
      target_name: surface_5
      inflow_user_data:
        user_function_name:
         velocity: VariableDensity
         mixture_fraction: VariableDensity

    - inflow_boundary_condition: bc_back
      target_name: surface_6
      inflow_user_data:
        user_function_name:
          velocity: VariableDensity
          mixture_fraction: VariableDensity

    solution_options:
      name: myOptions
      turbulence_model: laminar
      divU_stress_scaling: 1.0
      interp_rhou_together_for_mdot: yes
      use_consolidated_solver_algorithm: no

      options:

        - hybrid_factor:
            velocity: 0.0
            mixture_fraction: 0.0

        - limiter:
            pressure: no
            velocity: no
            mixture_fraction: no

        - element_source_terms:
            momentum: [buoyancy, momentum_time_derivative, VariableDensity]
            continuity: [VariableDensity]
            mixture_fraction: [mixture_fraction_time_derivative, VariableDensity]

        - consistent_mass_matrix_png:
           pressure: yes
           velocity: no
           mixture_fraction: no

        - laminar_schmidt:
            mixture_fraction: 0.8

        - user_constants:
            gravity: [-5.0, 6.0, 7.0]
            reference_density: 1.0

    solution_norm:
      output_frequency: 11
      file_name: variableDensP5.dat
      spacing: 12
      percision: 6
      target_name: block_1
      dof_user_function_pair:
       - [velocity, VariableDensityVelocity]
       - [mixture_fraction, VariableDensityMixtureFraction]

    output:
      output_data_base_name: variableDensP5.e
      output_frequency: 11
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - mixture_fraction
       - density
       - mixture_fraction_exact
       - velocity_exact

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 11
      time_step: 1.5e-3
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
