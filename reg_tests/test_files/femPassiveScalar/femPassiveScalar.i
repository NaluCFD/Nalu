Simulation:
  name: NaluFemWorkshipSimulation

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

realms:

  - name: scalarRealm
    mesh: ../../mesh/hex8_20_10_1.g
    activate_fem: yes
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        mixture_fraction: solve_scalar
        dzdx: solve_scalar

      systems:

        - MixtureFractionFEM:
            name: myZ
            max_iterations: 1
            convergence_tolerance: 1.e-2
            compute_png: yes

    initial_conditions:

      - constant: ic_1
        target_name: block_1
        value:
          velocity: [1.0,0.0,0.0]

      - user_function: ic_2
        target_name: block_1
        user_function_name:
         mixture_fraction: workshop_mms

    material_properties:
      target_name: block_1

      specifications:

        - name: density
          type: mixture_fraction
          primary_value: 1.0
          secondary_value: 1.0

        - name: viscosity
          type: mixture_fraction
          primary_value: 1.0e-3
          secondary_value: 1.0e-3

    boundary_conditions:

    - periodic_boundary_condition: bc_front_back
      target_name: [surface_1, surface_2]
      periodic_user_data:
        search_tolerance: 0.0001 

    - periodic_boundary_condition: bc_left_right
      target_name: [surface_3, surface_4]
      periodic_user_data:
        search_tolerance: 0.0001 

    - periodic_boundary_condition: bc_top_bottom
      target_name: [surface_5, surface_6]
      periodic_user_data:
        search_tolerance: 0.0001 

    solution_options:
      name: myOptions
      turbulence_model: laminar

      options:

        - element_source_terms:
            mixture_fraction: [mixture_fraction_time_derivative, advection, diffusion]
            dzdx: interior_png

        - laminar_schmidt:
            mixture_fraction: 1.0

        - consistent_mass_matrix_png:
            mixture_fraction: yes

    output:
      output_data_base_name: femPassiveScalar.e
      output_frequency: 10
      output_variables:
       - velocity
       - mixture_fraction
       - density
       - dzdx

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 50
      time_step: 1.0e-2
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - scalarRealm
