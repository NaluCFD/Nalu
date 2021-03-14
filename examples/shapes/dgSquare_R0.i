Simulation:
  name: NaluSim

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
    muelu_xml_file_name: milestone.xml

realms:

  - name: realm_1
    mesh: mesh/square_R0.g
    use_edges: no       
    automatic_decomposition_type: rcb
   
    time_step_control:
     target_courant: 2.0
     time_step_change_factor: 1.15
   
    equation_systems:
      name: theEqSys
      max_iterations: 2
   
      solver_system_specification:
        pressure: solve_cont
        velocity: solve_scalar
        mixture_fraction: solve_scalar

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - MixtureFraction:
            name: myZ
            max_iterations: 1
            convergence_tolerance: 1.e-2
            output_clipping_diagnostic: yes

    initial_conditions:
      - constant: ic_1
        target_name: [block_1, block_2]
        value:
          pressure: 0
          velocity: [1.0,0.0]
          mixture_fraction: 0.0

    material_properties:
      target_name: [block_1, block_2]
      specifications:
        - name: density
          type: constant
          value: 1.2

        - name: viscosity
          type: constant
          value: 0.016

    boundary_conditions:

    - inflow_boundary_condition: bc_left
      target_name: surface_1
      inflow_user_data:
        velocity: [1.0,0.0]
        mixture_fraction: 0.0

    - open_boundary_condition: bc_right
      target_name: surface_2
      open_user_data:
        pressure: 0.0
        mixture_fraction: 0.0

    - symmetry_boundary_condition: bc_top
      target_name: surface_3
      symmetry_user_data:

    - symmetry_boundary_condition: bc_bot
      target_name: surface_4
      symmetry_user_data:

    - non_conformal_boundary_condition: bc_in_out
      current_target_name: [surface_5]
      opposing_target_name: [surface_55]
      non_conformal_user_data:
        expand_box_percentage: 5.0 
        search_tolerance: 0.01 
        activate_dynamic_search_algorithm: yes

    - non_conformal_boundary_condition: bc_out_in
      current_target_name: [surface_55]
      opposing_target_name: [surface_5]
      non_conformal_user_data:
        expand_box_percentage: 5.0 
        search_tolerance: 0.01 
        activate_dynamic_search_algorithm: yes

    - wall_boundary_condition: bc_bluff
      target_name: surface_6
      wall_user_data:
        user_function_name:
         velocity: wind_energy
        user_function_string_parameters:
         velocity: [mmFront_ss6]
        mixture_fraction: 1.0

    solution_options:
      name: myOptions

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      mesh_motion:

        - name: mmBackground
          target_name: [block_1]
          omega: 0.0

        - name: mmFront_ss6
          target_name: [block_2]
          omega: 0.785398
          unit_vector: [0.0,0.0,1.0]
          compute_centroid: yes

      options:

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, advection_diffusion, NSO_2ND]
            continuity: advection
            mixture_fraction: [lumped_mixture_fraction_time_derivative, advection_diffusion, NSO_2ND]

        - non_conformal:
            gauss_labatto_quadrature: no
            upwind_advection: yes
            current_normal: yes
            include_png_penalty: no

        - consistent_mass_matrix_png:
            pressure: no
            velocity: no
            mixture_fraction: no

    output:
      output_data_base_name: output/dgSquare_R0.e
      output_frequency: 10
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - mixture_fraction
       - mesh_displacement

    restart:
      restart_data_base_name: restart/dgSquare_R0.rst
      restart_frequency: 1000
      restart_start: 1000

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 1000
      time_step: 0.02
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
