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
    recompute_preconditioner: no
    muelu_xml_file_name: ../../../xml/matches_ml_default.xml

realms:

  - name: fluidRealm
    mesh: mesh/2d_quad9_couette.exo
    use_edges: no
    automatic_decomposition_type: rcb

    time_step_control:
     target_courant: 2.0
     time_step_change_factor: 1.2
   
    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        dpdx: solve_scalar

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:
      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0
          velocity: [0,0]

    material_properties:
      target_name: block_1
      specifications:

        - name: density
          type: constant
          value: 1000.0

        - name: viscosity
          type: constant
          value: 8.9e-4

    boundary_conditions:

    - periodic_boundary_condition: bc_left_right
      target_name: [surface_1, surface_2]
      periodic_user_data:
        search_tolerance: 0.0001 

    - wall_boundary_condition: bc_top
      target_name: surface_3
      wall_user_data:
        velocity: [1.0e-3,0]

    - wall_boundary_condition: bc_bot
      target_name: surface_4
      wall_user_data:
        velocity: [0,0]

    solution_options:
      name: myOptions
      turbulence_model: laminar

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      options:
      
        - element_source_terms:
            momentum: [momentum_time_derivative, advection_diffusion, NSO_2ND_ALT, body_force]
            continuity: advection

        - consistent_mass_matrix_png:
            pressure: yes

        - element_source_term_parameters:
            momentum: [0.0,0.0]

    output:
      output_data_base_name: output/2d_quad9_couette.e
      output_frequency: 10
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - dpdx
       - density
       - viscosity

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      time_step: 0.1
      termination_step_count: 20
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms: 
        - fluidRealm
