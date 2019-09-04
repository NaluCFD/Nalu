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
    max_iterations: 100
    kspace: 100
    output_level: 0
    recompute_preconditioner: no
    muelu_xml_file_name: ../../xml/milestone.xml

realms:

  - name: realm_1
    mesh: ../../mesh/1x2x10_tet10_0p2.g
    activate_fem: yes
    use_edges: no 
    automatic_decomposition_type: rib     

    equation_systems:
      name: theEqSys
      max_iterations: 1

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

      - user_function: ic_1
        target_name: block_1
        user_function_name:
         velocity: OneTwoTenVelocity

      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0.0

    material_properties:
      target_name: block_1

      specifications:

        - name: density
          type: constant
          value: 1.0e-3

        - name: viscosity
          type: constant
          value: 1.0e-4

    boundary_conditions:

    - open_boundary_condition: bc_left
      target_name: surface_1
      open_user_data:
        velocity: [0.0,0.0,0.0]
        pressure: 0.016

    - open_boundary_condition: bc_right
      target_name: surface_2
      open_user_data:
        velocity: [0.0,0.0,0.0]
        pressure: 0.0

    - wall_boundary_condition: bc_walls
      target_name: surface_3
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    solution_options:
      name: myOptions
      turbulence_model: laminar

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      options:

        - hybrid_factor:
            velocity: 0.0

        - element_source_terms:
            momentum: [momentum_time_derivative, advection, diffusion]
            continuity: advection
            dpdx: interior_png

        - consistent_mass_matrix_png:
            pressure: yes

    solution_norm:
      output_frequency: 100
      file_name: 1x2x10Tet10.dat
      spacing: 12
      percision: 6
      target_name: block_1
      dof_user_function_pair:
       - [velocity, OneTwoTenVelocity]

    output:
      output_data_base_name: 1x2x10Tet10.e
      output_frequency: 5
      output_node_set: no 
      output_variables:
       - pressure
       - velocity
       - dpdx
       - velocity_exact

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 20
      time_step: 0.05
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
