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
    muelu_xml_file_name: ../../xml/milestone.xml
    summarize_muelu_timer: no

realms:

  - name: realm_1
    mesh: ../../mesh/theRectangle.g
    use_edges: no
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:

      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0.0
          velocity: [2.0,0.0,0.0]

    material_properties:
      target_name: block_1
      specifications:
        - name: density
          type: constant
          value: 1.18e-3

        - name: viscosity
          type: constant
          value: 1.8e-4

    boundary_conditions:

    - inflow_boundary_condition: bc_1
      target_name: surface_4
      inflow_user_data:
        velocity: [2.0,0.0,0.0]

    - open_boundary_condition: bc_2
      target_name: surface_6
      open_user_data:
        pressure: 0.0
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_3
      target_name: surface_1
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_4
      target_name: surface_2
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_5
      target_name: surface_3
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_6
      target_name: surface_5
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    solution_options:
      name: myOptions
      use_consolidated_solver_algorithm: yes

      options:

        - hybrid_factor:
            velocity: 1.0

        - limiter:
            pressure: no
            velocity: no

        - projected_nodal_gradient:
            pressure: element
            velocity: element

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, upw_advection_diffusion, lumped_actuator]
            continuity: [advection]

    actuator:
      type: ActLinePointDrag
      search_method: boost_rtree
      search_target_part: block_1

      specifications:

        - turbine_name: machine_one
          radius: 0.5
          omega: 1.57
          gaussian_decay_radius: 1.0
          gaussian_decay_target: 0.01
          tip_coordinates: [2.5, 2.0, 0.0]
          tail_coordinates: [2.5, -2.0, 0.0]
          number_of_points: 11

        - turbine_name: machine_two
          radius: 0.5
          omega: -1.57
          gaussian_decay_radius: 1.0
          gaussian_decay_target: 0.01
          tip_coordinates: [-2.5, 0.0, 2.0]
          tail_coordinates: [-2.5, 0.0, -2.0]
          number_of_points: 11

    output:
      output_data_base_name: actuatorLine.e
      output_frequency: 5
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - actuator_line_source

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 20
      time_step: 0.01
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: no

      realms:
        - realm_1
