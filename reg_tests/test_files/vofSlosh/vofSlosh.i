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
    max_iterations: 75
    kspace: 75
    output_level: 0
    recompute_preconditioner: true
    muelu_xml_file_name: ../../xml/milestone_aspect_ratio_smooth_s.xml

realms:

  - name: fluidRealm
    mesh: ../../mesh/slosh_obs_0p005.exo
    use_edges: no
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        velocity: solve_scalar
        volume_of_fluid: solve_scalar
        pressure: solve_cont

      systems:

        - VolumeOfFluid:
            name: myV
            max_iterations: 1
            convergence_tolerance: 1.e-2
            activate_smoothing: yes
            fourier_number: 0.25
            smoothing_iterations: 5
            compression_constant: 0.1

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-2

    initial_conditions:

      - constant: ic_top
        target_name: block_1
        value:
          velocity: [0,0]
          volume_of_fluid: 1.0
          pressure: 0.0

      - constant: ic_bot
        target_name: block_2
        value:
          velocity: [0,0]
          volume_of_fluid: 0.0
          pressure: 0.0

    material_properties:
      target_name: [block_1, block_2]

      specifications:

        - name: density
          type: volume_of_fluid
          phase_one: 1000.0
          phase_two: 1.0

        - name: viscosity
          type: volume_of_fluid
          phase_one: 8.9e-4
          phase_two: 1.8e-5

        - name: surface_tension
          type: constant
          value: 0.07

    boundary_conditions:

    - wall_boundary_condition: bc_left
      target_name: surface_1
      wall_user_data:
        velocity: [0.0,0.0]

    - wall_boundary_condition: bc_right
      target_name: surface_2
      wall_user_data:
        velocity: [0.0,0.0]

    - open_boundary_condition: bc_top
      target_name: surface_3
      open_user_data:
        velocity: [0.0,0.0]
        pressure: 0.0
        volume_of_fluid: 0.0

    - wall_boundary_condition: bc_bot
      target_name: surface_4
      wall_user_data:
        velocity: [0.0,0.0]

    solution_options:
      name: myOptions

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      activate_balanced_force_algorithm: yes
      activate_buoyancy_pressure_stabilization: no

      options:

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, advection_diffusion, buoyancy, NSO_2ND_ALT]
            continuity: [vof_advection]
            volume_of_fluid: [lumped_mass, scs_advection, sucv_nso, sharpen]

        - user_constants:
            gravity: [0.0,-9.81]
            reference_density: 0.0

    output:
      output_data_base_name: vofSlosh.e
      output_frequency: 20
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - volume_of_fluid
       - volume_of_fluid_smoothed
       - interface_normal
       - interface_curvature
       - surface_tension
       - dvofdx
       - density
       - viscosity
 
Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 200
      time_step: 0.001
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - fluidRealm
