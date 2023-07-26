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
    muelu_xml_file_name: ../../../xml/milestone_aspect_ratio_smooth.xml

realms:

  - name: fluidRealm
    mesh: mesh/3d_hex8_dam_break.exo
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
            vof_density_phase_one: 1000.0
            vof_density_phase_two: 1.0

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-2

    initial_conditions:

      - user_function: icUser
        target_name: block_1
        user_function_name:
         volume_of_fluid: rectangle

      - constant: ic_top
        target_name: block_1
        value:
          velocity: [0,0,0]
          pressure: 0.0

    material_properties:
      target_name: [block_1]

      specifications:

        - name: density
          type: volume_of_fluid
          phase_one: 1000.0
          phase_two: 1.0

        - name: viscosity
          type: volume_of_fluid
          phase_one: 1.0e-3
          phase_two: 1.98e-5

        - name: surface_tension
          type: constant
          value: 0.07

    boundary_conditions:

    - wall_boundary_condition: bc_bottom
      target_name: surface_1
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_wrap
      target_name: surface_2
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - open_boundary_condition: bc_open
      target_name: surface_3
      open_user_data:
        velocity: [0.0,0.0,0.0]
        pressure: 0.0
        volume_of_fluid: 0.0

    - wall_boundary_condition: bc_front
      target_name: surface_4
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    solution_options:
      name: myOptions

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      activate_balanced_force_algorithm: yes
      activate_buoyancy_pressure_stabilization: yes

      local_vof_m: 1.0
      local_vof_n: 1.0
      local_vof_c: 6.0

      options:

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, advection_diffusion, NSO_2ND_ALT, sharpen, capillary]
            continuity: [vof_advection]
            volume_of_fluid: [lumped_mass, scs_advection, sucv_nso, sharpen]

        - user_constants:
            gravity: [0.0, 0.0, -9.81]
            reference_density: 0.0

    data_probes:

      output_frequency: 5

      search_tolerance: 1.0e-3
      search_expansion_factor: 2.0

      specifications:

        - name: probe_volume 
          from_target_part: block_1

          line_of_site_specifications:
            - name: H3_Pstab_R0 
              number_of_points: 101
              tip_coordinates: [0.582, 0.0, 0.8]
              tail_coordinates: [0.582, 0.0, 0.0]

            - name: H2_Pstab_R0 
              number_of_points: 101
              tip_coordinates: [1.732, 0.0, 0.8]
              tail_coordinates: [1.732, 0.0, 0.0]

            - name: H1_Pstab_R0 
              number_of_points: 101
              tip_coordinates: [2.228, 0.0, 0.8]
              tail_coordinates: [2.228, 0.0, 0.0]

          output_variables:
            - field_name: volume_of_fluid
              field_size: 1
            - field_name: volume_of_fluid_smoothed
              field_size: 1

    output:
      output_data_base_name: output/3d_hex8_dam_break.e
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
      termination_time: 0.10
      time_step: 0.00250
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - fluidRealm
