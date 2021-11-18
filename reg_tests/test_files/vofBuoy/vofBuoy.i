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
    muelu_xml_file_name: ../../xml/vofBuoy.xml

realms:

  - name: realm_1
    mesh: ../../mesh/buoy.g
    use_edges: no       
    automatic_decomposition_type: rcb

    time_step_control:
     target_courant: 2.0
     time_step_change_factor: 1.15
   
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
            convergence_tolerance: 1.e-3
            activate_smoothing: yes
            fourier_number: 0.25
            smoothing_iterations: 5
            compression_constant: 0.1

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1 
            convergence_tolerance: 1e-3

    initial_conditions:
      - constant: ic_1
        target_name: [oset-blk, bg-blk]
        value:
          pressure: 0
          velocity: [0.0,0.0,0.0]

      - user_function: ic_2
        target_name: [oset-blk, bg-blk]
        user_function_name:
          volume_of_fluid: FixedHeight
        user_function_parameters:
          volume_of_fluid: [0.0, 0.07]

    material_properties:
      target_name: [oset-blk, bg-blk]
      specifications:
        - name: density
          type: volume_of_fluid
          phase_one: 1.0
          phase_two: 1000.0

        - name: viscosity
          type: volume_of_fluid
          phase_one: 1.8e-5
          phase_two: 8.9e-4

        - name: surface_tension
          type: constant
          value: 0.07

    boundary_conditions:
    
    - wall_boundary_condition: bg_bc
      target_name: bgwall
      wall_user_data:
        velocity: [0,0,0]

    - open_boundary_condition: bg_bc_open
      target_name: bgopen
      open_user_data:
        velocity: [0,0,0]
        pressure: 0
        volume_of_fluid: 1.0

    - overset_boundary_condition: buoy_edge
      overset_user_data:
        percent_overlap: 15.0
        percent_overlap_inner: 35.0
        background_block: bg-blk
        overset_block: oset-blk
        overset_surface: oset-edge
        background_cut_block: block_10
        background_cut_surface: surface_20
        cutting_shape: obb
        detailed_output: no

    - wall_boundary_condition: bc_buoy
      target_name: buoy
      wall_user_data:
        user_function_name:
         velocity: mesh_motion
        user_function_string_parameters:
         velocity: [mmSphere_ss5]

    solution_options:
      name: myOptions


      mesh_motion:

        - name: mmBackground
          target_name: [bg-blk]
          omega: 0.0

        - name: mmSphere_ss5
          target_name: [oset-blk]
          include_six_dof: yes
          body_omega: [0.0,0.0,0.0]
          body_mass: 167.88
          body_velocity: [-1.0,0,0]
          body_angle: [0.0,0.0,1.3]
          body_cc_disp: [1.0,0.5,0]
          principal_moments_inertia: [0.0,0.0,50.0]
          forcing_surface: [buoy]
          applied_force: [0.0,-1646.90,0.0]
          compute_centroid: yes

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      activate_balanced_force_algorithm: yes

      options:

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, advection_diffusion, NSO_2ND_ALT]
            continuity: [vof_advection]
            volume_of_fluid: [lumped_mass, scs_advection, sucv_nso, sharpen]

        - user_constants:
            gravity: [0.0, -9.81]
            reference_density: 0.0

    output:
      output_data_base_name: buoy_drop.e
      output_frequency: 20
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - mesh_displacement
       - fringe_node
       - dpdx
       - dudx
       - intersected_element
       - mesh_velocity
       - density
       - viscosity
       - volume_of_fluid_smoothed
       - volume_of_fluid

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 20
      time_step: 0.002
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
