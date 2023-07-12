Simulation:
  name: NaluSim

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-5
    max_iterations: 200
    kspace: 200
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0
    muelu_xml_file_name: ../../xml/milestone_aspect_ratio_smooth_s.xml

transfers:

# Fluid to ....
          
# MM to ....

  - name: xfer_mm_fluid_surf
    type: geometric
    realm_pair: [mmRealm, fluidRealm]
    mesh_part_pair: [surface_5, surface_5]
    transfer_variables:
      - [mesh_velocity, velocity_bc]

  - name: xfer_mm_fluid_vol
    type: geometric
    realm_pair: [mmRealm, fluidRealm]
    from_target_name: [block_1, block_2]
    to_target_name: [block_1, block_2]
    transfer_variables:
      - [mesh_velocity, mesh_velocity]
      - [mesh_displacement, mesh_displacement]
      - [div_mesh_velocity, div_mesh_velocity]

realms:

  - name: mmRealm
    mesh: ../../mesh/waveGenerator_R0.exo
    use_edges: no
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 2
    
      solver_system_specification:
        mesh_displacement: solve_scalar

      systems:
        - MeshDisplacement:
            name: myMM
            max_iterations: 1 
            convergence_tolerance: 1.e-3
            deform_wrt_model_coordinates: no

    boundary_conditions:

    - wall_boundary_condition: bc_left
      target_name: surface_1
      wall_user_data:
        mesh_displacement: [0.0,0.0]

    - wall_boundary_condition: bc_right
      target_name: surface_2
      wall_user_data:
        mesh_displacement: [0.0,0.0]

    - wall_boundary_condition: bc_bottom
      target_name: surface_3
      wall_user_data:
        mesh_displacement: [0.0,0.0]

    - wall_boundary_condition: bc_top
      target_name: surface_4
      wall_user_data:
        mesh_displacement: [0.0,0.0]

    - wall_boundary_condition: bc_cyliner
      target_name: surface_5
      wall_user_data:
        user_function_name:
         mesh_displacement: sinusoidal
        user_function_parameters:
         mesh_displacement: [0.08]

    solution_options:
      name: myOptions

      options:
        - projected_nodal_gradient:
            mesh_velocity: element

        - element_source_terms:
            mesh_displacement: [mesh_disp_lumped, elastic_stress]
    
    initial_conditions:

      - constant: ic_1
        target_name: [block_1, block_2]
        value:
          mesh_displacement: [0.0,0.0]

    material_properties:
      target_name: [block_1, block_2]
      specifications:
        - name: lame_mu
          type: geometric
        - name: lame_lambda
          type: geometric
        - name: density
          type: constant
          value: 1.0e3

    output:
      output_data_base_name: output/waveGenerator_mm.e
      output_frequency: 25
      output_node_set: no
      output_variables:
       - mesh_displacement
       - mesh_displacement_bc
       - dxTmp
       - mesh_velocity
       - dvdx

  - name: fluidRealm
    mesh: ../../mesh/waveGenerator_R0.exo
    use_edges: no 
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        volume_of_fluid: solve_scalar

      systems:

        - VolumeOfFluid:
            name: myV
            max_iterations: 1
            convergence_tolerance: 1.e-2
            activate_smoothing: yes
            fourier_number: 0.25
            smoothing_iterations: 5
            compression_constant: 0.25

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:

      - constant: ic_top
        target_name: block_1
        value:
          velocity: [0,0]
          volume_of_fluid: 1.0
          pressure: 101325.0

      - constant: ic_bot
        target_name: block_2
        value:
          velocity: [0,0]
          volume_of_fluid: 0.0
          pressure: 101325.0

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

    - wall_boundary_condition: bc_1
      target_name: surface_1
      wall_user_data:
        velocity: [0.0,0.0]

    - wall_boundary_condition: bc_2
      target_name: surface_2
      wall_user_data:
        velocity: [0.0,0.0]

    - wall_boundary_condition: bc_3
      target_name: surface_3
      wall_user_data:
         velocity: [0.0,0.0]

    - wall_boundary_condition: bc_4
      target_name: surface_4
      wall_user_data:
         velocity: [0.0,0.0]

    - wall_boundary_condition: bc_5
      target_name: surface_5
      wall_user_data:
         velocity: [0.0,0.0]
         fsi_interface: yes

    solution_options:
      name: myOptions
      externally_provided_mesh_deformation: yes

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      activate_balanced_force_algorithm: yes
      activate_buoyancy_pressure_stabilization: yes

      options:

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, advection_diffusion, lumped_gcl, NSO_2ND_ALT]
            continuity: [vof_advection, lumped_gcl]
            volume_of_fluid: [lumped_mass, scs_advection, sucv_nso, sharpen, lumped_gcl]

        - projected_nodal_gradient:
            momentum: element
            continuity: element
            volume_of_fluid: element

        - user_constants:
            gravity: [9.81, 0.0]
            reference_density: 0.0
    
    output:
      output_data_base_name: output/waveGenerator_f.e
      output_frequency: 10
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
       - mesh_displacement

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 100
      time_step: 0.002
      time_stepping_type: fixed
      time_step_count: 0
      nonlinear_iterations: 1
      second_order_accuracy: yes

      realms:
        - mmRealm
        - fluidRealm

      transfers:
        - xfer_mm_fluid_surf
        - xfer_mm_fluid_vol
