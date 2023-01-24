Simulation:
  name: NaluSim

linear_solvers:

  - name: solve_adv_diff
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
    recompute_preconditioner: false
    muelu_xml_file_name: ../../../xml/milestone_aspect_ratio_smooth.xml

realms:

  - name: realm_1
    mesh: mesh/2d_quad9_helium.exo
    use_edges: no 
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 4

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
            output_clipping_diagnostic: yes

    initial_conditions:
      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0
          velocity: [0,0]
          mixture_fraction: 0.0

    material_properties:
      target_name: block_1

      specifications:

        - name: density
          type: mixture_fraction
          primary_value: 0.163
          secondary_value: 1.18

        - name: viscosity
          type: mixture_fraction
          primary_value: 1.967e-5
          secondary_value: 1.85e-5

    boundary_conditions:

    - inflow_boundary_condition: bc_inflow
      target_name: surface_1
      inflow_user_data:
        velocity: [0.0,0.340,0.0]
        mixture_fraction: 1.0

    - wall_boundary_condition: bc_bottom
      target_name: surface_2
      wall_user_data:
        velocity: [0,0,0]
        
    - open_boundary_condition: bc_side
      target_name: surface_3
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        mixture_fraction: 0.0
        use_total_pressure: yes

    - open_boundary_condition: bc_top
      target_name: surface_4
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        mixture_fraction: 0.0
        use_total_pressure: yes

    solution_options:
      name: myOptions
      turbulence_model: laminar

      divU_stress_scaling: 1.0
   
      reduced_sens_cvfem_poisson: no

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      options:

        - laminar_schmidt:
            mixture_fraction: 0.9

        - turbulent_schmidt:
            mixture_fraction: 1.0

        - element_source_terms:
            momentum: [momentum_time_derivative, advection_diffusion, buoyancy, NSO_2ND_ALT]
            continuity: [density_time_derivative, advection]
            mixture_fraction: [mixture_fraction_time_derivative, advection_diffusion, NSO_2ND]

        - user_constants:
            gravity: [0.0,-9.81,0.0]
            reference_density: 1.18

        - consistent_mass_matrix_png:
            pressure: yes

    turbulence_averaging:
      time_filter_interval: 10.0
      specifications:
        - name: one
          target_name: block_1
          reynolds_averaged_variables:
            - velocity
            - mixture_fraction
          favre_averaged_variables:
            - velocity
            - mixture_fraction

    explicit_filtering:

      search_target_part: block_1
      filter_size: [0.05, 0.05]
      debug_output: no

      explicit_specifications:

        - field_name: velocity
          explicit_field_name: explicit_velocity
          field_size: 2
          field_type: node_rank
          field_state_size: 3

        - field_name: dudx
          explicit_field_name: explicit_dudx
          field_size: 4
          field_type: node_rank
          field_state_size: 1

        - field_name: pressure
          explicit_field_name: explicit_pressure
          field_size: 1
          field_type: node_rank
          field_state_size: 1

        - field_name: dpdx
          explicit_field_name: explicit_dpdx
          field_size: 2
          field_type: node_rank
          field_state_size: 1

        - field_name: mixture_fraction
          explicit_field_name: explicit_mixture_fraction
          field_size: 1
          field_type: node_rank
          field_state_size: 3

        - field_name: density
          explicit_field_name: explicit_density
          field_size: 1
          field_type: node_rank
          field_state_size: 3

        - field_name: viscosity
          explicit_field_name: explicit_viscosity
          field_size: 1
          field_type: node_rank
          field_state_size: 1

        - field_name: density_ra_one
          explicit_field_name: explicit_density_ra_one
          field_size: 1
          field_type: node_rank
          field_state_size: 1

      residual_specifications:

        - field_name: residual
          field_size: 2
          field_type: node_rank

        - field_name: explicit_residual
          field_size: 2
          field_type: node_rank

    output:
      output_data_base_name: output/2d_quad9_helium.e
      output_frequency: 50
      output_node_set: no
      output_variables:
       - velocity
       - dudx
       - pressure
       - mixture_fraction
       - density
       - density_ra_one
       - velocity_ra_one
       - velocity_fa_one
       - mixture_fraction_ra_one
       - mixture_fraction_fa_one
       - explicit_velocity
       - explicit_dudx
       - explicit_pressure
       - explicit_dpdx
       - explicit_mixture_fraction
       - explicit_density
       - explicit_viscosity
       - explicit_density_ra_one
       - residual
       - explicit_residual
       - explicit_filter
       - dual_nodal_volume


Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 0.1
      time_step: 0.001
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
