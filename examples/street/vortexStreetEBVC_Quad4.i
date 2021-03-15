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
    recompute_preconditioner: false
    muelu_xml_file_name: milestone_aspect_ratio_smooth.xml

realms:

  - name: realm_1
    mesh: circular_cylinder_coarse.exo
    use_edges: yes 
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 2 

      solver_system_specification:
        velocity: solve_scalar
        mixture_fraction: solve_scalar
        pressure: solve_cont

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
        target_name: [Unspecified-2-QUAD, Unspecified-3-TRIANGLE]
        value:
          pressure: 101325.0
          velocity: [0.0,0.0]
          mixture_fraction: 0.0

    material_properties:
      target_name: [Unspecified-2-QUAD, Unspecified-3-TRIANGLE]

      specifications:

        - name: density
          type: constant 
          value: 1.2

        - name: viscosity
          type: constant 
          value: 0.008

    boundary_conditions:

    - inflow_boundary_condition: bc_left
      target_name: Inflow
      inflow_user_data:
        velocity: [1.0,0.0]
        mixture_fraction: 0.0
        
    - open_boundary_condition: bc_right
      target_name: Outflow
      open_user_data:
        velocity: [0.0,0.0]
        pressure: 101325.0
        mixture_fraction: 0.0

    - symmetry_boundary_condition: bc_top_bottom
      target_name: Sides
      symmetry_user_data:

    - wall_boundary_condition: bc_cylinder
      target_name: Cylinder
      wall_user_data:
        velocity: [0.0,0.0]
        mixture_fraction: 1.0

    solution_options:
      name: myOptions
      turbulence_model: laminar 

      options:

        - hybrid_factor:
            velocity: 1.0
            mixture_fraction: 0.0

        - limiter:
            velocity: yes
            mixture_fraction: no

        - laminar_schmidt:
            mixture_fraction: 0.9

        - peclet_function_form:
            velocity: tanh
            mixture_fraction: tanh

        - peclet_function_tanh_transition:
            velocity: 2.0
            mixture_fraction: 1000.0

        - peclet_function_tanh_width:
            velocity: 4.0
            mixture_fraction: 100.0

        - noc_correction:
            pressure: yes

    post_processing:
    
    - type: surface
      physics: surface_force_and_moment
      output_file_name: vortexStreetEBVC_Quad4.dat
      frequency: 10
      parameters: [0.0,0.0]
      target_name: Cylinder

    turbulence_averaging:
      time_filter_interval: 1000.0
      specifications:
        - name: one
          target_name: [Unspecified-2-QUAD, Unspecified-3-TRIANGLE]
          reynolds_averaged_variables:
            - velocity
            - pressure
            - mixture_fraction

          compute_vorticity: yes 
          compute_q_criterion: yes

    output:
      output_data_base_name: output/vortexStreetEBVC_Quad4.e
      output_frequency: 20 
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - mixture_fraction
       - velocity_ra_one
       - pressure_ra_one
       - mixture_fraction_ra_one
       - vorticity
       - q_criterion

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 2000
      time_step: 1.0e-3
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
