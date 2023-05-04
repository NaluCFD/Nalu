#  Notes on target part names
#   [Unspecified-2-QUAD, Unspecified-3-TRIANGLE] = [block_202, block_303]
#   Inflow   = surface_1
#   Outflow  = surface_2
#   Sides    = surface_3
#   Cylinder = surface_4

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
    muelu_xml_file_name: ../../../xml/milestone_aspect_ratio_smooth.xml

realms:

  - name: realm_1
    mesh: mesh/2d_tri3_quad4_street.exo
    use_edges: no 
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
            mixture_fraction: 2.0

        - peclet_function_tanh_width:
            velocity: 4.0
            mixture_fraction: 4.0

        - projected_nodal_gradient:
            velocity: element
            pressure: element
            mixture_fraction: element

        - noc_correction: # active for edge-based only
            velocity: yes
            pressure: yes
            mixture_fraction: yes

    post_processing:
    
    - type: surface
      physics: surface_force_and_moment
      output_file_name: 2d_tri3_quad4_street.dat
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
      output_data_base_name: output/2d_tri3_quad4_street.e
      output_frequency: 20 
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - mixture_fraction
       - uf_mixture_fraction
       - velocity_ra_one
       - pressure_ra_one
       - mixture_fraction_ra_one
       - vorticity
       - q_criterion

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 100
      time_step: 1.0e-3
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
