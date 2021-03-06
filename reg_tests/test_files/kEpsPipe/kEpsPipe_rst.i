Simulation:
  name: kEpsPipe_rst

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs 
    tolerance: 1e-5
    max_iterations: 75 
    kspace: 75 
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 75
    kspace: 75
    muelu_xml_file_name: ../../xml/milestone.xml
    output_level: 0
    recompute_preconditioner: no

realms:

  - name: realm_1
    mesh: kEpsPipe.rst
    use_edges: yes
   
    time_step_control:
     target_courant: 10.0
     time_step_change_factor: 1.2

    equation_systems:
      name: theEqSys
      max_iterations: 1

      solver_system_specification:
        velocity: solve_scalar
        turbulent_ke: solve_scalar
        turbulent_dissipation: solve_scalar
        mixture_fraction: solve_scalar
        pressure: solve_cont

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - KEpsilon:
            name: myKEps
            max_iterations: 1
            convergence_tolerance: 1e-5
            output_clipping_diagnostic: yes

        - MixtureFraction:
            name: myZ
            max_iterations: 1
            convergence_tolerance: 1.e-5

    initial_conditions:
      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0
          velocity: [0,0,10000.0]
          turbulent_ke: 1.5e4
          turbulent_dissipation: 2.15e7
          mixture_fraction: 0.0
  
    material_properties:
      target_name: block_1
      specifications:
        - name: density
          type: constant
          value: 1.1814e-3
        - name: viscosity
          type: constant
          value: 1.79e-4

    boundary_conditions:

    - inflow_boundary_condition: bc_inflow
      target_name: surface_1
      inflow_user_data:
        velocity: [0.0,0.0,10000.0]
        turbulent_ke: 1.5e6
        turbulent_dissipation: 2.15e9
        mixture_fraction: 1.0

    - open_boundary_condition: bc_open
      target_name: surface_2
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        turbulent_ke: 1.0e-8
        turbulent_dissipation: 1.0e-8
        mixture_fraction: 0.0

    - wall_boundary_condition: bc_bottom
      target_name: surface_3
      wall_user_data:
        velocity: [0,0,0]
        mixture_fraction: 0.0
        use_wall_function: yes

    solution_options:
      name: myOptions
      turbulence_model: k_epsilon

      use_consolidated_solver_algorithm: no
      use_consolidated_face_elem_bc_algorithm: no

      options:

        - hybrid_factor:
            velocity: 1.0 
            turbulent_ke: 1.0
            turbulent_dissipation: 1.0
            mixture_fraction: 1.0

        - limiter:
            velocity: yes 
            turbulent_ke: yes
            turbulent_dissipation: yes
            mixture_fraction: yes

        - turbulent_schmidt:
            turbulent_ke: 1.0
            turbulent_dissipation: 1.3
            mixture_fraction: 1.0

        - laminar_schmidt:
            turbulent_ke: 1.0
            turbulent_dissipation: 1.0
            mixture_fraction: 1.0

        - upw_factor:
            velocity: 1.0
            turbulent_ke: 1.0
            turbulent_dissipation: 1.0
            mixture_fraction: 1.0

        - projected_nodal_gradient:
            velocity: element 
            pressure: element
            turbulent_ke: element 
            turbulent_dissipation: element 
            mixture_fraction: element 

        - turbulence_model_constants:
            tkeProdLimitRatio: 1000.0
            cMu: 0.09001
            cEpsOne: 1.44001
            cEpsTwo: 1.92001

    post_processing:
    
    - type: surface
      physics: surface_force_and_moment_wall_function
      output_file_name: kEpsPipe_rst.dat
      frequency: 1
      parameters: [0,0,0]
      target_name: [surface_3]

    output:
      output_data_base_name: kEpsPipe_rst.e
      output_frequency: 25
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - turbulent_ke
       - turbulent_dissipation
       - mixture_fraction
       - scalar_variance
       - scalar_dissipation
       - tau_wall
       - yplus

    restart:
      restart_data_base_name: kEpsPipe_rst.rst
      restart_start: 25
      restart_frequency: 25
      restart_time: 10000.0
    
Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 100
      time_step: 1.0e-5
      time_stepping_type: adaptive 
      time_step_count: 0
      second_order_accuracy: no

      realms: 
        - realm_1
