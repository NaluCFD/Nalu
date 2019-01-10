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
    tolerance: 1e-3
    max_iterations: 50
    kspace: 50
    output_level: 0
    muelu_xml_file_name: ../../xml/milestone.xml
    summarize_muelu_timer: no
    recompute_preconditioner: no

realms:

  - name: realm_1
    mesh: ../../mesh/hybrid_sphere.g
    use_edges: no 

    time_step_control:
     target_courant: 2.0
     time_step_change_factor: 1.2
   
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
            convergence_tolerance: 1e-5

        - MixtureFraction:
            name: myMixFrac
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:

      - constant: ic_1
        target_name: [block_1, block_2]
        value:
          pressure: 0.0
          velocity: [0.0,0.0,0.0]
          mixture_fraction: 0.0

    material_properties:
      target_name: [block_1, block_2]
      specifications:
        - name: density
          type: constant
          value: 1.18e-3

        - name: viscosity
          type: constant
          value: 1.8e-4

    boundary_conditions:

    - inflow_boundary_condition: bc_1
      target_name: surface_1
      inflow_user_data:
        velocity: [50.0,0.0,0.0]
        mixture_fraction: 0.0

    - open_boundary_condition: bc_2
      target_name: surface_2
      open_user_data:
        pressure: 0.0
        velocity: [0.0,0.0,0.0]
        mixture_fraction: 0.0

    - symmetry_boundary_condition: bc_3
      target_name: surface_3
      symmetry_user_data:

    - wall_boundary_condition: bc_5
      target_name: surface_5
      wall_user_data:
        velocity: [0.0,0.0]
        mixture_fraction: 1.0

    - overset_boundary_condition: bc_left
      overset_user_data:
        percent_overlap: 10.0 
        background_block: block_1
        overset_block: block_2
        overset_surface: surface_4
        background_cut_block: block_3
        background_cut_surface: surface_101

    solution_options:
      name: myOptions

      options:

        - hybrid_factor:
            velocity: 0.1
            mixture_fraction: 1.0

        - limiter:
            pressure: no
            velocity: no
            mixture_fraction: yes

        - projected_nodal_gradient:
            pressure: element
            velocity: element
            mixture_fraction: element
 
    output:
      output_data_base_name: hyrbid_sphere.e
      output_frequency: 10
      output_node_set: no 
      output_variables:
       - dual_nodal_volume
       - velocity
       - pressure
       - mixture_fraction
       - intersected_element

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 10
      time_step: 0.0001
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
