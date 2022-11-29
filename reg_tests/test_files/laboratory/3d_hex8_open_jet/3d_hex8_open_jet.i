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
    max_iterations: 50
    kspace: 50
    output_level: 0
    recompute_preconditioner: no
    muelu_xml_file_name: ../../../xml/matches_ml_default.xml

realms:

  - name: fluidRealm
    mesh: mesh/3d_hex8_open_jet.exo
    use_edges: yes
    automatic_decomposition_type: rcb

    time_step_control:
     target_courant: 10.0
     time_step_change_factor: 1.2
   
    equation_systems:
      name: theEqSys
      max_iterations: 1

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        mixture_fraction: solve_scalar

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - MixtureFraction:
            name: myZ
            max_iterations: 1
            convergence_tolerance: 1.e-2

    initial_conditions:
      - constant: ic_1
        target_name: block_1
        value:
          pressure: 101325.0
          velocity: [0.0, 0.0, 0.0]
          mixture_fraction: 0.0

    material_properties:
      target_name: block_1
      specifications:

        - name: density
          type: constant
          value: 1.0e-3

        - name: viscosity
          type: constant
          value: 2.0e-4

    boundary_conditions:

    - wall_boundary_condition: bc_bottom
      target_name: surface_1
      wall_user_data:
        velocity: [0.0, 0.0, 0.0]

    - inflow_boundary_condition: bc_jet
      target_name: surface_2
      inflow_user_data:
        velocity: [0.0, 0.0, 20.0]
        mixture_fraction: 1.0

    - open_boundary_condition: bc_wrap
      target_name: surface_3
      open_user_data:
        velocity: [0.0, 0.0, 0.0]
        pressure: 101325.0
        mixture_fraction: 0.0

    - open_boundary_condition: bc_top
      target_name: surface_4
      open_user_data:
        velocity: [0.0, 0.0, 0.0]
        pressure: 101325.0
        mixture_fraction: 0.0

    solution_options:
      name: myOptions
      turbulence_model: laminar

      options:

        - peclet_function_form:
            velocity: tanh
            mixture_fraction: tanh

        - peclet_function_tanh_transition:
            velocity: 60.0
            mixture_fraction: 2.0

        - peclet_function_tanh_width:
            velocity: 20.0
            mixture_fraction: 4.02

        - laminar_schmidt:
            mixture_fraction: 0.9

        - limiter:
            velocity: no
            mixture_fraction: yes

        - projected_nodal_gradient:
            velocity: element
            pressure: element
            mixture_fraction: element

        - noc_correction:
            velocity: no
            pressure: no
            mixture_fraction: no
          
    output:
      output_data_base_name: output/3d_hex8_open_jet.e
      output_frequency: 10
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - dpdx
       - density
       - viscosity
       - mixture_fraction

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      time_step: 1.0e-3
      termination_step_count: 20
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: no

      realms: 
        - fluidRealm
