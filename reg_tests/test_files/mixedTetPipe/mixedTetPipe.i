Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1

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
    tolerance: 1e-6
    max_iterations: 50
    kspace: 50
    output_level: 0
    muelu_xml_file_name: ../../xml/matches_ml_default.xml

realms:

  - name: realm_1
    mesh: ../../mesh/pipeTet.g
    use_edges: yes 

    time_step_control:
     target_courant: 2.0
     time_step_change_factor: 1.2
   
    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        turbulent_ke: solve_scalar

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            element_continuity_eqs: yes 
            convergence_tolerance: 1e-5

        - TurbKineticEnergy:
            name: myTke
            max_iterations: 1
            convergence_tolerance: 1.e-2

    initial_conditions:
      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0
          velocity: [0,0,0]
          turbulent_ke: 1.0e-6

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
        velocity: [0.0,0.0,100.0]
        turbulent_ke: 3.5

    - open_boundary_condition: bc_open
      target_name: surface_2
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        turbulent_ke: 0.0

    - wall_boundary_condition: bc_bottom
      target_name: surface_3
      wall_user_data:
        velocity: [0,0,0]
        use_wall_function: yes

    solution_options:
      name: myOptions
      turbulence_model: ksgs

      options:

        - turbulence_model_constants:
            kappa: 0.41
            cEps: 0.844
            cMuEps: 0.0855

        - hybrid_factor:
            velocity: 1.0 

        - alpha_upw:
            velocity: 1.0 
      
        - noc_correction:
            pressure: yes 

        - projected_nodal_gradient:
            pressure: element
            velocity: edge
            turbulent_ke: element
  
    output:
      output_data_base_name: mixedTetPipe.e
      output_frequency: 2 
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - turbulent_ke
       - turbulent_viscosity

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      time_step: 0.0005
      termination_time: 0.01
      time_stepping_type: adaptive
      time_step_count: 0

      realms: 
        - realm_1
