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

      constant_specification:
       universal_gas_constant: 8314.4621e3
       reference_pressure: 101325.0
       reference_temperature: 298.15

      reference_quantities:
        - species_name: N2
          mw: 28.0
          mass_fraction: 0.767

        - species_name: O2
          mw: 32.0
          mass_fraction: 0.233

      specifications:

        - name: density
          type: ideal_gas

        - name: viscosity
          type: polynomial

          coefficient_declaration:
          - species_name: O2 
            coefficients: [1.7894e-4, 273.11, 110.56]

          - species_name: N2 
            coefficients: [1.7894e-4, 273.11, 110.56]

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
       - density
       - viscosity
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
