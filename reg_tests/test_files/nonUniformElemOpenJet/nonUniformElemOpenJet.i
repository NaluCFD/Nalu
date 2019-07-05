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
    muelu_xml_file_name: ../../xml/matches_ml_default.xml

realms:

  - name: realm_1
    mesh: ../../mesh/2cm_ped_35K_mks.g
    use_edges: no 

    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        mass_fraction: solve_scalar

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - MassFraction:
            name: myYk
            max_iterations: 1
            convergence_tolerance: 1e-5
            number_of_species: 2

    material_properties:

      target_name: block_1

      constant_specification:
       universal_gas_constant: 8314.4621
       reference_pressure: 101325.0
       reference_temperature: 298.15

      reference_quantities:
        - species_name: aO2
          mw: 32.0

        - species_name: bN2
          mw: 28.0

      specifications:
 
        - name: density
          type: ideal_gas

        - name: viscosity
          type: polynomial

          coefficient_declaration:
          - species_name: aO2
            coefficients: [1.7894e-5, 273.11, 110.56]

          - species_name: bN2
            coefficients: [1.7894e-5, 273.11, 110.56]

    initial_conditions:
      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0
          velocity: [0,0]  
          mass_fraction: [0.22,0.78]

    boundary_conditions:

    - wall_boundary_condition: bc_bottom
      target_name: surface_1
      wall_user_data:
        velocity: [0,0,0]

    - inflow_boundary_condition: bc_inflow
      target_name: surface_2
      inflow_user_data:
        velocity: [0,0,10.0]
        mass_fraction: [0.85, 0.15]

    - wall_boundary_condition: bc_pipe
      target_name: surface_3
      wall_user_data:
        velocity: [0,0,0]

    - open_boundary_condition: bc_side
      target_name: surface_4
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        mass_fraction: [0.22, 0.78]

    - open_boundary_condition: bc_top
      target_name: surface_5
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        mass_fraction: [0.22, 0.78]

    solution_options:
      name: myOptions
      turbulence_model: wale
      shift_cvfem_mdot: no
      reduced_sens_cvfem_poisson: yes 

      error_indicator:
        type: pstab
        frequency: 1
    
      options:
        - hybrid_factor:
            velocity: 1.0
            mass_fraction: 1.0

        - laminar_schmidt:
           mass_fraction: 0.85

        - turbulent_schmidt:
           mass_fraction: 1.0

        - source_terms:
            continuity: density_time_derivative

        - limiter:
            pressure: no
            velocity: no
            mass_fraction: yes

        - shifted_gradient_operator:
            velocity: no
            pressure: no


    turbulence_averaging:
      time_filter_interval: 100.0
      averaging_type: moving_exponential

      specifications:
        - name: one
          target_name: block_1
          compute_mean_error_indicator: yes

          moving_averaged_variables:
            - pressure
            - velocity
            - mass_fraction
            - turbulent_viscosity

    output:
      output_data_base_name: nonUniformElemOpenJet.e
      output_frequency: 5
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - turbulent_viscosity
       - mass_fraction
       - error_indicator
       - mean_error_indicator
       - pressure_ma
       - velocity_ma
       - mass_fraction_ma
       - turbulent_viscosity_ma

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 0.01
      time_step: 0.001
      time_stepping_type: adaptive 
      time_step_count: 0

      realms:
        - realm_1
