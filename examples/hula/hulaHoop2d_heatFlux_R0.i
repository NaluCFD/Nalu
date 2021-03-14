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
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

realms:

  - name: realm_1
    mesh: hulaHoop2d_mks_R0.g
    use_edges: no
    automatic_decomposition_type: rcb

    time_step_control:
     target_courant: 1.0
     time_step_change_factor: 1.05

    equation_systems:
      name: theEqSys
      max_iterations: 2

      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont
        enthalpy: solve_scalar

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - Enthalpy:
            name: myEnth
            max_iterations: 1
            convergence_tolerance: 1e-5
            minimum_temperature: 299.999
     
    material_properties:

      target_name: [block_1, block_2]

      specifications:
 
# properties at T = 300

        - name: density
          type: constant
          value: 994.59

        - name: viscosity
          type: constant
          value: 8.62e-4

        - name: specific_heat
          type: constant
          value: 4.142e3

        - name: thermal_conductivity
          type: constant
          value: 6.103e-1

    initial_conditions:

      - constant: ic_top
        target_name: block_1
        value:
          pressure: 0
          velocity: [0,0]  
          temperature: 300.0

      - constant: ic_bottom
        target_name: block_2
        value:
          pressure: 0
          velocity: [0,0]  
          temperature: 300.0
  
    boundary_conditions:

    - wall_boundary_condition: bc_top
      target_name: surface_1
      wall_user_data:
        velocity: [0,0,0]
        heat_flux: -200.0

    - wall_boundary_condition: bc_bottom
      target_name: surface_2
      wall_user_data:
        velocity: [0,0,0]
        heat_flux: +200.0

    solution_options:
      name: myOptions
      turbulence_model: laminar 

      use_consolidated_solver_algorithm: yes

      options:

        - turbulent_prandtl:
            enthalpy: 0.90

        - user_constants:
            gravity: [0.0,-9.81]
            reference_density: 994.59
            reference_temperature: 300.0
            thermal_expansion_coefficient: 322.55e-6

        - limiter:
            pressure: no
            velocity: no
            enthalpy: yes 

        - peclet_function_form:
            velocity: tanh
            enthalpy: tanh

        - peclet_function_tanh_transition:
            velocity: 5000.0
            enthalpy: 2.0

        - peclet_function_tanh_width:
            velocity: 200.0
            enthalpy: 1.0

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, upw_advection_diffusion, buoyancy_boussinesq]
            continuity: [advection]
            enthalpy: [lumped_enthalpy_time_derivative, upw_advection_diffusion]

    turbulence_averaging:
      time_filter_interval: 100000.0

      specifications:

        - name: one
          target_name: [block_1, block_2]
          reynolds_averaged_variables:
            - velocity
            - temperature

    output:
      output_data_base_name: output_R0/hulaHoop2d_heatFlux_R0.e
      output_frequency: 25
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - enthalpy
       - temperature
       - specific_heat
       - thermal_conductivity
       - viscosity
       - velocity_ra_one
       - temperature_ra_one

    restart:
      restart_data_base_name: restart_R0/hulaHoop2d_heatFlux_R0.rst
      restart_frequency: 2000
      restart_start: 2000

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 5000.0
      time_step: 0.1
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
