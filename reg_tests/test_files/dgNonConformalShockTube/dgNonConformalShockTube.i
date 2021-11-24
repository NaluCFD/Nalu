Simulation:
  name: NaluSim

linear_solvers:

  - name: solve_REMOVE_REQUIRED
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-4
    max_iterations: 50
    kspace: 50
    output_level: 0

realms:

  - name: fluidRealm
    mesh: ../../mesh/shock_tube_0p001_quad4_nc.exo
    use_edges: yes 
    automatic_decomposition_type: rcb

    time_step_control:
     target_courant: 0.75
     time_step_change_factor: 1.2

    equation_systems:
      name: theEqSys
      max_iterations: 1

      solver_system_specification:
        gas_dynamics: solve_REMOVE_REQUIRED

      systems:
        - GasDynamics:
            name: myGD
            max_iterations: 1
            convergence_tolerance: 1e-2

    initial_conditions:

      - constant: ic_1
        target_name: block_1
        value:
          pressure: 101325.0
          velocity: [0.0,0.0]
          temperature: 298.15

      - constant: ic_2
        target_name: block_2
        value:
          pressure: 10132.50
          velocity: [0.0,0.0]
          temperature: 298.15

    material_properties:
      target_name: [block_1, block_2]

      constant_specification:
       universal_gas_constant: 8314.4621
       reference_temperature: 0.0

      reference_quantities:
        - species_name: FakeAir
          mw: 28.0
          mass_fraction: 1.0

      specifications:

        - name: density
          type: ideal_gas

        - name: viscosity
          type: constant
          value: 0.0

        - name: thermal_conductivity
          type: constant
          value: 0.0

        - name: gamma
          type: constant
          value: 1.4

    boundary_conditions:

    - symmetry_boundary_condition: bc_left
      target_name: surface_1
      symmetry_user_data:

    - symmetry_boundary_condition: bc_right
      target_name: surface_2
      symmetry_user_data:

    - symmetry_boundary_condition: bc_top
      target_name: surface_3
      symmetry_user_data:

    - symmetry_boundary_condition: bc_bot
      target_name: surface_4
      symmetry_user_data:

    - non_conformal_boundary_condition: bc_left
      target_name: [surface_5, surface_55]
      non_conformal_user_data:
        expand_box_percentage: 10.0 

    - non_conformal_boundary_condition: bc_right
      target_name: [surface_55, surface_5]
      non_conformal_user_data:
        expand_box_percentage: 10.0 

    solution_options:
      name: myOptions
      turbulence_model: laminar

      use_accoustically_compressible_algorithm: yes

      options:

        - non_conformal:
            gauss_labatto_quadrature: yes
            current_normal: yes

    data_probes:

      output_frequency: 100

      search_tolerance: 1.0e-3
      search_expansion_factor: 2.0

      specifications:

        - name: probe_volume 
          from_target_part: [block_1, block_2]

          line_of_site_specifications:

            - name: dgNonConformalShockTube
              number_of_points: 1001
              tip_coordinates: [+0.5, 0.0]
              tail_coordinates: [-0.5, 0.0]

          output_variables:
            - field_name: density
              field_size: 1
            - field_name: pressure
              field_size: 1
            - field_name: temperature
              field_size: 1
            - field_name: mach_number
              field_size: 1
            - field_name: velocity
              field_size: 2

    output:
      output_data_base_name: dgNonConformalShockTube.e
      output_frequency: 20
      output_node_set: no
      output_variables:
       - density
       - momentum
       - total_energy
       - velocity
       - total_enthalpy
       - enthalpy
       - pressure
       - temperature
       - mach_number
       - speed_of_sound
       - viscosity
       - specific_heat
       - thermal_conductivity
       - gamma

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 10.0e-5
      time_step: 5.0e-7
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: no

      realms:
        - fluidRealm
