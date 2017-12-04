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

transfers:

  - name: xfer_thermal_rte_inner_outer
    type: geometric
    realm_pair: [realmThermal, realmRte]
    from_target_name: [surface_2, surface_3]
    to_target_name: [surface_1, surface_2]
    transfer_variables:
      - [temperature, temperature_bc]

  - name: xfer_rte_thermal_inner_outer
    type: geometric
    realm_pair: [realmRte, realmThermal]
    from_target_name: [surface_1, surface_2]
    to_target_name: [surface_2, surface_3]
    transfer_variables:
      - [irradiation, irradiation]
     
realms:

  - name: realmThermal
    mesh: ../../mesh/twoshells_coarse_mks.g
    use_edges: no
   
    boundary_conditions:

    - wall_boundary_condition: bc_inner
      target_name: surface_1
      wall_user_data:
        temperature: 300.0

    - wall_boundary_condition: bc_inner_int
      target_name: surface_2
      wall_user_data:
        emissivity: 0.8
        irradiation: 0.0
        interface: yes

    - wall_boundary_condition: bc_outer_int
      target_name: surface_3
      wall_user_data:
        emissivity: 0.8
        irradiation: 0.0
        interface: yes

    - wall_boundary_condition: bc_outer
      target_name: surface_4
      wall_user_data:
        temperature: 600.0

    solution_options:
      name: myOptionsHC
      options:
        - projected_nodal_gradient:
            temperature: element
        - user_constants:
            stefan_boltzmann: 5.6704e-8
                
    initial_conditions:

      - constant: ic_1
        target_name: block_1
        value:
          temperature: 300

    material_properties:
      target_name: block_1
      specifications:
        - name: density
          type: constant
          value: 10.0
        - name: specific_heat
          type: constant
          value: 10.0
        - name: thermal_conductivity
          type: constant
          value: 52.0

    equation_systems:
      name: theEqSys
      max_iterations: 1

      solver_system_specification:
        temperature: solve_scalar 

      systems:
        - HeatConduction:
            name: myHC
            max_iterations: 1 
            convergence_tolerance: 1.e-5

    output:
      output_data_base_name: thermalNew.e
      output_frequency: 5
      output_node_set: no 
      output_variables:
       - temperature
       - dtdx
       - density
       - thermal_conductivity
       - specific_heat
       - irradiation

  - name: realmRte
    mesh: ../../mesh/ether_conform_coarse_mks.g
    use_edges: no 
   
    boundary_conditions:

    - wall_boundary_condition: bc_1
      target_name: surface_1
      wall_user_data:
        temperature: 300.0
        transmissivity: 0.0
        emissivity: 0.80 
        interface: yes

    - wall_boundary_condition: bc_2
      target_name: surface_2
      wall_user_data:
        temperature: 300.0
        transmissivity: 0.0
        emissivity: 0.80 
        interface: yes

    initial_conditions:

      - constant: ic_1
        target_name: block_1
        value:
          temperature: 0.0

    material_properties:
      target_name: block_1
      specifications:
        - name: absorption_coefficient
          type: constant
          value: 0.0 
        - name: scattering_coefficient
          type: constant
          value: 0.0

    equation_systems:
      name: theEqSys
      max_iterations: 1
    
      solver_system_specification:
        intensity: solve_scalar

      systems:
        - RadiativeTransport:
            name: myRTE
            max_iterations: 1 
            convergence_tolerance: 1.e-16
            quadrature_order: 4 
            activate_scattering: yes

    solution_options:
      name: myOptions
      options:
        - user_constants:
            stefan_boltzmann: 5.6704e-8

        - element_source_terms:
            intensity: [advection_sucv, absorption_black_body, isotropic_scattering]

    output:
      output_data_base_name: pmrNew.e
      output_frequency: 5
      output_node_set: no
      output_variables:
       - dual_nodal_volume
       - temperature
       - absorption_coefficient
       - intensity_bc
       - scalar_flux
       - radiative_heat_flux
       - div_radiative_heat_flux
       - irradiation

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 20.0 
      time_step: 1.0
      time_stepping_type: fixed
      time_step_count: 0
      nonlinear_iterations: 1 

      realms:
        - realmThermal
        - realmRte

      transfers:
        - xfer_thermal_rte_inner
        - xfer_thermal_rte_outer
        - xfer_rte_thermal_inner
        - xfer_rte_thermal_outer
