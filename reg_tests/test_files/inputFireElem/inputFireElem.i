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
    max_iterations: 100
    kspace: 100
    output_level: 0
    write_matrix_files: no

transfers:

# io to ....

  - name: xfer_io_pmr
    type: geometric
    realm_pair: [ioRealm, pmrRealm]
    from_target_name: [surface_1, surface_2, surface_3, surface_4]
    to_target_name: [surface_1, surface_2, surface_3, surface_4]
    objective: initialization
    transfer_variables:
      - [temperature, temperature_bc]

realms:

  - name: pmrRealm
    mesh: ../../mesh/13K_oneWay.e
    use_edges: no
   
    boundary_conditions:

    - wall_boundary_condition: bc_entrainment
      target_name: surface_4
      wall_user_data:
        temperature: 298.0
        transmissivity: 0.0
        emissivity: 1.0
        interface: yes

    - wall_boundary_condition: bc_open
      target_name: surface_3
      wall_user_data:
        temperature: 298.0
        transmissivity: 0.0
        emissivity: 1.0 
        interface: yes 

    - wall_boundary_condition: bc_wall
      target_name: surface_2
      wall_user_data:
        temperature: 298.0
        transmissivity: 0.0
        emissivity: 1.0 
        interface: yes

    - wall_boundary_condition: bc_pool
      target_name: surface_1
      wall_user_data:
        transmissivity: 0.0
        emissivity: 0.8
        temperature: 500.0
        interface: yes

    initial_conditions:

      - constant: ic_1
        target_name: block_1
        value:
          temperature: 298.0

    material_properties:
      target_name: block_1
      specifications:
        - name: absorption_coefficient
          type: constant
          value: 0.001 
        - name: scattering_coefficient
          type: constant
          value: 0.0

    equation_systems:
      name: theEqSys
      max_iterations: 20 
    
      solver_system_specification:
        intensity: solve_scalar

      systems:
        - RadiativeTransport:
            name: myRTE
            max_iterations: 1 
            convergence_tolerance: 1.e-8
            quadrature_order: 4 
            activate_scattering: no
            external_coupling: yes
            activate_upwind: no
            deactivate_sucv: no

    solution_options:
      name: myOptions

      options:

        - element_source_terms:
            intensity: [advection_sucv, absorption_black_body]

        - input_variables_from_file:
            temperature: tnd
            absorption_coefficient: and
            radiation_source: empnd

    output:
      output_data_base_name: inputFireElemSucv.e
      output_frequency: 1
      output_node_set: no
      output_variables:
       - dual_nodal_volume
       - temperature
       - absorption_coefficient
       - intensity_bc
       - intensity_0
       - intensity_5
       - intensity_10
       - intensity_15
       - intensity_20
       - intensity_25
       - intensity_30
       - intensity_35
       - scalar_flux
       - radiative_heat_flux
       - radiation_source
       - div_radiative_heat_flux
       - irradiation

  - name: ioRealm
    mesh: ../../mesh/13K_oneWay.e
    type: initialization

    field_registration:
      specifications:
        - field_name: temperature
          target_name: [surface_1, surface_2, surface_3, surface_4]
          field_size: 1
          field_type: node_rank

    solution_options:
      name: myOptions
      input_variables_from_file_restoration_time: 1000.0

      options:    

        - input_variables_from_file:
            temperature: tnd
 
    output:
      output_data_base_name: IO.e
      output_frequency: 2
      output_node_set: no
      output_variables:
       - temperature_bc

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 1 
      time_step: 0.5
      time_stepping_type: fixed
      time_step_count: 0
     
      realms: 
        - pmrRealm
        - ioRealm
