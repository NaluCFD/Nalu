Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres 
    preconditioner: sgs 
    tolerance: 1e-3
    max_iterations: 75 
    kspace: 75 
    output_level: 0

realms:

  - name: realm_1
    mesh: ../../mesh/oversetMeshAligned.g
    use_edges: no
    automatic_decomposition_type: rcb
    check_for_missing_bcs: yes

    equation_systems:
      name: theEqSys
      max_iterations: 2 
  
      solver_system_specification:
        temperature: solve_scalar
   
      systems:
        - HeatConduction:
            name: myHC
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:

      - user_function: ic_2
        target_name: [block_1, block_2]
        user_function_name:
         temperature: steady_2d_thermal

    material_properties:
      target_name: [block_1, block_2]
      specifications:
        - name: density
          type: constant
          value: 1.0
        - name: thermal_conductivity
          type: constant
          value: 1.0
        - name: specific_heat
          type: constant
          value: 1.0

    boundary_conditions:

    - wall_boundary_condition: bc_1
      target_name: surface_1
      wall_user_data:
        temperature: 298.0

    - wall_boundary_condition: bc_2
      target_name: surface_2
      wall_user_data:
        temperature: 298.0

    - wall_boundary_condition: bc_3
      target_name: surface_3
      wall_user_data:
        temperature: 298.0

    - wall_boundary_condition: bc_4
      target_name: surface_4
      wall_user_data:
        temperature: 298.0

    - wall_boundary_condition: bc_5
      target_name: surface_5
      wall_user_data:
        temperature: 500.0

    - overset_boundary_condition: bc_left
      overset_user_data:
        percent_overlap: 10.0 
        background_block: block_1
        overset_block: block_2
        overset_surface: surface_6
        background_cut_block: block_3
        background_cut_surface: surface_101

    solution_options:
      name: myOptions

      options:

        - projected_nodal_gradient:
            temperature: element 

    output:
      output_data_base_name: overset.e
      output_frequency: 2
      output_node_set: no 
      output_variables:
       - dual_nodal_volume
       - temperature
       - dtdx
       - intersected_element

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 20
      time_step: 0.5
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: no

      realms:
        - realm_1
