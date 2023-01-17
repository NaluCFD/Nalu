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

realms:

  - name: fluidRealm
    mesh: mesh/1d_quad4_adv_diff.exo
    use_edges: no
    automatic_decomposition_type: rcb

    equation_systems:
      name: theEqSys
      max_iterations: 5

      solver_system_specification:
        volume_of_fluid: solve_scalar

      systems:

        - VolumeOfFluid:
            name: myV
            max_iterations: 1
            convergence_tolerance: 1e-2
            stand_alone_equation_system: yes
            clipping_delta: 1.0

    initial_conditions:

      - constant: icConst
        target_name: block_1
        value:
          velocity: [1.0,0.0]

      - user_function: icUser
        target_name: block_1
        user_function_name:
         volume_of_fluid: gaussian
        user_function_parameters:
         volume_of_fluid: [0, 1.0, 0.0, 0.1]

    material_properties:
      target_name: [block_1]

      specifications:

        - name: viscosity
          type: constant
          value: 0.01

    boundary_conditions:

    - periodic_boundary_condition: bc_left_right
      target_name: [surface_1, surface_2]
      periodic_user_data:
        search_tolerance: 1.e-5

    - symmetry_boundary_condition: bc_wrap
      target_name: surface_3
      symmetry_user_data:

    solution_options:
      name: myOptions

      options:

        - element_source_terms:
           volume_of_fluid: [lumped_mass, scs_upw_advection_np, diffusion]

        - limiter:
            volume_of_fluid: yes

        - upw_factor:
            volume_of_fluid: 1.0
                
    output:
      output_data_base_name: output/1d_quad4_adv_diff.e
      output_frequency: 1
      output_node_set: no 
      output_variables:
       - volume_of_fluid
       - velocity
       - dvofdx
 
Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 200
      time_step: 0.01
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - fluidRealm
