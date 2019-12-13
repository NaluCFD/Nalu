Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1
    error_estimator: errest_1

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
    max_iterations: 100
    kspace: 100
    output_level: 0
    recompute_preconditioner: no
    muelu_xml_file_name: ../../xml/milestone.xml

realms:

  - name: realm_1
    mesh: ../../mesh/channel_Retau1000_yp31_dxp196_dzp131.g
    use_edges: no
    automatic_decomposition_type: rcb

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
            convergence_tolerance: 1e-5

        - TurbKineticEnergy:
            name: myTke
            max_iterations: 1
            convergence_tolerance: 1.0e-5

    initial_conditions:
      - user_function: icUser
        target_name: Unspecified-2-HEX
        user_function_name:
         velocity: ChannelFlowPerturbedPlug
        user_function_parameters:
         velocity: [1.0, 6.283185307, 2.0,2.094395102, 0, 1, 2]

      - constant: ic_const
        target_name: Unspecified-2-HEX
        value:
          turbulent_ke: 1.0e-8
          pressure: 101325.0
        
    material_properties:

      target_name: Unspecified-2-HEX

      specifications:
 
        - name: density
          type: constant
          value: 1.0
          
        - name: viscosity
          type: constant
          value: 5.0e-5

    boundary_conditions:

    - periodic_boundary_condition: bc_periodic_x
      target_name: [Inflow,Outflow]
      periodic_user_data:
        search_tolerance: 0.000001 

    - wall_boundary_condition: bc_wall_bot
      target_name: BottomWall
      wall_user_data:
        velocity: [0.0,0.0,0.0]
        use_wall_function_projected: yes
        projected_distance: 0.0625

    - wall_boundary_condition: bc_wall_top
      target_name: TopWall
      wall_user_data:
        velocity: [0.0,0.0,0.0]
        use_wall_function_projected: yes
        projected_distance: 0.0625

    - periodic_boundary_condition: bc_periodic_z
      target_name: [Left,Right]
      periodic_user_data:
        search_tolerance: 0.000001

    solution_options:
      name: myOptions
      turbulence_model: ksgs

      use_consolidated_solver_algorithm: yes
      use_consolidated_face_elem_bc_algorithm: yes

      options:

        - element_source_terms:
            momentum: [lumped_momentum_time_derivative, advection_diffusion]
            continuity: [advection]
            turbulent_ke: [lumped_turbulent_ke_time_derivative, upw_advection_diffusion, ksgs]

        - source_terms:
            momentum: body_force

        - source_term_parameters:
            momentum: [0.0025,0.0,0.0]

        - hybrid_factor:
            turbulent_ke: 1.0

        - limiter:
            turbulent_ke: yes

    post_processing:
    
    - type: surface
      physics: surface_force_and_moment_wall_function_projected
      output_file_name: re1kChannel.dat
      frequency: 1
      parameters: [0,0,0]
      target_name: [BottomWall, TopWall]

    turbulence_averaging:
      time_filter_interval: 10000000.0
      specifications:
        - name: one
          target_name: Unspecified-2-HEX
          reynolds_averaged_variables:
            - velocity
            - pressure
            - dkdx
            - turbulent_viscosity
            - resolved_turbulent_ke

          compute_tke: yes 
          compute_reynolds_stress: yes

        - name: two
          target_name: [BottomWall,TopWall]
          reynolds_averaged_variables:
            - yplus
            - tau_wall

    output:
      output_data_base_name: re1kChannel.e
      output_frequency: 25
      output_node_set: no 
      output_variables:
       - velocity
       - turbulent_viscosity
       - pressure
       - turbulent_ke
       - dpdx
       - dudx
       - dkdx_ra_one
       - tau_wall
       - yplus
       - velocity_ra_one
       - pressure_ra_one
       - turbulent_viscosity_ra_one
       - resolved_turbulent_ke_ra_one
       - reynolds_stress
       - yplus
       - tau_wall
       - yplus_ra_two
       - tau_wall_ra_two
       - assembled_area_force_moment_wfp

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 100
      time_step: 0.04
      time_stepping_type: fixed
      time_step_count: 0
      second_order_accuracy: yes 

      realms:
        - realm_1
