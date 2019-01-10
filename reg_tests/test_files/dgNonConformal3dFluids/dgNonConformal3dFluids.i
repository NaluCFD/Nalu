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
    muelu_xml_file_name: ../../xml/milestone.xml

realms:

  - name: realm_1
    mesh: ../../mesh/twoBlockMesh_cgs.g
    use_edges: no       

    equation_systems:
      name: theEqSys
      max_iterations: 2 
   
      solver_system_specification:
        pressure: solve_cont
        velocity: solve_scalar
     
      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:
      - constant: ic_1
        target_name: [block_1, block_2]
        value:
          pressure: 0
          velocity: [100.0,0.0]

    material_properties:
      target_name: [block_1, block_2]
      specifications:
        - name: density
          type: constant
          value: 1.0e-3

        - name: viscosity
          type: constant
          value: 1.8e-4

    boundary_conditions:

    - inflow_boundary_condition: bc_front
      target_name: surface_1
      inflow_user_data:
        velocity: [100.0,0.0,0.0]

    - open_boundary_condition: bc_back
      target_name: surface_2
      open_user_data:
        pressure: 0.0

    - wall_boundary_condition: bc_top
      target_name: surface_3
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_bot
      target_name: surface_4
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_right
      target_name: surface_5
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_left
      target_name: surface_6
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    - non_conformal_boundary_condition: bc_in_out
      target_name: [surface_7, surface_77]
      non_conformal_user_data:
        expand_box_percentage: 5.0 

    - non_conformal_boundary_condition: bc_out_in
      target_name: [surface_77, surface_7]
      non_conformal_user_data:
        expand_box_percentage: 5.0 

    solution_options:
      name: myOptions

      options:
        - hybrid_factor:
            velocity: 1.0
          
        - alpha_upw:
            velocity: 1.0

        - non_conformal:
            gauss_labatto_quadrature: no
            upwind_advection: yes
 
    output:
      output_data_base_name: dgElem.e
      output_frequency: 5
      output_node_set: no 
      output_variables:
       - velocity
       - pressure

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 100 
      time_step: 0.0010
      time_stepping_type: adaptive
      time_step_count: 0

      realms:
        - realm_1
