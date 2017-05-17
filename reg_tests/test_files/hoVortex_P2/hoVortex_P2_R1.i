# default parameters for the various aux function
defaults:
  - &centroidX -2.5  # initial location of vortex in x direction
  - &centroidY 0.0   # initial location of vortex in y direction
  - &rVortex  0.25   # initial vortex radius
  - &beta     15.0   # initial vortex strength
  - &uInf     10.0   # x-velocity speed of the vortex
  - &dens     1.0e-3 # density
  - &visc     1.0e-4 # viscosity

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
    max_iterations: 100
    kspace: 100
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
    muelu_xml_file_name: ../../xml/muelu_STV_HO.xml

realms:

  - name: realm_1
    mesh: ../../mesh/tquad4_80.g
    use_edges: no
    automatic_decomposition_type: rib
    polynomial_order: 2

    equation_systems:
      name: theEqSys
      max_iterations: 4
   
      solver_system_specification:
        pressure: solve_cont
        velocity: solve_scalar
        dpdx: solve_scalar
        duidx: solve_scalar

      systems:
        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-8

    initial_conditions:

      - user_function: ic_1
        target_name: block_1
        user_function_name:
         velocity: wind_energy_taylor_vortex
         pressure: wind_energy_taylor_vortex
        user_function_parameters:
          velocity: 
            - *centroidX
            - *centroidY
            - *rVortex
            - *beta
            - *uInf
            - *dens
            - *visc

          pressure:
            - *centroidX
            - *centroidY
            - *rVortex
            - *beta
            - *uInf
            - *dens

    material_properties:
      target_name: block_1

      specifications:
        - name: density
          type: constant
          value: *dens

        - name: viscosity
          type: constant
          value: *visc

    boundary_conditions:

    - inflow_boundary_condition: bc_left
      target_name: surface_1
      inflow_user_data:
        velocity: [*uInf, 0.0, 0.0]

    - open_boundary_condition: bc_right
      target_name: surface_2
      open_user_data:
        pressure: 0

    - symmetry_boundary_condition: bc_top
      target_name: surface_3
      symmetry_user_data:

    - symmetry_boundary_condition: bc_bottom
      target_name: surface_4
      symmetry_user_data:

    solution_options:
      name: myOptions
      turbulence_model: laminar

      options:
        - hybrid_factor:
            velocity: 0.0

        - limiter:
            pressure: no
            velocity: no

        - element_source_terms:
            momentum: momentum_time_derivative

        - consistent_mass_matrix_png:
            pressure: yes
            velocity: no

    output:
      output_data_base_name: hoVortex_P2_R1.e
      output_frequency: 10
      output_node_set: no 
      output_variables:
       - dual_nodal_volume
       - velocity
       - velocity_exact
       - pressure
       - dpdx
       - dpdx_exact

    solution_norm:
      output_frequency: 1
      file_name: hoVortex_P2_R1.dat
      spacing: 12
      percision: 6
      target_name: block_1
      dof_user_function_pair:
       - [velocity, wind_energy_taylor_vortex]
       - [dpdx, wind_energy_taylor_vortex_dpdx]

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_time: 0.01
      time_step: 2.5e-4
      time_stepping_type: fixed 
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
