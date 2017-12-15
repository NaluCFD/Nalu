
Wind Farm and Atmospheric Boundary Layer (ABL) Modeling
=======================================================

Wind energy analysis is the primary application area for the Nalu development
team and, therefore, Nalu includes all necessary models that would necessary to
simulate a wind farm on complex terrain under realistic atmospheric boundary
layer conditions. This section describes the theoretical basis of Nalu from a
wind energy perspective.

In order to evaluate the energy output and the structural loading on wind
turbines, the code must model: 1. the incoming turbulent wind field across the
entire wind farm, 2. the evolution of turbine wakes in turbulent inflow
conditions and their interaction with the downstream turbines, and 3. the
effects of complex terrain on the evolution of the incoming wind field as well
as the turbine wakes. Furthermore, the effects of Coriolis forces and the
buoyancy forces resulting from surface heating must be accounted for to get an
accurate estimate of the flow conditions throughout the wind farm.

Wind farm simulations can be broken down into two categories: *precursor*
simulations, and simulations with prescribed turbulent inflow velocity and
temperature profiles. *Precursor* simulations are used to trigger turbulence
generation and generate inflow velocity profiles that are used as inlet
conditions subsequent wind farm simulations. An alternative approach to
*precursor* simulations is to use inflow conditions from a mesoscale simulation
model, e.g., WRF.

Governing Equations
-------------------

We begin with a review of the momentum and enthalpy conservation equations
within the context of wind farm modeling. Equation :eq:`ablmom` shows the
Favre-filtered momentum conservation equation (Eq. :eq:`favmom`) reproduced here with
all the terms required to model a wind farm.

.. math::
   :label: ablmom

   \underbrace{\frac{\partial}{\partial t} \left(\bar{\rho}\, \widetilde{u}_i\right)}_\mathbf{I} +
   \underbrace{\frac{\partial}{\partial x_j} \left( \bar{\rho}\, \widetilde{u}_i \widetilde{u}_j \right)}_\mathbf{II} =
     - \underbrace{\frac{\partial p'}{\partial x_j} \delta_{ij}}_\mathbf{III}
     - \underbrace{\frac{\partial \tau_{ij}}{\partial x_j}}_\mathbf{IV}
     - \underbrace{2\bar{\rho}\,\epsilon_{ijk}\,\Omega_ju_k}_\mathbf{V}
     + \underbrace{\left(\bar{\rho} - \rho_\circ \right) g_i}_\mathbf{VI}
     + \underbrace{S^{u}_{i}}_\mathbf{VII} + \underbrace{f^{T}_i}_\mathbf{VIII}

Term :math:`\mathbf{I}` represents the time rate of change of momentum (inertia);

Term :math:`\mathbf{II}` represents advection;

Term :math:`\mathbf{III}` represents the pressure gradient forces (deviation from
hydrostatic and horizontal mean gradient);

Term :math:`\mathbf{IV}` represents stresses (both viscous and sub-filter scale
(SFS)/Reynolds stresses);

Term :math:`\mathbf{V}` describes the influence Coriolis forces due to earth's rotation -- see  Sec. :numref:`earth_coriolis_force`;

Term :math:`\mathbf{VI}` describes the effects of buoyancy using the Boussinesq approximation -- see :numref:`boussinesq_buoyancy_model`;

Term :math:`\mathbf{VII}` represents the source term used to drive the flow to a
horizontal mean velocity at desired height(s) -- see :numref:`abl_forcing_term`; and

Term :math:`\mathbf{VIII}` is an optional term representing body forces when
modeling turbine with actuator disk or line representations -- see :numref:`theory_actuator_wind_turbine_models`.

In wind energy applications, the energy conservation equation is often written
in terms of the Favre-filtered potential temperature equation, as shown below

.. math::
   :label: abl_pottemp

   \frac{\partial}{\partial t} \left(\bar{\rho}\, \widetilde{\theta}\right) +
   \frac{\partial}{\partial t} \left(\bar{\rho}\, \widetilde{u}_j \widetilde{\theta} \right) = - \frac{\partial}{\partial x_j} \hat{q}_j

where, :math:`q_j` represents the temperature transport due to molecular and SFS
turbulence effects. Due to the high Reynolds number associated with ABL flows,
the molecular effects are neglected everywhere except near the terrain.
Potential temperature is related to absolute temperature by the following
equation

.. math::

   \theta = T \left ( \frac{\bar{p}}{p_\circ} \right)^{-\left(\frac{R}{c_p}\right)}

Under the assumption of ideal gas conditions and constant :math:`c_p`, the gradients in
potential temperature are proportional to the gradients in absolute temperature,
i.e.,

.. math::

   \left[ \frac{\partial T}{\partial t}, \frac{\partial T}{\partial x}, \frac{\partial T}{\partial y} \right] =
   \left( \frac{\bar{p}}{p_\circ} \right)^\left(\frac{R}{c_p}\right) \left[ \frac{\partial \theta}{\partial t}, \frac{\partial \theta}{\partial x}, \frac{\partial \theta}{\partial y} \right]

Furthermore, ignoring the pressure and viscous work terms in Eq. :eq:`fav-enth`
and assuming constant density (incompressible flow), it can be shown that
solving the enthalpy equation is equivalent to solving the potential temperature
equation. Care must be taken to scale the SFS flux terms appropriately in the
equations, and appropriate initial conditions and boundary conditions for
potential temperature must be provided. The resulting solution can then be
interpreted as the variation of potential temperature field in the computational
domain.

Initial & Boundary Conditions
-----------------------------

This section briefly describes the boundary conditions available in Nalu for
modeling wind farm problems. The terrain and top boundary conditions are
described first as they are common to precusor and wind farm simulations.

Initial conditions
~~~~~~~~~~~~~~~~~~

Nalu has the ability to initialize the internal flow fields to uniform
conditions for all pressure, velocity, temperature, and TKE (:math:`k`) in the
:inpfile:`input file <initial_conditions.constant>`. Nalu also provides a *user
function* to add perturbations to the velocity field to trigger turbulence
generation during precursor simulations. To specify more complex flow field
conditions, a temperature profile with a capping inversion for example, users
are referred to pre-processing utilities available in `NaluWindUtils
<http://naluwindutils.readthedocs.io/en/latest/>`_ library.

Terrain (Wall) boundary condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users are referred to :numref:`abl_surface_conditions` for the treatment of the
terrain BC using roughness models. For enthalpy, users can provide a surface heat
flux for modeling stratified flows.

Top boundary condition
~~~~~~~~~~~~~~~~~~~~~~

For momentum, a :ref:`symmetry BC <theory_symmetry_bc>` is used when modeling
wind farm problems. For enthalpy equation, a normal temperature gradient can be
specified to drive the flow to a desired temperaure profile, e.g., capping
inversion temperature profile.

Inlet conditions
~~~~~~~~~~~~~~~~

Time histories of inflow velocity and temperaure profiles can be provided as
inputs (via I/O transfer) to drive the wind farm simulation with the desired
flow conditions. See :numref:`verification_abl_prescribed_inflow` for more
details on this capability. Driving a wind farm simulation using velocity and
temperature fields from a mesoscale (WRF) simulation would require an additional
pre-processing steps with the `wrftonalu
<http://naluwindutils.readthedocs.io/en/latest/user/wrftonalu.html>`_ utility.

Outlet conditions
~~~~~~~~~~~~~~~~~

See the description of :ref:`open BC <theory_open_bc>` for detailed description
of the outlet BC implementation. For wind energy problems, it is necessary to
activate the global mass correction as a single value of pressure across the
boundary layer is not apprpriate in the presence of buoyancy effects. It might
also be necessary to fix the reference pressure at an interior node in order to
ensure that the Pressure Poisson solver is well conditioned.
