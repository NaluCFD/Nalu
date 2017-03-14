Low Mach Number Derivation
--------------------------

The low Mach number equations are a subset of the fully compressible
equations of motion (momentum, continuity and energy), admitting large
variations in gas density while remaining acoustically incompressible.
The low Mach number equations are preferred over the full compressible
equations for low speed flow problems as the accoustics are of little
consequence to the overall simulation accuracy. The technique avoids the
need to resolve fast-moving acoustic signals. Derivations of the low
Mach number equations can be found in found in Rehm and
Baum, :cite:`Rehm:1978`, or
Paolucci, :cite:`Paolucci:1982`.

The equations are derived from the compressible equations using a
perturbation expansion in terms of the lower limit of the Mach number
squared; hence the name. The asymptotic expansion leads to a splitting
of pressure into a spatially constant thermodynamic pressure and a
locally varying dynamic pressure. The dynamic pressure is decoupled from
the thermodynamic state and cannot propagate acoustic waves. The
thermodynamic pressure is used in the equation of state and to determine
thermophysical properties. The thermodynamic pressure can vary in time
and can be calculated using a global energy balance.

Asymptotic Expansion
++++++++++++++++++++

The asymptotic expansion for the low Mach number equations begins with
the full compressible equations in Cartesian coordinates. The equations
are the minimum set required to propagate acoustic waves. The equations
are written in divergence form using Einstein notation (summation over
repeated indices):

.. math::

   {{\partial \rho} \over {\partial t}} + {{\partial \rho u_j} 
                                  \over {\partial x_j}} & =  0 , \\
   {{\partial \rho u_i} \over {\partial t}} + {{\partial \rho u_j u_i}
       \over {\partial x_j}} + {{\partial P} \over {\partial x_i}} & = 
       {{\partial \tau_{ij}} \over {\partial x_j}} + \rho g_i , \\
   {{\partial \rho E} \over {\partial t}} + {{\partial \rho u_j H}
       \over {\partial x_j}}  & = 
       - {{\partial q_j} \over {\partial x_j}}
       + {{\partial u_i \tau_{ij}} \over {\partial x_j}} + \rho u_i g_i .


The primitive variables are the velocity components, :math:`u_i`, the
pressure, :math:`P`, and the temperature :math:`T`. The viscous shear
stress tensor is :math:`\tau_{ij}`, the heat conduction is :math:`q_i`,
the total enthalpy is :math:`H`, the total internal energy is :math:`E`,
the density is :math:`\rho`, and the gravity vector is :math:`g_i`. The
total internal energy and total enthalpy contain the kinetic energy
contributions. The equations are closed using the following models and
definitions:

.. math::

   P & = \rho {R \over W} T , \\
   E & = H - P/\rho , \\
   H & = h + {1 \over 2} u_k u_k , \\
   \tau_{ij} & = \mu \left( {{\partial u_i} \over {\partial x_j}}
                     +        {{\partial u_j} \over {\partial x_i}} \right)
               - {2 \over 3} \mu {{\partial u_k} \over {\partial x_k}} 
                  \delta_{ij} , \\
   q_i & = - k {{\partial T} \over {\partial x_i}} .

The mean molecular weight of the gas is :math:`W`, the molecular
viscosity is :math:`\mu`, and the thermal conductivity is :math:`k`. A
Newtonian fluid is assumed along with the Stokes hypothesis for the
stress tensor.

The equations are scaled so that the variables are all of order one. The
velocities, lengths, and times are nondimensionalized by a
characteristic velocity, :math:`U_\infty`, and a length scale,
:math:`L`. The pressure, density, and temperature are nondimensionalized
by :math:`P_\infty`, :math:`\rho_\infty`, and :math:`T_\infty`. The
enthalpy and energy are nondimensionalized by
:math:`C_{p,\infty} T_\infty`. Dimensionless variables are noted by
overbars. The dimensionless equations are:

.. math::
   
   {{\partial \bar{\rho}} \over {\partial \bar{t}}} 
       + {{\partial \bar{\rho} \bar{u}_j} 
          \over {\partial \bar{x}_j}} & = 0 , \\
   {{\partial \bar{\rho} \bar{u}_i} \over {\partial \bar{t}}} 
     + {{\partial \bar{\rho} \bar{u}_j \bar{u}_i}
       \over {\partial \bar{x}_j}} + {1 \over {\gamma {\rm Ma}^2}}
         {{\partial \bar{P}} \over {\partial \bar{x}_i}} & =
       {1 \over {\rm Re}}{{\partial \bar{\tau}_{ij}} 
               \over {\partial \bar{x}_j}} 
        +  {1 \over {\rm Fr}_i} \bar{\rho} , \\
   {{\partial \bar{\rho} \bar{h}} \over {\partial \bar{t}}} 
      + {{\partial \bar{\rho} \bar{u}_j \bar{h}}
       \over {\partial \bar{x}_j}}  
     & =
       - {1 \over {\rm Pr}} {1 \over {\rm Re}} 
                       {{\partial \bar{q}_j} \over {\partial \bar{x}_j}}
    + {{\gamma - 1} \over \gamma} {{\partial \bar{P}} \over  {\partial \bar{t}}} \\
      & + {{\gamma - 1} \over \gamma} {{{\rm Ma}^2} \over {\rm Re}}
     {{\partial \bar{u}_i \bar{\tau}_{ij}} \over {\partial \bar{x}_j}}  
    +  \bar{\rho} \bar{u}_i {{\gamma - 1} \over \gamma} 
               {{{\rm Ma}^2} \over {\rm Fr}_i} \nonumber \\
     & - {{\gamma - 1} \over 2} {\rm Ma}^2
    \left( {{\partial \bar{\rho} \bar{u}_k \bar{u}_k} \over {\partial \bar{t}}} 
       +   {{\partial \bar{\rho} \bar{u}_j \bar{u}_k \bar{u}_k} 
           \over {\partial \bar{x}_j}} \right) . \nonumber

The groupings of characteristic scaling terms are:

.. math::
   

   
   {\rm Re} & = {{\rho_\infty U_\infty L} \over {\mu_\infty}},
       \quad \quad \phantom{xxx} {\rm Reynolds number}, \\
   {\rm Pr} & =  {{C_{p,\infty} \mu_\infty} \over {k_\infty}},
       \quad \quad \phantom{xxx} {\rm Prandtl number}, \\
   {\rm Fr}_i & = {{u_\infty^2} \over {g_i L}},
       \quad \quad \phantom{xxxxxxi} {\rm Froude number}, \quad g_i \ne 0, \\
   {\rm Ma} & = \sqrt{{u^2_\infty} \over {\gamma R T_\infty /W}},
       \quad \quad {\rm Mach number},

where :math:`\gamma` is the ratio of specific heats.

For small Mach numbers, :math:`{\rm Ma} \ll 1`, the kinetic energy,
viscous work, and gravity work terms can be neglected in the energy
equation since those terms are scaled by the square of the Mach number.
The inverse of Mach number squared remains in the momentum equations,
suggesting singular behavior. In order to explore the singularity, the
pressure, velocity and temperature are expanded as asymptotic series in
terms of the parameter :math:`\epsilon`:

.. math::
   
      \bar{P} & = \bar{P}_0     + \bar{P}_1 \epsilon     + \bar{P}_2 \epsilon^2 \ldots \\
    \bar{u}_i & = \bar{u}_{i,0} + \bar{u}_{i,1} \epsilon + \bar{u}_{i,2} \epsilon^2 \ldots \\
      \bar{T} & = \bar{T}_0     + \bar{T}_1 \epsilon     + \bar{T}_2 \epsilon^2 \ldots

The zeroeth-order terms are collected together in each of the
equations. The form of the continuity equation stays the same. The
gradient of the pressure in the zeroeth-order momentum equations can
become singular since it is divided by the characteristic Mach number
squared. In order for the zeroeth-order momentum equations to remain
well-behaved, the spatial variation of the :math:`\bar{P}_0` term must
be zero. If the magnitude of the expansion parameter is selected to be
proportional to the square of the characteristic Mach number,
:math:`\epsilon = \gamma {\rm Ma}^2`, then the :math:`\bar{P}_1` term
can be included in the zeroeth-order momentum equation.

.. math::

   {1 \over {\gamma {\rm Ma}^2}}
      {{\partial \bar{P}} \over {\partial x_i}}  =
      {{\partial} \over {\partial x_i}} \left( {1 \over {\gamma {\rm Ma}^2}} \bar{P}_0
          + {\epsilon \over {\gamma {\rm Ma}^2}} \bar{P}_1 + \ldots \right) =
      {{\partial} \over {\partial x_i}} \left( \bar{P}_1 + \epsilon \bar{P}_2 + \ldots 
      \phantom{1 \over {\gamma {\rm Ma}^2}} \right)

The form of the energy equation remains the same, less the kinetic
energy, viscous work and gravity work terms. The :math:`P_0` term
remains in the energy equation as a time derivative. The low Mach number
equations are the zeroeth-order equations in the expansion including the
:math:`P_1` term in the momentum equations. The expansion results in two
different types of pressure and they are considered to be split into a
thermodynamic component and a dynamic component. The thermodynamic
pressure is constant in space, but can change in time. The thermodynamic
pressure is used in the equation of state. The dynamic pressure only
arises as a gradient term in the momentum equation and acts to enforce
continuity. The unsplit dimensional pressure is

.. math:: P = P_{th} + \gamma {\rm Ma}^2 P_1,

where the dynamic pressure, :math:`p=P-P_{th}`, is related to a
pressure coefficient

.. math:: \bar{P}_1 = {{P - P_{th}} \over {\rho_\infty u^2_\infty}} P_{th}.

The resulting unscaled low Mach number equations are:

.. math::

   {{\partial \rho} \over {\partial t}} + {{\partial \rho u_j} 
                                  \over {\partial x_j}} & = 0, \label{lmcon} \\
   {{\partial \rho u_i} \over {\partial t}} + {{\partial \rho u_j u_i}
       \over {\partial x_j}} + {{\partial P} \over {\partial x_i}} & =
       {{\partial \tau_{ij}} \over {\partial x_j}} 
        +  \left( \rho - \rho_{\circ} \right) g_i, \label{lmmom} \\
   {{\partial \rho h} \over {\partial t}} + {{\partial \rho u_j h}
       \over {\partial x_j}}  & =
       - {{\partial q_j} \over {\partial x_j}}
       + {{\partial P_{th}} \over {\partial t}}, \label{lmenrg}

where the ideal gas law becomes

.. math:: P_{th}  =  \rho {R \over W} T.

The hydrostatic pressure gradient has been subtracted from the momentum
equation, assuming an ambient density of :math:`\rho_{\circ}`. The
stress tensor and heat conduction remain the same as in the original
equations.
