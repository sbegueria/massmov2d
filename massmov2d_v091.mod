# _Important preliminar information_

# Introduction
# MassMov2D is a simple fluid mechanics model for simulating the kinematics
# (runout and deposition) of fluid-like mass movements such as landslides,
# rock slides, debris flows and mud flows, over a complex topography. The
#flow is modelled as a 2D continuum by using a depth-integrated form of the
# Navier-Stokes equations under the shallow water assumption (Saint-Venant
# equation). The model is formulated as a system of hyperbolic differential
#equations, its state variables being the flow thickness and velocity. The
# momentum conservation equation considers the net acceleration due to the
# weight of the material, internal pressure differences, flow resisting
# forces and the convective acceleration (i.e. the time rate of change in
# velocity due to changes in position in the field). The flow behaviour is
# controlled by the resisting forces, for which a set of alternative rheolo-
# gical models can be used. The constitutive equations are solved in a two-
# dimensional finite differences (Eulerian) mesh, using a two-step scheme
# with numerical regularization. The time step is one second, although inter-
# nally the model used fractional time steps which vary depending on the
# flow characteristics, based on the Courant-Friedrichs-Levy condition.
# More details on the implementation of MassMov2D can be found on:
# Beguería S., Van Asch T.W.J., Malet J.P., Grondahl S. (2009), A GIS-based
# numerical model for simulating the kinematics of mud, debris flows over
# complex terrain. Natural Hazards and Earth System Sciences 9, 1897-1909.

# Use
# Download this script, documentation and example data files from:
# http://digital.csic.es/handle/10261/11804.
#
# Install the last version (3.0) of PCRaster:
# http://mailman.geo.uu.nl/pipermail/pcraster-info/2010-January/000403.html.
# See also: http://pcraster.geo.uu.nl.
#
# To run the model use the following syntax:
#
# pcrcalc -f massmov_v091.mod [rheology] [rho] [yield_stress] [viscosity]
# [basal_friction_angle] [internal_friction_angle] [fluidization_rate]
# [timesteps]
#
# For example: pcrcalc -f MassMov_v091.mod 1 2000 0 300 11.5 15 10 250

# Disclaimer
# This script is provided by the author 'as is' and any express or implied
# warranties, including, but not limited to, the implied warranties of fitness
# for a particular purpose are disclaimed. In no event shall the author be
# liable for any direct, indirect, incidental, special, exemplary or conse-
# quential damages arising in any way out of the use of this software.

# License
# MassMov2D is copyright © 2006-2011 of Santiago Begueria-Portugues.
# Some rights reserved
# BY NC SA - This work is distributed under an 'Attribution Non-Commercial
# Share Alike' license. This means that anyone is allowed to use it and build
# upon it work non-commercially, as long as credit is given to the author and
# the new creations are licensed under the identical terms. For more informa-
# tion refer to: http://creativecommons.org/licenses/by-nc-sa/3.0/

# Contact
# santiago.begueria@csic.es


binding

 # Input maps
 DEM          = "../dem.map";	# Elevation map (sliding surface, m)
 H_ini        = "../h_ini.map";	# Landslide body thickness (m, vertical)
 Outlet       = "../outlet.map";# Domain boundary (binary)
 Fluiddist    = "../dist.map";	# Distance to fluidized toe (m)

 # Mud rheology
 Rheol        = scalar($1);		# Rheology type (1 = frictional, 2 = viscous)
 Rho          = scalar($2);		# Density of the df (kg.m-3)
 YieldStress  = scalar($3);		# Apparent yield stress (Pa), for Bingham
 Visco        = scalar($4);		# Bingham viscosity (Pa.s), for Bingham and Coulomb
 Chezy        = scalar($4);		# Chézy roughness coefficient (m.s-2), for Voellmy
 HBExp   	  = 1;				# Herschel-Bulkley exponent (=1 for Bingham)
 Friction     = scalar($5);		# Angle of friction basal (º)
 IntFriction  = scalar($6);		# Angle of internal friction (º)

 # Inlet boundary condition:
 # Only for landslides: rate of fluidization of the failed mass
 Fluidrate    = scalar($7);		# Transition from solid to fluid (m.s-1) 

 # Timestep control
 dT           = 1;				# Timeslice (seconds)
 Timesteps    = $8;		        # Number of timesteps of the simulation

 # Numerical stability control
 CFLlimsup    = 0.5;			# Higher value of the Courant-Friedrichs-Levy
 CFLliminf    = 0.3;			# Lower value of the Courant-Friedrichs-Levy
 MinNLoops    = 1;				# Minimum number of internal loops
 MaxNLoops    = 124;			# Maximum number of internal loops
 InitialLoops = 1;				# Initial number of internal loops
 LaxFactor    = 0.5;			# Low pass filtering (0 = no filtering)
 dTLeap   	  = 0.3;			# Fraction of dT for time extrapolation
 verysmall    = 0.01;			# threshold to avoid div by very small values

 # Other constants
 Grav 		  = 9.8;			# Gravity acceleration, Earth (m.s-2)
 #Grav 		  = 3.69;			# Gravity acceleration, Mars (m.s-2)



areamap

# calculation domain
 DEM;


timer

 1 Timesteps 1;					# repetitions: from t=A to t=B by incr=1
 rep_a = 1, 1 + 1..endtime;		# reporting intervals (at t=1, then each
 rep_b = 1, 25 + 25..endtime;	# n+n until endtime)
 end   = endtime;


initial

 # Fixed terms (constants)
 Sqrt2 = sqrt(2);                                 # just the square root of two
 CL = celllength();                               # cell length (plain)
 TanPhi = tan(Friction);                          # tangent of the friction angle

 # Fixed maps
 # calculation domain
 Clone = if (defined(DEM), boolean(1));           # the computation domain
 Zeros = if(Clone, scalar(0));                    # a map full of zeros

 # basal surface gradient and vertical curvature;
 # central diffs. within the DEM domain, upwind diffs. at the boundaries
 Z = DEM;
 # gradient (-x and -y components) and laplacian of Z
 dZ\dx = gradx(Z);                                # gradient in x: <0 to the E
 dZ\dy = grady(Z);                                # gradient in y: <0 to the N
 d2Z\dx2 = lax(laplacian(Z), 1);		  		  # surface curvature
 d2Z\dy2 = d2Z\dx2;
 # slope angles, negative towards E and N
 SinAlfa = sin(atan(slope(Z)));
 CosAlfa = cos(atan(slope(Z)));
 SinAlfa_x = sin(atan(dZ\dx));
 SinAlfa_y = sin(atan(dZ\dy));
 CosAlfa_x = cos(atan(dZ\dx));
 CosAlfa_y = cos(atan(dZ\dy));
 CL_x = CL / CosAlfa_x;
 CL_y = CL / CosAlfa_y;
 CA = CL*CL;                                      # cell area
 #CA = CL_x * CL_y;

 # earth pressure coefficient (Koch, 1998)
 SinPhi = sin(IntFriction);
 K_act = (1-SinPhi) / (1+SinPhi);                 # extension
 K_pas = (1+SinPhi) / (1-SinPhi);                 # compression
 #K_act = 0.6;
 #K_pas = 1.5;
 K = if(Clone, scalar(1));

 # Initial values
 H  = if(Clone, scalar(0));                       # body thickness
 i  = 0;    									  # a counter
 n  = 0;                                 		  # another counter
 U = if(Clone, scalar(0));                        # momentum, x-wise
 V = if(Clone, scalar(0));                        # momentum, y-wise
 H_loop = H;                                      # H, for internal loop
 U_loop = U;                                      # ''
 V_loop = V;                                      # ''
 H_inlt = H;                                      # H at the ghost inlet cells
 V_inlt = V;
 U_inlt = U;
 Vel_inlt = 0;
 Vol_expect = maptotal(H/CosAlfa * CA);			  # initial volume (m3)
 Incr_volexp = 0;
 Vol_gone   = 0;                                  # volume gone through outlet (m3)
 Loops = InitialLoops;                            # number of internal loops
 exit = boolean(0);
 T_x = 0;
 T_y = 0;
 T_x_b = 0;
 T_y_b = 0;


dynamic

 # boundary conditions (mass and momentum source and sinks)
 # only for landslides; comment out if not applicable
 H_fluid = if(Fluiddist le i*dT*Fluidrate and Fluiddist gt (i-1)*dT*Fluidrate,
	H_ini, Zeros) ; # newly fluidized mass
 H_land  = if(Fluiddist gt i*dT*Fluidrate,
	H_ini, Zeros); # remaining (solid) landslide mass
 H = H + H_fluid;
 Incr_volexp = maptotal(H_fluid/CosAlfa * CA);
 # do not comment the next line
 i = i + 1; # }}}


 repeat {        # main loop

   # initialize variables in the internal loop
   n = 0;
   exit = boolean(0);
   CFLmax = 0;
   Vel_aver = 0;
   H_loop = H;
   U_loop = U;
   V_loop = V;
   Vol_exp_loop = Vol_expect;
   Vol_gone_loop = Vol_gone;

   repeat {      # internal loop: subdivide dT into a number
                 # of smaller timesteps to ensure stability


     ## 1. DETERMINE FLOW VELOCITY (U,V) AT t+dt
     ## d/dt(u) + cosAlfa_x u d/dx(u) + cosAlfa_y v d/dy(u) = - G_x - P_x + F_x
     ## d/dt(v) + cosAlfa_y v d/dx(v) + cosAlfa_x u d/dy(v) = - G_y - P_y + F_y
     ## G, P and T conform the local or time acceleration term, representing the
     ## time rate of change at a fixed position for a non-steady flow.
     ## I is the convective acceleration, representing the time rate of change
     ## due to change in position in the field for a non-uniform flow.
     ## (Following Mangeney et al. 2007)

     # define debris flow domain
     is_df = if(H_loop gt verysmall, boolean(1), boolean(0));

     # calculate acceleration sources at time t (G, P, I, T)
     # G, gravitational acceleration = cosAlfa_x g tanAlfa_x = g sinAlfa_X (m.s-2)
     G_x  = if(is_df, -Grav*SinAlfa_x, 0);
     G_y  = if(is_df, -Grav*SinAlfa_y, 0);

     # P, pressure acceleration = cosAlfa_x k d/dx(g cosAlfa_x h) (m.s-2)
     K = 0.5 * (K + if(gradx(U_loop) + grady(V_loop) ge 0, K_act, K_pas));
     K = cover(lax(K, LaxFactor), Zeros);	# earth pressure ratio
     P_x  = if(is_df, -K*CosAlfa_x*gradx(Grav*(H_loop+H_land)*CosAlfa_x), 0);
     P_y  = if(is_df, -K*CosAlfa_y*grady(Grav*(H_loop+H_land)*CosAlfa_y), 0);

     # I, convective acceleration = cosAlfa_x u du/dx + cosAlfa_y v du/dy (m2.s-2)
     dU\dx = CosAlfa_x*gradx(abs(U_loop+U_inlt));
     dU\dy = CosAlfa_y*grady(abs(U_loop+U_inlt));
     dV\dx = CosAlfa_x*gradx(abs(V_loop+V_inlt));
     dV\dy = CosAlfa_y*grady(abs(V_loop+V_inlt));
     I_x = -(U_loop*dU\dx + V_loop*dU\dy);
     I_y = -(U_loop*dV\dx + V_loop*dV\dy);

     # T, negative acceleration due to basal shear stress = cosAlfa g S q_x (m.s-2)
	 Vel = sqrt(U_loop*U_loop + V_loop*V_loop);
	 # Purely frictional
     T = if(Rheol eq 1 and Chezy le verysmall and H_loop gt verysmall,
 		TanPhi, 5*TanPhi);
	 # Voellmy (turbulent frictional)
     T = if(Rheol eq 1 and Chezy gt verysmall and H_loop gt verysmall,
 		TanPhi + Vel*Vel/Chezy/H_loop, 5*TanPhi);
	 # Bingham/Herschel-Bulkley/Coulomb viscous (viscoplastic, no yield stress)
     T = if(Rheol eq 2 and YieldStress le verysmall, if(H_loop gt verysmall,
		TanPhi + 3*Visco*Vel**HBExp/Rho/H_loop/H_loop, 5*TanPhi), T);
	 # Bingham/Herschel-Bulkley/Coulomb viscous (viscoplastic, with yield stress)
     T = if(Rheol eq 2 and YieldStress gt verysmall, if(H_loop gt verysmall,
		TanPhi + (1.5*YieldStress + 3*Visco*Vel**HBExp/H_loop)/(Rho*H_loop),
		5*YieldStress), T);
	 # Compute x- and y- components
     T_x = if(is_df and Vel gt verysmall, abs(U_loop)/Vel, 1) * Grav*CosAlfa_x*T;
     T_y = if(is_df and Vel gt verysmall, abs(V_loop)/Vel, 1) * Grav*CosAlfa_y*T;

	 # Estimate flow at t = t + ddt
     ddt = dTLeap * dT/Loops; # fractional interval (s)
     # U
	 GP_x = G_x + P_x;
     U_loop_a = cover(lax(if(is_df, U_loop), LaxFactor), Zeros);
     U_loop_b = if(U_loop_a gt 0, max(0, U_loop_a + ddt * (GP_x + I_x - T_x)),
	 	if(U_loop_a lt 0, min(U_loop_a + ddt * (GP_x - I_x + T_x), 0), 0));
     U_loop_b = if(U_loop_a eq 0 and GP_x gt 0, max(0, ddt * (GP_x + I_x - T_x)),
	 	if(U_loop_a eq 0 and GP_x lt 0, min(ddt * (GP_x - I_x+T_x), 0), U_loop_b));
     # V
	 GP_y = G_y + P_y;
     V_loop_a = cover(lax(if(is_df, V_loop), LaxFactor), Zeros);
     V_loop_b = if(V_loop_a gt 0, max(0, V_loop_a + ddt * (GP_y + I_y - T_y)),
	 	if(V_loop_a lt 0, min(V_loop_a + ddt * (GP_y - I_y + T_y), 0), 0));
     V_loop_b = if(V_loop_a eq 0 and GP_y gt 0, max(0, ddt * (GP_y + I_y - T_y)),
	 	if(V_loop_a eq 0 and GP_y lt 0, min(ddt * (GP_y -I_y+T_y), 0), V_loop_b));

     # T, negative acceleration due to basal shear stress = cosAlfa g S q_x (m.s-2)
	 Vel_b = sqrt(U_loop_b*U_loop_b + V_loop_b*V_loop_b);
	 # Purely frictional
     T_b = T;
	 # Voellmy (turbulent frictional)
     T_b = if(Rheol eq 1 and Chezy gt verysmall and H_loop gt verysmall,
 		TanPhi + Vel_b*Vel_b/Chezy/H_loop, 5*TanPhi);
	 # Bingham/Herschel-Bulkley/Coulomb viscous (viscoplastic, no yield stress)
     T_b = if(Rheol eq 2 and YieldStress le verysmall, if(H_loop gt verysmall,
		TanPhi + 3*Visco*Vel_b**HBExp/Rho/H_loop/H_loop, 5*TanPhi), T);
	 # Bingham/Herschel-Bulkley/Coulomb viscous (viscoplastic, with yield stress)
     T_b = if(Rheol eq 2 and YieldStress gt verysmall, if(H_loop gt verysmall,
		TanPhi + (1.5*YieldStress + 3*Visco*Vel_b**HBExp/H_loop)/(Rho*H_loop),
		5*YieldStress), T);
	 # Compute x- and y- components
     T_x_b = if(is_df and Vel_b gt verysmall, abs(U_loop_b)/Vel_b, 1) *
		Grav*CosAlfa_x*T_b;
     T_y_b = if(is_df and Vel_b gt verysmall, abs(V_loop_b)/Vel_b, 1) *
		Grav*CosAlfa_y*T_b;

     # estimate flow at t = t + dt (using T at time t + ddt)
     # U
     U_loop = if(U_loop_a gt 0, max(0, U_loop_a + dT/Loops * (GP_x + I_x - T_x_b)),
      	if(U_loop_a lt 0, min(U_loop_a + dT/Loops * (GP_x - I_x + T_x_b), 0), 0));
     U_loop = if(U_loop_a eq 0 and GP_x gt 0, max(0, dT/Loops * (GP_x +I_x-T_x_b)),
   		if(U_loop_a eq 0 and GP_x lt 0, min(dT/Loops * (GP_x - I_x + T_x_b), 0),
		U_loop));
     # H
     V_loop = if(V_loop_a gt 0, max(0, V_loop_a + dT/Loops * (GP_y + I_y - T_y_b)),
       if(V_loop_a lt 0, min(V_loop_a + dT/Loops * (GP_y - I_y + T_y_b), 0), 0));
     V_loop = if(V_loop_a eq 0 and GP_y gt 0, max(0, dT/Loops * (GP_y +I_y-T_y_b)),
       if(V_loop_a eq 0 and GP_y lt 0, min(dT/Loops * (GP_y - I_y + T_y_b), 0),
	   V_loop));



     ## 2. DETERMINE D.F. THICKNESS (H) AT t = t + dt
     ## dh/dt + cosAlfa_x d/dx(hu) + cosAlfa_y d/dy(hv) = 0
	 ## dh = -dt * (cos Alfa_x d/dx(hu) + cosAlfa_y d/dy(hv))
     ## (this is: dh = -dT *  div(momentum_field))

     # rate of d.f. thickness variation (m.s-1); dH\dT > 0 net gain of material
     HU_loop = if(is_df, H_loop*U_loop, 0);
     HV_loop = if(is_df, H_loop*V_loop, 0);

	 dH\dT = -CosAlfa_x * (shift0(HU_loop, 0, -1) - shift0(HU_loop, 0, 1)) /
		(2*CL_x) - CosAlfa_y * (shift0(HV_loop, 1,  0) - shift0(HV_loop,-1, 0)) /
		(2*CL_y);

     # update body height
     H_loop_a = cover(lax(if(is_df, H_loop), LaxFactor), Zeros);
     H_loop_a = if(dH\dT eq 0, H_loop, H_loop_a); # avoid dispersion at d.f.borders
	 H_loop   = H_loop_a - dT/Loops * dH\dT;      # if there is no change
     # ~0.001 mass balance error in our numerical scheme;
     # hence, we need to correct for m. b. e. in the next lines
     Vol_gone = Vol_gone + maptotal(if(Outlet, H_loop/CosAlfa * CA));
     H_loop   = if(Outlet, 0, H_loop);
     Vol_actual = maptotal(H_loop/CosAlfa * CA) + Vol_gone;
     Vol_exp_loop = Vol_exp_loop + Incr_volexp/Loops;
     MassBalError = if(Vol_exp_loop gt verysmall, Vol_actual/Vol_exp_loop, 1);
     H_loop = if(MassBalError gt verysmall, H_loop/MassBalError, H_loop);



     ## 3. NUMERICAL STABILITY ANALYSIS
     ## solution is stable if CFL (Courant-Friedrichs-Levy) condition < 1

     # CFL
     CardVel = if(is_df,
		max(abs(CosAlfa_x*U_loop), abs(CosAlfa_y*V_loop)), 0);
	 #CFL = dT/Loops * Sqrt2 * mapmaximum(CardVel) / CL_x;
	 #CFL_inlt = dT/Loops * Sqrt2 * Vel_inlt / CL_y;
     CFL = dT/Loops * Sqrt2 * mapmaximum(CardVel) / CL;
     CFL_inlt = dT/Loops * Sqrt2 * Vel_inlt / CL;
     CFLmax = max(CFLmax, CFL, CFL_inlt);

     # update number of internal loops
     IsUnstable = if(CFLmax ge CFLlimsup,	# TRUE if the solution becomes
       boolean(1), boolean(0));				# unstable, FALSE if stable
     n = if(IsUnstable, n, n+1);			# prevents exiting if unst. & n eq Loops

	} until n eq Loops or IsUnstable;

   Loops = if (IsUnstable, Loops+1, Loops);	# try with more loops if unstable
   #Loops = if (IsUnstable, roundoff(Loops*1.5), Loops);


 # } exit the main loop
 } until not IsUnstable or Loops gt MaxNLoops;
 # some dummy identities needed to avoid an exception error due to a bug in pcrcalc!
 # see http://pcraster.geo.uu.nl/documentation/manual_updates/repeatUntil.html
 H = H;
 U = U;
 V = V;
 H_inlt = H_inlt;
 V_inlt = V_inlt;
 U_inlt = U_inlt;
 Vel_inlt = Vel_inlt;
 H_land = H_land;
 Vol_gone = Vol_gone;
 Vol_expect = Vol_expect;
 IsUnstable = IsUnstable;
 CFLmax = CFLmax;
 Incr_volexp = Incr_volexp;
 Loops = Loops;



 ## 4. UPDATE VARIABLES FOR TIMESTEP MARCHING

 # conservative variables
 H = H_loop;
 U = U_loop;
 V = V_loop;
 Vol_expect = Vol_exp_loop;
 Vol_gone = Vol_gone_loop;
 Loops = if(CFLmax le CFLliminf,			# don't let CFL drop below 0.5 to
    max(Loops-1, MinNLoops), Loops);		# control numerical dispersion
   #max(roundoff(Loops/1.5), MinNLoops), Loops);
 # }}}


 ## 5. MODEL OUTPUTS
 # Uncomment the lines according to what outputs are needed. Note that hard drive
 # accessing for reading / writing takes a substantial amount of time so enabling
 # many reports will increase the computation time.

 Vel = if(H gt verysmall,						# velocity modulus (m.s-1)
	sqrt(U_loop*U_loop + V_loop*V_loop), 0);


 # maps at endtime
 #report (end) s${4}_${5}_h.map = if(H gt 0,H);	# body thickness (m)
 #report (end) s${4}_${5}_b.map = if(H gt 0,scalar(1));	# body thickness (m)
 #report (end) s${4}_${5}_v.map = Vel;			# velocity last loop (m.s-1)
 #report (end) hh = if(H gt 0, H);				# body thickness (m)
 #report (end) vv = Vel;						# velocity last loop (m.s-1)
 #report (end) u = U;
 #report (end) v = V;

 # map series
 report (rep_a) h = if(H+H_land gt 0, H+H_land);# body thickness, no zeros (m)
 #report (rep_a) h = if(H+H_inlt gt 0,H+H_inlt);# body thickness, no zeros (m)
 #report (rep_a) h = H;                     	# body thickness (m)
 #report (rep_a) zh = Z + 1*H+H_land;           # free surface (with vert.exag.)
 report (rep_a) vel = if(H+H_land gt 0, Vel);   # flux velocity (m.s-1)
 #report (rep_a) q = if(H gt 0, H*CL*Vel);
 #report (rep_a) q;
 #report (rep_a) uh = if(UH ne 0, UH);
 #report (rep_a) vh = if (VH ne 0, VH);
 #report (rep_a) u = if(UH/H ne 0, UH/H);
 #report (rep_a) v = if(VH/H ne 0, VH/H);
 #report (rep_a) k = K;
 #report (rep_a) gp = sqrt(GP_x*GP_x + GP_y*GP_y);
 #report (rep_a) i  = sqrt(I_x*I_x + I_y*I_y);
 #report (rep_a) t = sqrt(T_x*T_x + T_y*T_y);
 #report (rep_a) gx = if(G_x ne 0, G_x);
 #report (rep_a) px = if(P_x ne 0, P_x);
 #report (rep_a) gpix = GPI_x;
 #report (rep_a) ix = if(I_x ne 0, I_x);
 #report (rep_a) tx = if(T_x ne 0, T_x);
 #report (rep_a) txb = if(T_x_b ne 0, T_x_b);
 #report (rep_a) gy = if(G_y ne 0, G_y);
 #report (rep_a) py = if(P_y ne 0, P_y);
 #report (rep_a) gpy = GP_y;
 #report (rep_a) iy = I_y;
 #report (rep_a) ty = T_y;

 #time series
 #report maxh.tss    = mapmaximum(H);       # maximum d.f. velocity (m)
 #report maxv.tss    = mapmaximum(Vel);     # maximum d.f. velocity (m.-1)
 #report maxh.tss    = mapmaximum(CosAlfa*H);# maximum d.f. velocity (m)
 #report maxv.tss    = mapmaximum(Vel);     # maximum d.f. velocity (m.-1)
 #report maxvh.tss   = mapmaximum(Vel*H);   # maximum d.f. momentum (m.m.s-1)
 #report maxq.tss    = mapmaximum(q);       # maximum d.f. discharge (m3.s-1)
 #report cfl.tss     = CFLmax;              # maximum Courant condition value
 #report loops.tss   = n;                   # nº intern loops whithin timestep
 #report vol.tss	 = Vol_actual;          # total volume (m3)
 #report volexp.tss	 = Vol_expect;          # expected volume (m3)
 #report massbal.tss = MassBalError;        # mass balance error (m3.m-3)
 #report massout.tss = Vol_gone;            # mass gone through outlet (m3)
 #report H.tss = timeoutput(SamplePnts,H);  # d.f. thickness at sample points
 #report V.tss = timeoutput(SamplePnts,     # d.f. velocity at sample points
 #  Vel);
 #report frc.tss = timeoutput(SamplePnts,   # force exerted by the flow
 # if(H gt verysmall, Rho*H*CL *            # against a rigid obstacle (kN)
 # (CL*CosAlfa_y + UH*VH/(2*Grav*H*H)), 0));
 #report inlth.tss = H_inlt;
 #report inltv.tss = V_inlt;
 # model params
 #report v091_${1}_${2}_${3}_${4}_${5}_${6}_${7}_${8}.tss =
 # boolean(0);                               # just the simulation params.


 # input maps
 #report (end) dem.map    = Z;
 #report (end) h_ini.map  = H_ini;
 #report (end) outlet.map = Outlet; # }}}

