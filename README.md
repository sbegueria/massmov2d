# massmov2d
A numerical model for mass movements over complex topography

massmov2d is a fluid mechanics model for simulating the kinematics (runout and deposition) of fluid-like mass movements such as landslides, rock slides, debris flows and mud flows, over a complex topography. The flow is modelled as a 2D continuum by using a depth-integrated form of the Navier-Stokes equations under the shallow water assumption (Saint-Venant equation). The flow behaviour is controlled by the resisting forces, for which a set of alternative rheological models can be used.

massmov2d has been implemented in the PCRaster dynamic modelling language. PCRaster is a free GIS package and dynamic modeling system, and provides standard tools for editing the input maps and visualizing the results through map animations, time series, etc. Besides that, being an open text script written in a easy to learn language allows the user of massmov2d to dig into the code and eventually incorporate new modeling concepts such as alternative rheologies.

## Usage

pcrcalc -f massmov_v091.mod [rheology] [rho] [yield_stress] [viscosity/Chèzy] [b_friction_angle]  [i_friction_angle] [fluid_rate]
[timesteps]

## Arguments

rheology An integer, defining the rheological law to be used. Use 1 for pure frictional and Voellmy frictional, and 2 for Bingham viscous or Coulomb viscous.
rho Density of the flow (kg.m-3). A scalar.
yield_stress Apparent yield stress (Pa), used only on viscous flows. A scalar.
viscosity/Chèzy Chézy roughness coefficient (m.s-2), for Voellmy rheology, or the dynamic viscosity (Pa.s), for viscous rheologies. A scalar.
b_friction_angle Angle of friction between the flow and the base (o). Used for both frictional and Coulomb viscous rheologies. An angle.
i_friction_angle Angle of internal friction (o). An angle.
fluid_rate Upward velocity of transition from solid to fluid of the landsliding mass (m.s-1). A scalar.
timesteps Number of time steps of the simulation (-). An integer. 1 of 5

## Details

massmov2d is a two-dimensional model of fluid-like mass movements kinematics over a complex topography. It is based on a depth-averaged form of the equation of motion for a fluid continuum, and controlled by rheology. Mass movements can be defined as gravity- driven flows of a mixture of rocks, soil, and water, with a behavior exhibiting properties of viscous and turbulent flows.

The flow is modeled as a 2D continuum by using a depth-integrated form of the Navier- Stokes equations under the shallow water assumption (Saint-Venant equation), an approach that has become classical for debris-flow modeling (e.g., Savage and Hutter, 1989; Mangeney-Castelnau et al., 2005). The model is based on a system of hyperbolic differential equations, its state variables being the flow thickness and velocity. The momentum conservation equation considers the net acceleration due to the weight of the material, internal pressure differences, flow resisting forces and the convective acceleration (i.e. the time rate of change in velocity due to change in position in the field). The flow behavior is controlled by the resisting forces, for which a set of alternative rheological models can be used. The constitutive equations are solved in a two-dimensional finite differences (Eulerian) mesh, using a two-step scheme with numerical regularization. The time step is one second, although internally the model used fractional time steps which vary depending on the flow characteristics, based on the Courant-Friedrichs-Levy condition. More details on the implementation of massmov2d can be found on Beguería et al. (2009).

The model has been implemented as a PCRaster script. PCRaster is a free GIS package and dynamic modeling system, and provides standard tools for editing the input maps and visualizing the results through map animations, time series, etc. Besides that, being an open text script written in a easy to learn language allows the user of massmov2d to dig into the code and eventually incorporate new modeling concepts such as alternative rheologies. The PCRaster software and manuals can be downloaded from http://pcraster.geo.uu.nl/. It is recommended that the latest version plus updates are installed, as well as the aguila visualizing software. PCRaster includes a set of functions allowing basic vector calculus operations since version 3.0 . This functions are needed in order to operate massmov2d, so please make sure that you have downloaded the latest version. As of March 2011 version 3.0 is a beta release, and you can find it at http://mailman.geo.uu.nl/pipermail/pcraster-info/2010-January/000403.html.

A number of input files must be provided to massmov2d. These maps need to be stored in the same directory than the script, and their names must be identical to the ones listed here:

dem.map Surface elevation (m). A scalar map. 2 of 5  
h_ini.map Landslide initial thickness, measured normal to the surface (m). A scalar map.
outlet.map Open boundaries of the computation domain (if any). A boolean map.
fluiddist.map Distance to the toe of the landslide (m). A scalar map.

The surface elevation map is a common digital elevation model (DEM), from which the initial landslide body has been already subtracted. This map defines not only the basal topography for the flow, but also the computation domain, so any missing value pixels in this map will be kept out of the calculations.

The landslide initial thickness map is measured in the direction normal to the surface, and not in the Cartesian vertical direction; it can be obtained as h = h_vertical * cos(alpha), alpha being the slope gradient.

The outlet map is usually just a map of (boolean) zeros, in which case the boundaries of the spatial domain are considered closed (the flow can not escape from the spatial domain). In some cases, though, it is convenient to allow the flow to escape from the computation domain, for example if there is a river at the bottom of the hill slope which is able to remove all the material that enters the stream. In this case, an open boundary can be indicated as a line of (boolean) ones in the outlet map.

The fluidization distance map indicates the distance to the toe of the landslide. The toe of the landslide, which is defined as having zero distance, is considered to be fluid at the beginning of the simulation. From that area, fluidization advances at a time rate defined by the fluid_rate parameter.

## Output

The model produces a series of standard outputs, and many others can be produced by uncommenting the appropriate lines of code in the script. Reporting is done for each second of the simulation time, although other reporting frequencies can be set by modifying the rep_a and rep_b variables in the timer section of the script. The standard outputs are the following files:

h0000000.xxx The flow thickness, in the direction normal to the surface (m). A scalar map time series.
vel00000.xxx The flow velocity, in the direction parallel to the surface (m.s-1). A scalar map time series.
v091_x1_x2_x3_x4_x5_x6_x7_x8.tss
A dummy time series containing just zeroes. The name of the file provides information on the simulation parameters, such as the script version number and the input values passed to the model (x1 to x8).

## Disclaimer

This script is provided 'as is' and any express or implied warranties, including, but not limited to, the implied warranties of fitness for a particular purpose are disclaimed. In no event shall the author be liable for any direct, indirect, incidental, special, exemplary or consequential damages arising in any way out of the use of this software.

## Examples

A particular simulation scenario is defined by the basal topography (dem.map) and the initial landslide body thickness (h_ini.map), plus the parameters defining the characteristics of the flow which are passed to MassMov2D as command-line parameters. For example, given flow having a density of 1850 kg.m-3 and angle of internal friction 15o, considering a fluidization rate of 10 m.s-1 from the toe of the landslide and a total simulation time of 250 seconds, simulations could be launched for the different rheological laws implemented in massmov2d v0.9.1 by:

• Pure frictional: pcrcalc -f massmov_v091.mod 1 1850 0 0 11.5 15 10 250
• Voellmy fictional: pcrcalc -f massmov_v091.mod 1 1850 0 300 11.5 15 10
250
• Bingham viscous: pcrcalc -f massmov_v091.mod 2 1850 50 100 0 15 10 250
• Coulomb viscous: pcrcalc -f massmov_v091.mod 2 1850 50 100 11.5 15 10 250

The output map series and time series generated can be viewed with aguila, using the following command:
aguila -s[1,250,1] h vel

## References

* Beguería, S., Van Asch, T.W.J., Malet, J.P. & Gröndahl S., 'A GIS based numerical model for simulating the kinematics of mud and debris flows over complex terrain', Nat Hazards Earth Syst Sci 9, 1897-1909, 2009.
4 of 5
* Mangeney-Castelnau, A., Bouchut, F., Vilotte, J.-P., Lajeunesse, E., Aubertin, A. amd Pirulli, M., 'On the use of Saint Venant equations to simulate the spreading of a granular mass', J Geophys Res 110(B9), B09103, 2005.
* Savage, S.B. and Hutter, K., 'The motion of a finite mass of granular material down a rough incline', J Fluid Mech 199, 177–215, 1989.
