#pragma once
#include <string>

//! defines.h has been cut down in favor of defining the some
//! variables here in parameters.h
//! parameters.cpp has also been reduced since all the parameters defined there
//! were constants anyway. Only functions which have been relocated from elsewhere are
//! now implemented there, where all constants are now defined here for compiler access,
//! which improves performance.

//! Some defines.h definitions come first for use in variables previously from parameters.cpp
//! Variables you set begin with the GENERAL REMARKS below

//! Set these with #defines because we're unable to fully escape the proprocessor right now
//! These are used in header files, where we have to use #if
#define USE_NEUTRAL_SPECIES_AS_FIELD	true
#define USE_ION_PRODUCTION_AS_FIELD	true
#define USE_DUST_SPECIES_AS_FIELD	false

//! But I'm going to use constexpr where I can in the hopes that one day we can
//! These are used in functions, where we can use 'if constexpr'
constexpr bool use_neutral_species_as_field = USE_NEUTRAL_SPECIES_AS_FIELD;
constexpr bool use_ion_production_as_field = USE_ION_PRODUCTION_AS_FIELD;
constexpr bool use_dust_species_as_field = USE_DUST_SPECIES_AS_FIELD;

//! Turn atmosphere on and off
constexpr bool use_atmosphere = true;

//! Turn dipole on and off
//! Be sure and set USE_CFBG_FIELD below accordingly
constexpr bool use_dipole = true;

//! Defines are needed for array sizing, regular variables (below) are used elsewhere
//! number of neutral species
//! define needed for chemical reactions array
constexpr int NUM_NEUTRAL_SPECIES = 2 * use_atmosphere;

//! number of species that are treated as particles
constexpr int NUM_PARTICLE_SPECIES = 1 + NUM_NEUTRAL_SPECIES;

//! number of ion species, may be different from NUM_PARTICLE_SPECIES, see above
//! likely the same for Ganymede
constexpr int NUM_CHARGED_SPECIES = 1 + NUM_NEUTRAL_SPECIES;

//! number of ion production fields for external ionosphere file
constexpr int num_ion_prod_fields = NUM_NEUTRAL_SPECIES;

//! specify number of fields to average
constexpr int num_average_fields = 0;

//! specify number of extern density fields for Ion injection
constexpr int num_externRhoVelocityFields = 0;


//! -------------- FIELD RELATED -----------------------------------------

//! use a dipole field as background
#define USE_CFBG_BFIELD		true

//! use a better smooth_Field method
#define use_new_smooth_Field_method  true

constexpr int NUM_SUB_CYCLE = 10;


//! switch off/on physical parts of the E/BField Equation
#define CONVEC_TERM		true
#define HALL_TERM		true
#define gradPE_TERM		false

//! switch on/off use eta in Leap Frog method.
//! Must be false if advance_obstacle_B below is set true
#define ETA_TERM_BField	false

//! switch on/off use eta in EField method.
#define ETA_TERM_EField	false

//! see up top about the prepocessor
#define NONADIABATIC_GRADPE_TERM			false
#define nonadiabatic_gradPE_COLLISION_TERM	false
#define nonadiabatic_gradPE_SOURCE_TERM		false
#define nonadiabatic_gradPE_DRAIN_TERM		false
constexpr bool nonadiabatic_gradPE_TERM = NONADIABATIC_GRADPE_TERM;

//! Set nonadiabatic_gradPE_TERM_derivation_form to
//! 0 : non-conservative form
//! 1 : conservative form
//! 2 : conservative form using extra field (id_scratch_scalar) for pow(PE, 1./kappa_electron)
#define nonadiabatic_gradPE_TERM_derivation_form 2

#define nonadiabatic_gradPE_TERM_set_divUe_to_zero_DEBUG	false

//! Set nonadiabatic_gradPE_TERM_curl_derivation_form to
//! 0 : Compute rot( 1/rho * grad(Pe) )
//! 1 : Compute - rot( Pe * grad(1/rho) )
//! 2 : Compute grad(1/rho) x grad(Pe)
//! If rot(grad(...)) = 0, all terms are equivalent
#define nonadiabatic_gradPE_TERM_curl_derivation_form	1


//! set if you want to smooth PE similar to E and B
constexpr bool nonadiabatic_gradPE_TERM_smooth_PE = true;

//! define this to ensure PE >= 0
#define nonadiabatic_gradPE_TERM_ensure_PE_not_negative		true

//! Only use in case of Debug
//! If you use conservative form from above
//! set nonadiabatic_gradPE_TERM_adiabatic_term_only to true
//! to use only the adiabatic term in the time advancement
#define nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only	false
#define nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only	false

//! Set nonadiabatic_gradPE_initial_condition to
//! 1 if you want to have homogenous electron pressure as initial value
//!   the initial value for the electron pressure will be set via Electron_Betas[0]
//! 0 if you want to have homogenous electron temperature as initial value
//!   the initial value for the electron pressure will be computed
//!   using the ion density and Electron_Betas[0]
//!   PE[u_v_w] = RHO[u_v_w] * Electron_Betas[0]
//! NOTE: good old AIKEF with adiabatic electron pressure assumes
//!       homogeous electron pressure defined by Electron_Betas
//!       as initial value
constexpr int nonadiabatic_gradPE_initial_condition = 0;

//! If you use nonadiabatic_gradPE_TERM
//! you have to decide, whether Pe shall be advanced
//! together with B or in an extra loop using 
//! the explicit midpoint method.
//! In the first case, define nonadiabatic_gradPE_TERM_advance_Pe_with_B
#define nonadiabatic_gradPE_TERM_advance_Pe_with_B	true



//! If kappa_electron is 2.0, you should set this to true
//! to get better performance, especially with vectorclass defined below
#define kappa_electron_is_2 true

//!------------------------------------------------------------------------------!//
//! -----------------------------------------------------------------------------!//
//!------------------------------------------------------------------------------!//

//! For gathering ... 
constexpr int ALL_SPECIES = -1;

constexpr bool noVTH2  = false;
constexpr bool getVTH2 = true;

//! using short for INT32 does not speed up the code
//! (maybe 1% at most)

//!------------- Geometry Field Flags -------------------//
constexpr int FILE_TYPE_INFO  	= 200;
constexpr int FILE_TYPE_BLOCKS 	= 201;

//!------------- Field IDs -------------------//

constexpr int id_notDefined    = -1;


//! -------------- PARTICLE RELATED --------------------------------------
//! JUST FOR STATISTICS:
//! FOR SOME REASON CAN SLOW SIGNIFICANTLY DOWN THE CODE !!!
constexpr bool ESTIMATE_FASTEST_PARTICLE = true;

//! MAX_LEVEL_DEF
//! to determine the fastest particle per level set min. MAX LEVEL 
constexpr int MAX_LEVEL_DEF  = 7;

//! Delete particle at very high velocities which
//! do not fulfill the courant criteria and thereby
//! will make the move method crash
//! -> this should rather be used for debugging
//!    since ANY particle should fulfill the CC
//! Use 0 to turn off
constexpr int DELETE_PARTICLE_IF_V_EXCEEDS = 20;


//! define if particle shall be sorted by weight into cell in order to
//! accelerate merging procedure. Since code does barely decelerate
//! with sorting, this should be always true.
constexpr bool SORT_BY_WEIGHT = true;

//! mark and trace particles
constexpr bool TRACK_PARTICLE = true;

//! ----------------------------------------------------------------------


//! Includes are here above variables can be used as needed in defines.h
#include "defines.h"
#include "CBlk.h"
#include "types.h"
#include "constexpr_sqrt.h"

#include <math.h>

//! below has to be included on SuSE 11
#include <string.h>
#include <stdlib.h>


//! GENERAL REMARKS:
//! - Since virtually all parameters are given in normalized units corresponding
//!   to the common Hybrid Code Normalization (See eg. Thorsten Bagdonats (TB) PhD Thesis pp.9),
//!   it is not mentioned additionally at each variables describtion.
//! - In case parameters use different units, it is explained in the respective comment
//! - There are two further files, which can be modified by the user:
//!   "CHybrid_FieldNames.cpp"  to name Fields as they shall appear in visualisation
//!   "CHybrid_IonProfiles.cpp" to define obstacle ion profiles


//! SOME GENERAL TERMS:
//! - Node:                 Point on the numerical Mesh on which physical fields are defined.
//!                         (abbriviated by Nd).

//! - Ghost Node:    Node at the boundary of the numarical mesh (abbriviated by GN).
//!                         The values on GN are not calculated, but copied or fixed.

//! - Resolution:    Nodes per Length.

//! - Level:        Always meaning the level of refinement. Incrementing
//!                         the Level by one results in doubling the resouloution.

//! - Cell:                 The cuboid Volume surrounded by 8 Nodes

//! - Zone:                 Identical to cell.

//! - Block:        Numerical Mesh, consisting of a defined number of Nodes
//!                         for each dimension. Each Block is surrounded by a layer
//!                         of Ghost Nodes.
//!                         (Block is abbriviated by Blk, its Nodes by BlkNds).

//! - Root Block:     Block at Level 0 (abbriviated by RB).
//!                          Number of Root Blocks remains constant within the entire
//!                          simulation run. Root Blocks CANNOT be added or removed after
//!                          run has started.

//! - Oct:                  One eighth of a block the is build by halfing the block
//!                         length in each direction

//! - Simulation-Box: Collectivity of all Blocks, defining the entire discretised scenario.
//!                          The initial Simulation Box is exclusively build up by Root Blocks.
//!                          Even though the RootBlocks are refined during the simulation run,
//!                          the length (size/Volume) of the Simulation Box does NOT change, only
//!                          resolution is spacially increased.
//!                          Initially the number of Nodes in one dimension (eg. X) of the
//!                          Simulation Box is given by by RB_X *(BlkNds_X-2).

//! - Macro Particle: Physical Ions are represented by Macro Particles (MP). Each MP
//!                         represents a certain number of Ions depending on the MP's weight
//!                         Do not confuse Weight and (Ion)Mass!

//! - Cross Section:  Cut through simulation box data of (CS)

//! - Child Block:    Each Block can be refined in ...


//! Parameter are sorted in terms of
//! 0) Normalization and Background Parameter
//! 1) Geometry Parameter
//! 2) Output Parameter
//! 3) Numerical Parameter
//! 4) Plasma Parameter
//! 5) Solar Wind Background Parameter (removed in current version)
//! 6) Obstacle Parameter



//!-------------------------------------------------------------//
//!-- 0) Normalization and Solar Wind (sw) Background Parameter //
//!-------------------------------------------------------------//


//! Natural and fundamental constants
constexpr D_REAL pi = M_PI;
constexpr D_REAL mu_0 = 4 * pi * 1.e-7;
constexpr D_REAL m_p = 1.672631e-27;
constexpr D_REAL e = 1.602e-19;
constexpr D_REAL kB = 1.38065e-23;
constexpr D_REAL amu = 1.660538782e-27;
constexpr D_REAL c_SI = 299792458;
constexpr D_REAL eps_0 = 8.854187871e-12;

//! parameters of plasma

//! Normalisation quantities
//! Magnetic Field in Telsa
constexpr D_REAL fieldmag = ce_sqrt(15*15 + 24*24 + 75*75);
constexpr D_REAL SI_B0 = fieldmag*1.e-9;
//! Number Density in Particle per Cubicmeter
constexpr D_REAL SI_n0 = 12.e+6;
//! Mass of normalization species in amu
constexpr D_REAL SI_m0 = 14.*m_p;
//! Alfvén Velocity in Meters per Seconds
constexpr D_REAL SI_v0 = SI_B0/ce_sqrt(mu_0*SI_n0*SI_m0);
//! Ion Inertia Length in Meters
constexpr D_REAL SI_x0 = SI_v0*SI_m0/(e*SI_B0);
//! Inverse Ion Gyration Frequency in Seconds
constexpr D_REAL SI_t0 = SI_m0/(e*SI_B0);

//! IMF (Magnetic Background BField) in normalized units
constexpr D_REAL B_sw[3] = {-15.0/fieldmag, 24.0/fieldmag, -75.0/fieldmag};

//! Angle of magnetic field to solar wind velocity in deg
constexpr D_REAL B_angle = 0;

//! Velocity of Solar Wind / upstream plasma
//! in Metre Per Second;
constexpr D_REAL SW_v = 130.e+3;

//! Alfvénic Mach Number
constexpr D_REAL MA = SW_v/SI_v0;
//! upstream plasma velocity in normalized units
constexpr D_REAL V_sw[3] = {MA, 0., 0.};
// constexpr D_REAL V_sw[3] = {8.2, 0., 0.};

//! Magnetic moment normalization
constexpr D_REAL SI_M0 = 4*pi*SI_B0*SI_x0*SI_x0*SI_x0/mu_0;

//! Calculation for Plasma Betas
//! Beta = thermalPressure / magnetic Pressure
//! Beta = calcBeta*Temperatur (in Kelvin)
//! calcBeta = SI_n0 * kB /(SI_B0*SI_B0/(2*mu_0))
constexpr D_REAL calcBeta = SI_n0*kB*2*mu_0/(SI_B0*SI_B0);
//! normalized electron mass for electron pressure calculation
constexpr D_REAL mass_electron_norm = 9.10938215e-31/SI_m0;

//! ---------------------------------------------------------------------------
//! hendrik's Enceladus parameters
//! ---------------------------------------------------------------------------
//! use realistic Enceladus plume model or simple test geometry
constexpr bool use_test_scenario=false;
//! which Enceladus flyby
constexpr short flyby = 13; 
//! electron impact ionization frequency in 1/s
constexpr D_REAL Ence_nu_e_total = 1.5e-09; 
//! upstream ion density of all flybys (cf. table 2.2 of phd thesis H.Kriegel)
constexpr D_REAL var_ion_dens[] = {90e+06,60e+06,65e+06,90e+06,55e+06,90e+06,45e+06,60e+06,45e+06,
					75.e+06,65e+06,60e+06,55e+06,65e+06,85e+06,90e+06,105e+06,110e+06,100e+06,70e+06};
  //! SI normalization values
constexpr D_REAL frac_of_species[5]={0.225 + 0.775*use_test_scenario,0.225,0.225,0.225,0.1}; //! O, OH, H2O, H3O, H  NOTE: sum has to be equal to one
constexpr D_REAL average_ion_mass = (frac_of_species[0]*16.+frac_of_species[1]*17.+frac_of_species[2]*18.+frac_of_species[3]*19.+frac_of_species[4]*1)
				*!use_test_scenario + 18.*use_test_scenario;
//! Parameter analytical Plume (if used)
constexpr D_REAL base_density = 1e+15;	//! neutral peak density in SI (m^3)
constexpr D_REAL xi_nD_0 = 1.e+09; 		//! dust peak density in SI (m^3)
constexpr D_REAL H_d = 0.25*252.e+03/SI_x0;
constexpr D_REAL H_d_dust = 948.e+03/SI_x0;
constexpr D_REAL opening_angle_gas = 10;	//! in degrees
constexpr D_REAL opening_angle_dust = 20.;
constexpr D_REAL theta0=0;			//! tilt of plume
constexpr D_REAL phi0= 0;
//! ---------------------------------------------------------------------------


//! ---------------------------------------------------------------------------
//! christoph's comet parameters
//! ---------------------------------------------------------------------------
//! outgassing velocity in Meters per seconds
constexpr D_REAL COM_v0 = 1.e+3;
//! Gas Production Rate of the Comet in 1 per Second
constexpr D_REAL COM_Q[] = {5.e27,0.,0.};
constexpr D_REAL COM_Q_Total = 5.e27;
//! Coll. freq. neutrals/ions in m^3/s
// constexpr D_REAL SI_nu = 1.1e-15;
//! Radius of Cometary Ion Injection
constexpr D_REAL COM_a_SI = 2000.;
//! Raw Density for Neutral Density
//! n_N = COM_RawDensity*1/r^2
constexpr D_REAL COM_RawDensity = 1/(4*M_PI*COM_v0*SI_x0*SI_x0*SI_n0);
//! Calc Com_a
constexpr D_REAL COM_a = COM_a_SI/SI_x0;
//! ---------------------------------------------------------------------------


//!-------------------------------------------------------------//
//!------------- 1) Geometry Parameter -------------------------//
//!-------------------------------------------------------------//

//! Number of Nodes in each Block:
//! - only even values allowed (otw. false interpolation)
//! - Ghost Nodes are included
//!   -> only BlkNds_-2 Nodes are used for physical Field !
//!   NOTE:
//! ------------------------------------------------------------------
//! - USE EQUAL NUMBER IN EACH DIRECITON IN CASE REFINEMENT IS USED, -
//! - ELSE INTERPOLATION AT BOUNDARIES FAILS FOR SOME REASON !!!     -
//!   (maybe child to parent field method ?!?)
//! ------------------------------------------------------------------
constexpr INT32 BlkNds_X = 16;
constexpr INT32 BlkNds_Y = 16;
constexpr INT32 BlkNds_Z = 22;

//! Number of Root Blocks in Szenario
constexpr INT32 RB_X = 16;
constexpr INT32 RB_Y = 16;
constexpr INT32 RB_Z = 16;

//! Number of maximal refinement-levels used
constexpr INT32 MAX_LEVEL = 2;

//! Radius of moon
//! (is used instead of R_Obstacle for geometry of neutral cloud since some
//!  test simulations had R_Obstable=0 but boxsize and axis should still
//!  be given in Enceladus radii)
constexpr D_REAL R_Moon = 2634100./SI_x0;

//! Size of simulation box
constexpr D_REAL LX = 30*R_Moon;
constexpr D_REAL LY = 30*R_Moon;
constexpr D_REAL LZ = 40*R_Moon;
//! Origin of simulationbox (= Position of Obstacle)
constexpr D_REAL Box_Origin[3] = {0.4001 * LX, 0.5001 * LY, 0.5001 * LZ};

//! Radius of Obstacle
constexpr D_REAL R_Obstacle = R_Moon;
constexpr D_REAL squared_R_Obstacle = R_Obstacle*R_Obstacle;

//! use periodic Boundaries;
//! every particle leaving the simulation box at any side is inserted
//! at the opposite side of the simulation box. Field derivatives at
//! boundaties are calculated using values from the oposite side of the
//! simulation box.
constexpr bool use_periodic_bounds[3] = {false, false, false};



//! specify number of cells respective boundary set are
//! set to eta value defined below
constexpr INT32 resistive_bound_cells[]  = {0,0,
                                        0,0,
                                        0,0};

//! specify value of eta at boundary
constexpr F_REAL resistive_bound_eta[]  = { 0., 0.,
                                        0., 0.,
                                        0., 0.};

//! specify distance of respective boundary for eta to go
//! linearly from eta value above to zero
constexpr F_REAL resistive_bound_dist[]  = {0.,0.,
                                        0.,0.,
                                        0.,0.}; 

//! alfven wings are strongly reflected at the boundaries
//! therefore the wings have to be damped by using eta
//! over long distances (even longer than one block)
//! also optimal_MPiC is reduced in that block(s)
//! NOTE: works so far only for z-direction
constexpr bool eta_Alfven_Wing_boundary = false;

//! parameters used in calculating eta for the Alven wing boundaries
//! see implementation in CBlk_EtaProfile.cpp
//! Only used if eta_Alfven_Wing_boundary is true above
constexpr D_REAL eta_Alfven_xhalfaxis = 6*R_Moon;
constexpr D_REAL eta_Alfven_zhalfaxis = 7*R_Moon;
constexpr D_REAL superellipseexponent = 6;


//! General remark on Boundary conditions:
//! 1. entry: x_min-side        2. entry: x_max-side
//! 3. entry: y_min-side        4. entry: y_max-side
//! 5. entry: z_min-side        6. entry: z_max-side
//! in case "use_periodic_bounds" is set true, this parameter
//! will be ignored.


//!---------------------------------------------------------------------
//!-------------- Particle Boundary Conditions -------------------------
//!---------------------------------------------------------------------
//! 1) First decide, whether particle boundaries have to homogenous
//!      or not
//! 2) In case
//!      a) true  -> further specify whether in- or outflow shall be used
//!      b) false -> inflow & outflow specification will be ignored and
//!         inhomogenous boundaries will be used that are specified
//!         in "utils_boundaries.cpp"
//!NOTE: 7th entry is used for simulation initialization
constexpr bool use_hom_particle_bounds[7] = {
    true,  true,
    true,  true,
    true,  true,
    true
  };



//! 3) In case homogenous boundaries are used, set use_particle_inflow_bounds to
//!    a) true  ->  inject particle at respective side (inflow)
//!    b) false -> do not inject particle at respective side (outflow)
//! NOTE: Up to now applying outflow boundaries for particle at min-sides
//!       does not work well (Due to assymetry in gather method).
constexpr bool use_particle_inflow_bounds[6] = {
    true , false,
    false , false,
    true , true
  };


//! 3) In case outflow boundaries are used, set use_particle_copy_bounds to
//!    a) true  ->  copy particles from penultimate cell to ultimate
//!    b) false -> to eareas all particle in boundary cell
//!        In case respective value in "use_particle_inflow_bounds" is
//!        set true as well, values set here will be ignored.
constexpr bool use_particle_copy_bounds[6] = {
    false, false,
    false, false,
    false, false
  };


//!---------------------------------------------------------------------
//!----------------- Field Boundary Conditions -------------------------
//!---------------------------------------------------------------------



//! 1) First decide, whether BField boundaries have to homogenous
//!      or not
//! 2) In case
//!      a) true  -> further specify whether in or outflow shall be used
//!      b) false -> inflow & outflow specification will be ignored and
//!         inhomogenous boundaries will be used that are specified
//!         in "utils_boundaries.cpp"
//!NOTE: 7th entry is used for simulation initialization
constexpr bool use_hom_B_bounds[7] = {
    true, true,
    true, true,
    true, true,
    true
  };
  

//! set new IMF after TL_new_Bsw (e.g. sector border, magnetopause crossing)
constexpr D_REAL new_B_sw[3] = {0., 0., 1.};

//! TL when new BField boundaries shall be set
//! (set zero to switch off)
constexpr INT32 TL_new_Bsw = 0;

//! in order to avois sharp transition rotate old to
//! new field within finite time
constexpr INT32 num_TL_adapt =  -1;
  
  

//! - in case inhomogenous BField boundaries are used, they can be set
//!   depending on whether plasma is streaming into or out of the boundary
//!   at the respective cell
//! -> it is tested wether the boundary-normal is parallel to the velocity
//!    that is specified in "utils_boundaries.cpp"
//! NOTE:
//! This only makes sense in case inhomogenous particle conditions are
//! specified at the respective boundary
//! NOTE:
//! works for X-MIN/MAX boundaries only
constexpr bool set_BFieldInflow_at_VelocityInflow = false;



//! In case homogenous boundaries are used:
//! For each side of the simulation box has to be specified, which
//! type of boundary conditions shall be used for the EM-Fields:
//! true:  set fields constant to background value (inflow)
//! false: set fields derivative to zero (outflow)
//! in case "use_periodic_bounds" is set true, this parameter
//! will be ignored:

//! specify boundaries for Bfield
constexpr bool use_B_inflow_bounds[6] = {
    true, false,
    true, true,
    true, true
  };

//! specify boundaries for EField
  constexpr bool use_E_inflow_bounds[6] = {
    true, false,
    true, true,
    true, true
  };
  
//! specify boundaries for UField
  constexpr bool use_U_inflow_bounds[6] = {
    true, false,
    true, true,
    true, true
  };
  
//! Relevant if nonadiabatic_gradPE_TERM is true
  constexpr bool use_PE_inflow_bounds[6] = {
    true, false,
    true, true,
    true, true
  };
//!
				      
//! In case inhomogenous boundaries are used they must be specied in:
//! "utils_boundaries.cpp"
//! Either
//! 1) An anylytic function can be specified.
//! 2) A one dimensional file can be applied.

//! In case 2) following values have to be specified
//! It is implicitely assumed that extern boundary conditions are 1 dimensional
//! which means that they are either a function of X or Y or Z, but not a
//! function of eg. X and Y at the same time.

//! ----- EBF = Extern Boundary Field -----------------------------------

//! specify how many files shall be read
constexpr INT32 num_EBF = 0;

//! specify how many coordinate values have to be read for each data point
//! (this indicates whether the dataset boundary values are a 1D, 2D or 3D function)
constexpr INT32 num_coords_EBF[] = {3,3,3};

//! specify how many components the fields inside the files have
constexpr INT32 num_comps_EBF[] = {1,1,1};

//! specify name of the respective file
constexpr char file_names_EBF[][RUN_NAME_SIZE] = {"nSW_bounds.txt" , "nCI_bounds", "u_Bounds.txt"};



//!-------------------------------------------------------------//
//!------------- 2) Output Parameter ---------------------------//
//!-------------------------------------------------------------//

//! name of simulation run (used as a prefix for each file)
constexpr char Run_Name[RUN_NAME_SIZE] = "Ganymede";

//! set output path relative to directory of executable
constexpr char data_output_path[RUN_NAME_SIZE] = "../data";


//! End Run when TL_MAX is reached
constexpr INT32 TL_MAX = 60000;

//! Perform in-situ timing for logging and performance evaluation
//! NOTE: This has significant performance impact, it's better to use
//!		  an external tool such as vtune
constexpr bool do_timing = false;

//! SET 0 TO SWITCH OFF RESPECTIVE OUTPUT 

//! vtk    is fileformat for Paraview
//! silo   is fileformat for VisIt
//! native is fileformat for own Visualization

//! zonal Data: Each Zone/Cell is rendered using the color value of the
//!                 lower left Node, resulting in a mosaic-like plot.
//!                 Even though VisIt offers th possibility do visualize
//!                 nodal data in a zonal fashion, some strange interpolatoion is used.
//!                 Hence it is impossible to identify which value is defined on a
//!                 certain node (-> not usefull for e.g. debugging).
//! nodal Data: Each Node is asociated with one color value,
//!                 color of areas inbetween nodes are interpolated.
//! set true for zonal style.
//! set false for nodal style.
constexpr bool SILO_STYLE_ZONAL = false;


constexpr INT32  TL_OUTPUT_2D_SILO    =    500;
constexpr INT32  TL_OUTPUT_2D_NATIVE  =      0;
constexpr INT32  TL_OUTPUT_2D_GNUPLOT =      0;

constexpr INT32  TL_OUTPUT_3D_SILO	 =     0;
constexpr INT32  TL_OUTPUT_3D_NATIVE  =    0;
constexpr INT32  TL_OUTPUT_3D_uniform_grid =   10000;
constexpr INT32  TL_OUTPUT_3D_ASCII   =    0;

//! Enable HDF5 GZIP Compresion for SILO-Output
//! Approx. 30% smaller output files but needs slightly more computation time

constexpr bool COMPRESS_SILO = true;

//! write file with total enegy, momentum and mass of ions.
constexpr INT32  TL_OUTPUT_ErgMomMass =    0;

//! write Line Output

constexpr INT32  TL_OUTPUT_LINEOUT    =   0;
//!  Particle Detctor Output
constexpr INT32 TL_OUTPUT_PARTICLE_DETECTOR =  0;
//! track path of trajectory and write 1dim field data
constexpr INT32  TL_OUTPUT_TRAJECTORY = 500;

constexpr INT32  TL_OUTPUT_PARTICLE_TRACKS = 0;

constexpr INT32  TL_OUTPUT_GETMESH = 0;

//!-----------------------------------------------------------------
//!---------- statefile related -------------------------------
//!-----------------------------------------------------------------

//! save state every xx time level (set to 0 to switch off)
//! NOTE:
//! For each block the responsible MPI process will be stored, hence simulations that
//! were operated on N process must not be re-started with a number of processes
//! smaller than N, yet they can be re-started with a number larger than N.
constexpr INT32 TL_SAVE_STATE  = 2000;

//! OPTOIONS:
//! - LEAVE_UNCHANGED
//! - CONVERT_TO_TRACK
//! - CONVERT_TO_ORDINARY
constexpr INT32 convert_TRACK_PARTICLE_State = LEAVE_UNCHANGED;

//! in general all processes should write their state file at the same time, however
//! some filesystem cannot handle that much data at once. In such case this parameter 
//! should be set true in order to force the processes to write their state file 
//! one after the other.
constexpr bool serialize_writing_of_StateFile = true;

//! Secure State File
//! after a successfull restore the code mv the folder State to StateTLXXXXXX
constexpr bool secure_state_file = false;

//! resubmit of jobfile
//!resubing job if TL_Max is not reached at the end of the Walltime
//!Walltime must be in second
constexpr bool resubmit_enable = false;
constexpr INT32 resubmit_walltime =400*3600;
constexpr char resubmit_jobscript[100] = "msub job.sh";
constexpr INT32 resubmit_security_factor = 60*20;

//!-----------------------------------------------------------------
//!---------- silo output related -------------------------------
//!-----------------------------------------------------------------

//! use factor eg. to scale coordinated from normalized to SI units
//! in visit visualisation.
constexpr F_REAL factor_scale_mesh = 1./R_Moon;

//! Only for 2D output
//! Do not cut through GN Planes in depth direction
//! GN surroundig the Blks are always printed !!!
//! in case set to false all GN have to be
//! update before output.
//! In case set to true update GN are not required,
//! however the cross section might be slightly shifted
//! versus th correct position.
constexpr INT32  avoid_GhostCXS   =  true;

//! Cross Section (CS) to cut (in normed units with respect to origin=(0,0,0))
constexpr F_REAL CROSS_SECTION[3] = {0., 0., 0.};

//! specify label name for VisIt visualization
constexpr char visit_xlabel[20] = " x";
constexpr char visit_ylabel[20] = " y";
constexpr char visit_zlabel[20] = " z";

//! specify length unit of coordinate axis
constexpr char visit_length_unit[20] = "R_G";


//!-----------------------------------------------------------------
//!---------- trajectory output related -------------------------------
//!-----------------------------------------------------------------

//! Name of trajectory file
//! - Specify trajectory in normalized units relative to Box-Origin
//! - the following example plots the sub-solar line
//!   ending at the box origin
//! - fieldxyz means the data measured by a spacecraft 
//! (time_string)               -30.    0.      0.      (fieldx) (fieldy) (fieldz)
//! (time_string)               -20.    0.      0.      (fieldx) (fieldy) (fieldz)
//! (time_string)               -10.    0.      0.      (fieldx) (fieldy) (fieldz)
//! (time_string)                 0.    0.      0.      (fieldx) (fieldy) (fieldz)


//! specify whether additional time string should be read
constexpr bool do_read_timestring = true;

//! reading the measured field (column 5,6,7) and write into new file
constexpr bool do_read_SpaceCraftField = true;

//! specify how many trajectories shall be traced
constexpr INT32 num_trajectories = 2;

//! specify names of files in which positions of trajectories are strored
constexpr char Trajectory_FileName[][RUN_NAME_SIZE] = 
{"juno_x","juno_kernel","juno_kernel_noZ","juno_mag","juno_mag_noZ","trajx_-5.0","trajz_-1.1","trajz_-1.5","trajz_-2.0","trajz_-2.5","trajz_-4.0","trajz_-5.0",
"trajz_1.1","trajz_1.5","trajz_2.0","trajz_2.5","trajz_4.0","trajz_5.0"};


//! if trajctory file is not given in normalized units, the positions
//! have to be normalized (by ion inertia length)
constexpr D_REAL trajectory_pos_conversion_fact=1./R_Moon;

//! decide whether to trace trajectory or only to display in 2D out
//! (latter one applyies especially to marks such as MP or BS)
constexpr bool trajectory_2D_out_only[] = {false, false, false, false, false};

//! compose single file from all files of parallel trajectory output
//! -> will make the code somewhat slower but result in significantly
//!    less files
constexpr bool pack_parallel_output_to_single_file = true;

//! How many fields shall be traced each trajectory
constexpr INT32 num_fields_to_trace = 10;


//! NOTE:
//! In order to trace average fields, use
//! for avered field 1: id_average_Field1 +0
//! for avered field 1: id_average_Field1 +1
//! for avered field 1: id_average_Field1 +2
//! etc ...

//! specify which fields to trace via field id defined above
constexpr INT32 ID_of_Fields_to_trace[] = {id_BTotal, id_rhoSpecies1, id_rhoSpecies1+1, id_rhoSpecies1+2,
				       id_rhoSpecies1+3, id_rhoSpecies1+4, id_rhoSpecies1+5, id_rhoSpecies1+6, 
				       id_rhoSpecies1+7, id_rhoSpecies1+8};

//! decide whether new file shall be created each time level
//! or output shall be appended to existing file
//! (e.g. for time series)
constexpr bool create_new_file_each_TL = true;

//!-----------------------------------------------------------------
//!---------- Lineout Output related ----------------------------
//!-----------------------------------------------------------------


//! Number of lines for output
constexpr INT32 num_lines = 0;

//! Conversation Factor to normalised Units
//! if lines is not given in normalized units, the positions
//! have to be normalized (by ion inertia length)
constexpr D_REAL lineout_pos_conversion_fact=1.;

//! Names of Lines
//! these will be the names of the output files
//! since these names are also used as variable names for
//! silo ouput, only alphanumeric characters are allowed
//! (especially no "." )
constexpr char LineOut_FileName[][RUN_NAME_SIZE] = {
    "lineX",
    "lineY",
    "lineZ"
  };


//! Number of fields in Lineout Output                                              
constexpr INT32 num_fields_to_trace_lineout = 1;

//! specify which fields to trace via field id defined above
constexpr INT32 ID_of_Fields_to_trace_lineout[] = {id_numberdensity_neutralSpecies1,   id_BTotal,ENERGY_SPECTRUM,
    id_PEtotal
    , id_rho_np1
    , id_UI_plus
    , id_divU
    , id_rotB
    , id_average_Field1 + 0 /* id_PEtotal */
    , id_average_Field1 + 1 /* id_rho_np1 */
    , id_average_Field1 + 2 /* id_UI_plus+0 */
  };

//! compose single file from all files of parallel trajectory output
//! -> will make the code somewhat slower but result in significantly
//!    less files
constexpr bool pack_parallel_output_to_single_file_lineout = true;

//! decide whether new file shall be created each time level
//! or output shall be appended to existing file
//! (e.g. for time series)
constexpr bool create_new_file_each_TL_lineout = true;


//!-----------------------------------------------------------------
//!---------- Average fields related -------------------------------
//!-----------------------------------------------------------------

//! NOTE --- HOW TO TRACE AVERAGE FIELDS ---
//! When average fields shall be traced specify:
//! id_average_Field1 +0,
//! id_average_Field1 +1
//! id_average_Field1 +2
//! etc.
//! where id_average_Field1 is the first field specified
//! below in "IDs_of_Fields_to_average"


//! specify how many time levels averaging shall
//! start before "TL_FINISH_AVERAGE_FIELDS"
constexpr INT32  num_TL_average_Fields = TL_OUTPUT_2D_SILO-1;

//! specify prefix that is added to field name in order to
//! distinguish original vs avereged field
constexpr char average_field_name_prefix[] = "avrg";

//! specify when averaging shall be finished
//! (most likely equal to output)
//! set zero to turn off averaging
constexpr INT32  TL_FINISH_AVERAGE_FIELDS = TL_OUTPUT_2D_SILO;


//! eg. if TL_FINISH_AVERAGE_FIELDS = 20 and
//!        num_TL_average_Fields = 10
//!        fields of time levels [10:19], [30:39], ...
//! will be averaged.
//! For this example average-fields will be normalized
//! (ie. devided bynum_TL_average_Fields) at TL=19,39, ...


//! specify fields that shall be averaged via their ID
//! as specified above
constexpr INT32 IDs_of_Fields_to_average[] = {
    id_PEtotal
    , id_rho_np1
    , id_UI_plus
  };


//!-----------------------------------------------------------------
//!---------- uniform_grid Output related -------------------------------
//!-----------------------------------------------------------------

//! some specifications for uniform_grid output
constexpr F_REAL uniform_grid_B0SI = SI_B0;
constexpr F_REAL uniform_grid_n0SI = SI_n0;
constexpr F_REAL uniform_grid_x0SI = SI_x0;
constexpr F_REAL uniform_grid_v0SI = SI_v0;

//! Since the uniform grid cannot handle block adaptive meshes,
//! values of the simulation code's adaptive mesh are interpolated
//! to a catesian mesh for the uniform_grid Tool. The number of mesh nodes
//! of the cartesian mesh has to be specified in "numNds_uniform_gridMesh".
//! The resolution should be choosen in such way that it matches up
//! with the resolution of the adaptive mesh at the highest level
//! of refinement (RB_X*(BlkNds_X-2)*MAX_LEVEL), but one has to be 
//! careful about the size of the resulting file!

//! If the uniform_grid Output should be just a part of the simulation boxsize
//! (e.g. for smaller filesize), output box length could be scaled
constexpr F_REAL uniform_grid_shrink_domain[3] = {3., 3., 4.};

constexpr INT32 numNds_uniform_gridMesh[3] = { 	(int)(RB_X*(BlkNds_X-2)*(MAX_LEVEL+1)/uniform_grid_shrink_domain[0])	,  
					   (int)(RB_Y*(BlkNds_Y-2)*(MAX_LEVEL+1)/uniform_grid_shrink_domain[1])	,
					   (int)(RB_Z*(BlkNds_Z-2)*(MAX_LEVEL+1)/uniform_grid_shrink_domain[2]) };
					
//! number of fields which should appear in uniform_grid visualization
constexpr INT32 uniform_grid_num_fields_to_trace = 10;

//! id's of fields which should appear in uniform_grid visualization
constexpr INT32 uniform_grid_ID_of_Fields_to_trace[] = { id_BTotal, id_EField, id_rhoSpecies1, id_rhoSpecies1+1,
						     id_rhoSpecies1+2, id_rhoSpecies1+3, id_rhoSpecies1+4,
						     id_rhoSpecies1+5, id_rhoSpecies1+6, id_rhoSpecies1+7};

//! type of field for id's selected above, neccessary to assign correct
//! units to fields in visualization
//! possible values are uniform_grid_TYPE_BFIELD for nT, uniform_grid_TYPE_EFIELD
//! for V/km,
//! uniform_grid_TYPE_VELOCITY_FIELD for km/s and uniform_grid_TYPE_DENSITY_FIELD
//! for 1/cm3
constexpr INT32 uniform_grid_Field_type[] = {uniform_grid_TYPE_BFIELD, uniform_grid_TYPE_EFIELD,
					 uniform_grid_TYPE_DENSITY_FIELD, uniform_grid_TYPE_DENSITY_FIELD,
					 uniform_grid_TYPE_DENSITY_FIELD, uniform_grid_TYPE_DENSITY_FIELD,
					 uniform_grid_TYPE_VELOCITY_FIELD, uniform_grid_TYPE_VELOCITY_FIELD,
					 uniform_grid_TYPE_VELOCITY_FIELD, uniform_grid_TYPE_VELOCITY_FIELD};


//!--------------------------------------------------------------
//!------------ Track Particle (Test Particle Simulation) -------
//!--------------------------------------------------------------

//!--------------------------------------------
//! NOTE:                                     -
//! TRACK_PARTICLE has to be set above		  -
//!--------------------------------------------

//! In test particle simulations any function is deactivated except
//! for particle movement and acceleration
constexpr bool TestParticle_Simulation = false;


//! Trace Particles:
constexpr INT32  num_mark_particle_in_species[] = {0, 0, 0, 0,0,0,0,0};
//! NOTE: this number is only exact, if INJECT_PARTICLE_TO_TRACK
//!	  is set to true, else it is rather a proxy

//! volume to mark particles (xmin,xmax,ymin,ymax,zmin,zmax)
constexpr D_REAL mark_particle_volume[6] ={-5.*R_Moon,-3.*R_Moon ,
					-1.*R_Moon, 1.*R_Moon,
					-3.*R_Moon,-1.*R_Moon };
					
//! Trace Particles:
// constexpr INT32  num_mark_particle_in_species[] = {100, 0, 0, 0};

//! if particles shall be tracked in full selfconsistent hybrid
//! simulations, some particles can be "marked"
//! -> choose in function "mark_particle" which to mark
//! set TL_MARK_PARTICLE_TRACKS to -1 in case
//! particles shall not be marked
constexpr INT32  TL_MARK_PARTICLE_TRACKS   = -1;

//! in case a Test Particle Simulation is carried ou in which
//! particle are not initially available, inject particle via
//! function "inject_particle_to_track" right after fields are
//! restored
constexpr bool INJECT_PARTICLE_TO_TRACK = false; 

//! in case a Test Particle Simulation is carried ou in which
//! particle are not initially available, inject particle via
//! function "inject_particle_to_track" right after fields are
//! restored
constexpr INT32 TL_INJECT_PARTICLE_TO_TRACK = -1;

constexpr INT32 TL_CONVERT_PARTICLE_TRACKS_TO_SILO_MESH = 0;

constexpr bool binary_particle_tracks = false;

constexpr INT32 num_tracks_each_group = 100;

//!--------------------------------------------------------------

//!-----------------------------------------------------------------
//!---------- Particle Detector output related ---------------------
//!-----------------------------------------------------------------

//! Number of detector boxes
constexpr INT32 num_particle_detector = 0;
constexpr char Detector_FileName[][RUN_NAME_SIZE] = { "detector1",
                                                    "detector2",
                                                    "detector3"};

//! Detector Box Geometry

constexpr D_REAL detector_box_xmin[] = {-100000./SI_x0,0,0};
constexpr D_REAL detector_box_xmax[] = {-60000./SI_x0,0,0};
constexpr D_REAL detector_box_ymin[] = {-10000./SI_x0,0,0};
constexpr D_REAL detector_box_ymax[] = {10000./SI_x0,0,0};
constexpr D_REAL detector_box_zmin[] = {-10000./SI_x0,0,0};
constexpr D_REAL detector_box_zmax[] = {10000./SI_x0,0,0};

//!-------------------------------------------------------------//
//!-------------- 3a) Numerical Parameter ----------------------//
//!-------------------------------------------------------------//


//! Numerical time step at Level 0.
//! NOTE:
//! dt should be at least 5 times smaller than
//! Courant Criteria suggest. Otherwise density 
//! lacks at Blk Levl borders occur in flow direction.
//! TODO MAX_LEVEL is set BELOW this parameter!
// constexpr D_REAL  dt =  0.1* RB_X*(BlkNds_X-2)/LX / (SW_v/SI_v0);
constexpr D_REAL  dt =  0.009; //0.012 //! set at 0.0004 for symmetric 40nT run

//! specify mesh type
//! 1) UNIFORM
//! 2) STAGGERED
constexpr INT32 mesh_type = UNIFORM;


//! one step:
//! 1) dtB = rot(uxB-...)
//! two step
//! 1)   E = -uxB+...
//! 2) dtB = -rotE
constexpr bool LF_one_step = true;

//! Choose advance_B_Algorithm for Plasma Equations:
//! 0: Leap Frog   (fast)
//! 1: Runge Kutta (a bit more stable than LF)
constexpr INT32 advance_B_Algorithm = 1;

//! advance_B_loop includes:
//! - advance B Plasma 
//!   -> num_sub_cycle times, see above)
//! - advance B Obstacle
//! - divergence cleaning
constexpr INT32 num_advance_B_loops = 3;


//! plasma resistivity
//! could be used instead of smoothing E and B
//! similar numerical diffusion is achieved if
//! smooth_EB = Eta_sw * dt / (cell size)^2 (all in normalized units)
//! (see phd thesis kriegel, pages 71 - 73)
//! rule of thumb: Eta_sw = 0.1 * cell size 
constexpr D_REAL Eta_sw = 0.01;

//! strength of EM-Field smoothing
//! (set value for each refinement level)
constexpr bool several_smooth_in_levels = false; //! if true, smoothing is carried out 2^Level times
constexpr D_REAL smooth_Factor = 0.016;
constexpr D_REAL smooth_E[] = { smooth_Factor * 1., smooth_Factor * 2., smooth_Factor *3. , smooth_Factor * 8.};
constexpr D_REAL smooth_B[] = { smooth_Factor * 1., smooth_Factor * 2., smooth_Factor *3. , smooth_Factor * 8.};

//! number & strength of smoothing the Resistivity profile of the obstacle
//! values are used only if use_resistive_obstacle (section 6) is set to true 
constexpr INT32  num_smooth_eta = 5;
constexpr D_REAL smooth_eta[] = {2./8. *1.0,  1./4. *0.5, 1./4. *1.0, 1. * 0.5,  1.* 0.5};

//! Relevant if nonadiabatic_gradPE_TERM_smooth_PE is true
constexpr D_REAL smooth_PE[] = { smooth_Factor * 1., smooth_Factor * 2., smooth_Factor *4. , smooth_Factor * 8.};
//!



//! Minimal Charge Density (MCD_J2U) is used in "convert_Ji2Ui" function
//! as an lower limit for rho. Applying MCD_BField results in strongly
//! underestimated velocities.
//! (rho meaning the total density of all ion species)
//! NOTE: this has to be larger than the smallest density different from zero.
//! NOTE: Has to be MCD_BField or comparable for negative species
//!       since a low density cell can have large currents 
constexpr D_REAL MCD_J2U = 1e-4;

//! Minimal Charge Density (MCD_BField) is used 
//! in "magnetice Field" function as an lower limit for rho.
//! (rho meaning the total density of all ion species)
constexpr D_REAL MCD_BField = 0.1;



//!-------------------------------------------------------------//
//!-------------- 3b) Refinement Parameter ---------------------//
//!-------------------------------------------------------------//


//! Apply particle splitting and merging every x TL
//! set 0 to switch off
//! NOTE:
//! if merged particle enter split area "stripes" occur !?!
//! (at least for cold plasma)
constexpr INT32 TL_SPLIT = 1; 
constexpr INT32 TL_MERGE = 1;

//! optimal_MPiC is the aspired number of Macro Particle for a certain species in each cell.
//! In case of inflow_species, each cell is initialized with this number of
//! particles for the respective species.
//! During simulation the number of particles in a cell will change, so splitting
//! and merging procedures try to adjust the active number to optimal_MPiC.
constexpr INT32 optimal_MPiC[] = {16, 16, 16, 16, 10, 10, 10, 5, 20, 10, 30, 30};

//! NOTE:
//! In case to many particle are inside one cell (>8000), the computational
//! speed decreases exponentially hence the code gets AWFULLY slow!
//! Additionally using several thousands particle inside a single cell does
//! not make sense anyway, so "num_merge_each_cell" must be increased.

//! define how many merge/split processes should be carried out in each cell
//! -> the higher oPiC the higer these parameters should be set
constexpr int num_merge_each_cell = 300;
constexpr int num_split_each_cell = 40;

//! ---- DEFINE HOW TO SPLIT/MERGE -------------------------------
//! Decide wheter to split heaviest or most centred particle
//! -> in case most centred is chosen, CoM is conserved by
//!    default an parameters below do not affect splitting
constexpr bool split_heaviest_particle = true;


//! In order to speedup the merging process (which is O(N²)), the
//! best triple of particle is estimated beyond "num_particle_in_MergeList"
//! particle. Since all particle are weight sorted it is ensured, that only
//! the lightest particle on one cell are inside the mergelist
constexpr INT32 num_particle_in_MergeList = 0.5*(optimal_MPiC[0]+optimal_MPiC[1]+optimal_MPiC[2]);//+optimal_MPiC[3]);

//! In order to avoid very high macroparticle numbers in the tail
//! if cell position > merge_tail_distance => merge
constexpr D_REAL merge_tail_distance = 0.; //previously set at 3

//! To slightly violate centre of mass conservation drastically
//! reduces noise and avoids numerical artifacts
//! -> in summary this seems to reflect the physics much better
//! -> THIS SHOULD ALWAYS BE SET TO FALSE EXCEPT FOR TESTING !!!!
constexpr bool conserve_CoM_at_merging = false;
constexpr bool conserve_CoM_at_splitting = false;



//! new randomly generated positions may be out of cell
//! -> define how often to retry generating random position
//! -> this parameter will be ignored in case 'conserve_CoM'
//!    is set false
constexpr INT32 num_randomize_position_merge = 200;
constexpr INT32 num_randomize_position_split = 100;

//! level individual activation of splitting/merging
//! in general it should not be required to merge
//! or split in the highest level when "fac_oMPiC_at_LevBoundBlks"
//! is set to 8 or higher
constexpr bool do_merge_in_L[] = {true, true, true, false};
constexpr bool do_split_in_L[] = {true, true, true, false};

//! decide whether to split exclusively particle of certain
//! species, eg. do not split planetary Ion species
constexpr bool do_split_in_species[] = {true, true, true,true,true,true,true,true,true};
//constexpr bool do_split_in_species[] = {true, false, false, false, false, false, false, false};//true,true,true,true,true,true,true};
constexpr bool do_merge_in_species[] = {true, true, true,true,true,true,true,true,true}; //true,true,true,true,true,true,true};
//! --------------------------------------------------------------





//! ---- GENERAL REFINEMENT OPTIONS -------------------------------


//! increase oMPiC at boundary blocks neighbours in order to
//! OBTAIN smoother transition
constexpr F_REAL fac_oMPiC_at_LevBoundBlks = 2.;


//! decide wheter to smooth or to simply inject moments to parent:
//! Assume L0 is refined to L1 hence particle are in L1.

//! Using the smooth version gives exactly the same density profile for L0,
//! as if particle would have been collected in L0.
//! This only cannot be fullfilled at level boundaries:
//! At minus boundaries for first node inside refined block (SURE ???)
//! At plus boundaries for first node outside refined block
//! (can be seen in visu, follow reason using sketch)
constexpr bool smooth_moments_to_parent = false;


//! This in-accuricy described above appears with and without using gather blocks,
//! so there seems to be no real advantage in using gather blocks.
constexpr bool use_gather_blocks = false;
//! ----------------------------------------------------------------





//! ---- INITIAL HIERARCHICAL REFINEMENT -------------------------------

//! A) CUBOID
//! Specify minimal and maximal Coordinate relative to Box Origin
constexpr INT32 TL_REFINE_STATIC_CUBOID[] = {0, 0, -1};

//! "left, lower corner"
constexpr D_REAL minX_refBox_of_L[] = {-6.0* R_Moon, -4.0* R_Moon, -0.2* R_Obstacle};
constexpr D_REAL minY_refBox_of_L[] = {-8.0* R_Moon, -4.0* R_Moon, -0.2* R_Obstacle};
constexpr D_REAL minZ_refBox_of_L[] = {-8.0* R_Moon, -4.0* R_Moon, -1.3* R_Obstacle};
//! "right, upper corner"
constexpr D_REAL maxX_refBox_of_L[] = { 10.* R_Moon,  6.0* R_Moon,  0.2* R_Obstacle};
constexpr D_REAL maxY_refBox_of_L[] = { 8.0* R_Moon,  4.0* R_Moon,  0.2* R_Obstacle};
constexpr D_REAL maxZ_refBox_of_L[] = { 8.0* R_Moon,  4.0* R_Moon,  -1.0* R_Obstacle};


//! B) SPHERE
//! Specify radius of sphere in which Blks shall be refined

constexpr INT32 TL_REFINE_STATIC_SPHERE[] = {-1, -1, -1,};

constexpr D_REAL radius_refSphere_of_L[] = { R_Moon+(1.5*1560800./SI_x0), R_Moon+(1560800./SI_x0), R_Moon+(3121600./SI_x0), 0., 0.};
//! Define an inner boundary for the spherical refinement
//! (radius < inner boundary --> L0 refinement).
//! Set to 0 to ignore and refine the entire moon.
constexpr D_REAL SPHERE_INNER_BOUNDARY[] = {0.5, 0, 0};

//! C) Zylinder along x
constexpr INT32 TL_REFINE_STATIC_ZYLINDER_X[] = {-1, -1, -1, -1, -1};

constexpr D_REAL radiusX_refZylinder_of_L[] = {5.* R_Obstacle, 2.5* R_Obstacle};

constexpr D_REAL minX_refZylinder_of_L[] = { -LX, -LX};
constexpr D_REAL maxX_refZylinder_of_L[] = {  LX,  0.};

//! D) Zylinder along y
constexpr INT32 TL_REFINE_STATIC_ZYLINDER_Y[] = {-1, -1, -1, -1, -1};

constexpr D_REAL radiusY_refZylinder_of_L[] = { 0., 0., 0., 0.};

constexpr D_REAL minY_refZylinder_of_L[] = { 0., 0.};
constexpr D_REAL maxY_refZylinder_of_L[] = { 0., 0.};



//! E) Zylinder along z
constexpr INT32 TL_REFINE_STATIC_ZYLINDER_Z[] = {-1, -1, -1, -1, -1};

constexpr D_REAL radiusZ_refZylinder_of_L[] = { 0., 0.,  3.*R_Obstacle};

constexpr D_REAL minZ_refZylinder_of_L[] = {0., 0.,  -1.3333*R_Obstacle};
constexpr D_REAL maxZ_refZylinder_of_L[] = {0., 0.,  +1.3333*R_Obstacle};

//! forbid removal of initial refined blks
constexpr bool never_remove_static_refinement = true;






//! ----------------------------------------------------------------









//! ---- TIME ADAPTIVE REFINEMENT ----------------------------------
//! general remarks:
//! - refcrit = refinement criteria
//! - any field that has an id (see above) can be used
//!   as refinement criteria

//! REFINEMENT PROCEDURE
//! 1)
//! -  for every oct within one block it is tested, whether a
//!    certain threshhold is exceeded.
//! -> if this is true for any refcrit, the oct becomes "falgged for refinement"

//! 2)
//! - if "force_refine_environment" or "force_refine_entire_block" is set, additional
//!   octs will be falgged for refinement.

//! 3)
//! -  afterwards for any existing oct it is tested, whether a
//!    certain threshhold is undercut.
//! -> if this is true for any refcrit, the oct becomes "falgged for removal"

//! 4)
//! a) octs that are flagged for refinement will be refined
//! b) octs that are flagged for removal will be removed
//! c) octs that are not flagged at all will remain unaffected
//! -  In case refcrit A marks an oct for fefinement while refcrit B marks the
//!    oct for removal, the oct will become refinent


//! Apply Mesh Refinement Algorithm every x TL
//! set 0 to switch off
constexpr INT32 TL_REFINE_MESH =   0;


//! specify how many criterias should be used
constexpr INT32 num_refcrit = 5;


//! specify field IDs for all refcrits
constexpr INT32 refcrit_field_IDs[] = {id_rotB, id_UI_plus, id_UI_plus, id_UI_plus, id_UI_plus};

//! component of respective field that should be used as refcrit:
//!  X, Y, Z or MAGNITUDE
//! (0, 1, 2,   3)
constexpr INT32 refcrit_field_comp[] = {3,  2, 2, 1, 1};


//! Choose for each refcrit whether it is maximum or minimum based.

//! A)
//! In case maximum based is choosen, for each oct the following condition will be checked:
//! if(oct_average/global_oct_average_maximum > refine_threshold) flag for refinemnt
//! -> the refine_threshold should be a value between 0 and 1 since the value of
//!    "oct_average/global_oct_average_maximum" will be a value between 0 and 1

//! B)
//! In case maximum based is set false, for each oct the following condition will be checked:
//! if(oct_average/global_minimum < refine_threshold) flag for refinemnt.
constexpr bool refcrit_maximum_based[] = {true, false, true, false, true};


//! decide wether oct_average of respective field shall be normalized
//! by its simulation global extremum of wheter the absolute value
//! shall be used
constexpr bool refcrit_normalize_to_global_extremum[] = {false, false, false, false, false};


//! define refine_threshold for each every refcrit and every level:
//! -> each column lists the refine_thresholds for all levels of one refcrit
//! -> each   row  lists the refine_thresholds for all refcrits of one level
constexpr F_REAL refine_threshold[][NUM_CRITERIA] ={{      0.04,    -0.3,    0.3,   -0.6,    0.6},
                                                {      0.12,    -88.0,    88.0,   -113.0,    113.0},
                                                {      0.5,  -100.0,  100.0, -100.0,  100.0},
                                                {     10.00,  -100.0,  100.0, -100.0,  100.0},
                                                {     10.00,  -100.0,  100.0, -100.0,  100.0}};


//! define remove_threshold for each every refcrit and every level:
//! -> each column lists the remove_thresholds for all levels of one refcrit
//! -> each   row  lists the remove_thresholds for all refcrits of one level
constexpr F_REAL remove_threshold[][NUM_CRITERIA] ={{  0.8 *0.04, -0.8* 0.3, 0.8 *0.3, -0.8*0.6, 0.8*0.6},
                                                {  0.8 *0.12, -0.8* 88.0, 0.8 *88.0, -111.8* 3.0, 111.8 *3.0},
                                                {  0.8 *0.5,     100.0,   -100.0,     100.0,   -100.0},
                                                {  0.8 *1.00,     100.0,   -100.0,     100.0,   -100.0},
                                                {  0.8 *1.00,     100.0,   -100.0,     100.0,   -100.0}};




//! In order to avoid discontonuities on level boundaries environment
//! of flagged block s refined
//! -> In case set false lots of noise is introduced
//! NOTE:
//! -> THIS SHOULD ALWAYS BE SET TRUE, OTHERWISE REFINEMENT BOUNDARIES WILL
//!    BE TO CLOSE TO FEATURES OF INTEREST !!!!
//!    (only set false for testing or quick previews)
constexpr bool force_refine_environment =  true;


//! in case one Block is flagged, every Block of
//! the respective child_array will be flagged
//! resulting in standard block AMR
//! -> should always be set FALSE, except for comparison of
//!    "standard block AMR" vs. "hybrid block AMR"
//! NOTE:
//! WILL NOT BE APPLIED FOR INITIAL REFINEMENT
//! BUT ONLY REFINEMENT DURING RUNTIME
constexpr bool force_refine_entire_block  = false;
//! ----------------------------------------------------------------




//!-------------------------------------------------------------//
//!-------------- 3c) Optimization Parameter -------------------//
//!-------------------------------------------------------------//

//! decide whether to order blocks along
//! Space Filling Curve (SFC) or to use the
//! common linear schema
constexpr bool use_SFC = false;


//! in order to distribute blocks to intervals of SFC,
//! processing time of each block is recorded.
constexpr INT32 TL_reset_block_timing = 0;

//! in order to obtrain imroved massload distribute
//! blocks among mpi processes
//! Attention: Only some fields will be redistrubuted!!! 
//! Problems with Neutral gas velocity and electron pressure!!!
constexpr INT32 TL_REDISTRIBUTE_BLOCKS = 0;


//! distribute blocks with some time level offset
//! -> this should make the code redist blocks
//!    some time levels after refinement took
//!    place
//! (for some reason crashes are more likely
//! (when values different from zero are used ???)
//! (-> better set to zero)
constexpr INT32 OFFSET_REDISTRIBUTE_BLOCKS = 0;

//! decide to assign children to respective parent or not
constexpr bool distribute_RB_based = false;

//! force redistribution aof blocks after restore
constexpr bool redistribute_after_restore = true;

//! choose redistibrution criteria:
//! available are:
//! TOTAL_TIME
//! FIELD_TIME
//! BLOCK_NUMBER
//! PARTICLE_TIME
//! PARTICLE_NUMBER

//! below should be used for root block based distribution
//! TOTAL_TIME_INC_CHILDREN
//! FIELD_TIME_INC_CHILDREN
//! BLOCK_NUMBER_INC_CHILDREN
//! PARTICLE_TIME_INC_CHILDREN
//! PARTICLE_NUMBER_INC_CHILDREN
constexpr INT32 distribution_criteria = PARTICLE_NUMBER;

//! Synchronisation of MPI Sends and Recieves
//! true = send -> recieve -> wait
//! false = send -> recieve -> wait_all
constexpr bool sync_mpi_send_rec = true;

//! Some processes will operate many blocks with few particles each, other
//! processes will operate few blocks with many particles each.
//! 
//! It may happen that inside block XY are as many particles as in
//! 100 other blocks. Hence one process would deal with block XY
//! where another process would operate all 100 other blocks which leads
//! to a bad load balance at the magnetic field routine.

//! This can be avoided by limiting the maximal number of blocks/each_process
//! to a multiple of the average block number each process (num_blocks/num_processes)

//! (this parameter is ignored in case BLOCK_NUMBER is
//! specified as distribution_criteria)
constexpr F_REAL multiple_of_average = 2.;


//! it is iteratively estimate how many block a given process has to operate,
//! depending on the number of particles inside this block
//! The higher the ireration the better the distribution should be, of course
//! it will sature for very high numbers (>100000)
constexpr INT32 max_iteration_redistribute = 400000;


//! show global information (particles moved, collected, moved ...)
//! since all information have to be gathered collecively
//! together which is expensive, this SHOUld not be done to frequently
constexpr INT32  TL_LOGFILE_globalINFO = 1;

//! decide whether each process should write its own log file
//! (in general only useful for debugging)
constexpr bool  show_proc0_logfile_only = true;

//! one particle_array is allocated each cell.
//! -> In case particles enter the cell but the array is to small, a new array is
//!    allocated, particles are copied into it and the old array is deleted
//! -> in order to avoid to frequent re-allocation, min_pArray_Size_factor*num_particle
//!    array elemts are allocated when ever the particle array is to small
constexpr D_REAL min_pArray_Size_factor = 1.1;

//! - In order to avoid to high memory consumption, a small array is allocated in
//!   case num_elements*max_pArray_Size_factor > num_particle
//! - particles are copied and the large array deleted
constexpr D_REAL max_pArray_Size_factor = 1.3;

//! specify how often the it should be checked whether the particle arrays are too large
constexpr INT32 TL_RESIZE_pARRAYS = 1;


//!-------------------------------------------------------------//
//!-------------- 3d) Simulations Configuration ----------------//
//!-------------------------------------------------------------//

//! These variables allows a fast configuration of A.I.K.E.F. for simple test, e.g. test of the ionisation routine
//! This can be used for an automatic test routine. ...
constexpr bool run_calc_first_E = true;
constexpr bool run_CAM = true;
constexpr bool run_calc_second_E = true;
constexpr bool run_accelerate_Particle = true;
constexpr bool run_collect_Ui_minus = true;
constexpr bool run_move_Particle = true;
constexpr bool run_Split_Merge_Particle = true;
constexpr bool run_negative_Particles = false;
constexpr bool run_inject_obstacle_ions = use_atmosphere;
constexpr bool run_collect_Rho_prepare_Recombination_Density = true;
constexpr bool run_chemical_Reactions = false;
constexpr bool run_resize_pArrays = true;
constexpr bool run_collect_Rho_calc_Recombination_Density = true;
constexpr bool run_average_Ui_rho_setREZrho = true;
constexpr bool run_advanceB = true;


//!-------------------------------------------------------------//
//!-------------- 4) Ion / Ionospheric Parameter ---------------//
//!-------------------------------------------------------------//

//ionization switch
constexpr INT32 TLIO = 0;
constexpr bool ionizationonly = true;

//! Total Number of ion and dust species (NOT electrons), has to be at least one.
//! ("Total" meaning Inflow-Ions plus Obstacle-Ions)
constexpr INT32 num_Charged_Species = NUM_CHARGED_SPECIES;

//! Total number of species that are treated as PARTICLES
//! may be different from num_Charged_Species in case
//! a dust (or ion) species is added to the density from
//! some field
//! NOTE: Has to be set above
constexpr INT32 num_Particle_Species = NUM_PARTICLE_SPECIES;  

//! Number of ion species emmitted from each Box-Boundary.
//! In case set to zero, no plasma will be Initialzed.
//! Still obstacle ions can be inserted, which is a 
//! good (and especially fast) opportunity to find out,
//! how an obstacle ion profile develops.
constexpr INT32 num_Inflow_Species = 1;

//! Indices of inflow species.
//! Exactly "num_Inflow_Species" integers have to be provided.
//! Be carefull not to use the same species for inflow- and obstacle ion
//! at the same time when editing "CHybrid_IonProfiles.cpp".
//! Indices starting at 0 !!!
constexpr INT32 index_Inflow_Species[] = {0};


//! define when start to inject respective species
//! - will be ignored in case is no inflow species
//! - set to zero for continues inflow
constexpr D_REAL start_inject_species_each_xx_t0[] = {0 ,0, 0, 0, 0};

//! define duration of injection
//! - will be ignored in case start is zero
constexpr D_REAL duration_inject_species_each_xx_t0[] = {0 ,0, 0, 0, 0};


//! Special Velocity Distribution
//! true: ions which will injected at the boundaries or into empty cell will have a special distribution in velocity space. (c.f. CBlk_Part-Init.cpp -> fill_empty_Cell  for more details)
//! false: ions will have a normal maxwellian velocity distribution
constexpr bool special_Velocity_Distribution[] = {false, false, false, false};


//! CONCERNING BRACKETS BELOW:
//! Only the first entry up to num_Particle_Species will be used,
//! all other values will be ignored.

//! Specify Mass, Charge, density, beta/thermal velocity and Obstacle rho for all species:
//!   Mass of 1. denotes one ion Mass of the normalization species (eg. 1 proton mass)
//! Charge of 1. denotes one ion Charge of the normalization species (eg. 1 proton Charge)

//! SPECIES 0 is reference species (!!!)


constexpr D_REAL  Ion_Masses[] = {
    14./14.,
    32./SI_m0*m_p, //O2
    2./SI_m0*m_p,  // H2
    16./SI_m0*m_p // H2O
    //1.
};

constexpr D_REAL Ion_Charges[] = { 
    1.,
    1.,
    1.,
    1.
  };
			     
//! Plasma Densities for each inflow species.
//! (values at indices of obstacle species will be ignored)
constexpr D_REAL rho_sw[] = {0.1/0.1, 0.,  0.,  0.,  0., 0.,0.,0.,0.};



//! Ion plasma betas for each Species:
//! THIS VARIABLE MAY ONLY BE USED, IN CASE "rho_sw" AND "Ion_Masses" ARE
//! SPECIFIED FOR THE RESPECTIVE SPECIES (as mass-density is needed to derive
//! thermal velocity from plasma beta).
//! In case "rho_sw" or "Ion_Masses" are not set for respective
//! species (e.g. in case of obstacle ions), directly specify
//! temperature by using "Ti" variables.
//! Number is ion temperature
constexpr D_REAL Ion_Betas[]  =  { 
    calcBeta*100*e/kB,
    0.,//0.1116,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.
  };

//! allow for different temperature parallel and perpendicular to B
//! TODO currently it is assumed that B0 points in z direction...
//! "temperatures" have to be inserted in Joule
constexpr D_REAL Ti_para[] = {0*e,0*kB*180*e,0.,0.,0.,0.,0.,0.,0.,0.};
constexpr D_REAL Ti_perp[] = {0*e,0*kB*180*e,0.,0.,0.,0.,0.,0.,0.,0.};


//!electron temperature in K
constexpr D_REAL Te = 100*e/kB;
constexpr D_REAL Te_ionosphere = 0.5*e/kB;    
//! electron plasma betas for each Species.
// constexpr D_REAL Electron_Betas[]  = {calcBeta*Te, 0., 0., 0.};
constexpr D_REAL Electron_Betas[]  = {calcBeta*Te, 0., 0.,0.};

//! Adiabatic exponent for equation of state
constexpr D_REAL kappa_electron =  2.0;



////! Relevant if either nonadiabatic_gradPE_TERM or use_neutral_species_as_field are true
//! Number of neutral species:
constexpr INT32 num_Neutral_Species = NUM_NEUTRAL_SPECIES;

//! use analytical expression for neutral density
//! Profile is defined in function neutral_profile() in parameters.cpp
constexpr bool analytical_neutral_profile = true;

//! use neutral density profile from file
constexpr bool neutral_profile_from_file = false;
//! NOTE parameters above usually exclude each other. However, there may be
//! rare cases in which both, analytical profile and profile from file are used.
//! In that case, you should be careful which is set first...

//! Particle-masses of neutral species (total number of species is defined above)
// constexpr D_REAL Neutral_Masses[] = { 18*amu/SI_m0, 1., 1., 1. };
constexpr D_REAL Neutral_Masses[] = {Ion_Masses[1]/m_p*amu, Ion_Masses[2]/m_p*amu, Ion_Masses[3]/m_p*amu};//,1.};//{44.*amu/SI_m0, 32.*amu/SI_m0, 1.};

//! Beta of newly ionised electrons
constexpr D_REAL NewElectron_Betas[]  = {calcBeta*0*e/kB, 0*e/kB*calcBeta, 0*e/kB*calcBeta, 0.*calcBeta,0.*calcBeta,0.*calcBeta,0.*calcBeta,0.*calcBeta,0.*calcBeta,0.*calcBeta,0.*calcBeta};
// constexpr D_REAL NewElectron_Betas[]  = {0.*calcBeta*100., 0*calcBeta*180., calcBeta*180., calcBeta*180.};

//! temperature of neutral gas in Kelvin!
constexpr D_REAL Neutral_Temperature[] = {0,0,0,0};

//! degrees of freedom of neutral species
//! for example 3 for a pointmass
//! water = 12 see Page 120
constexpr D_REAL Neutral_DOF[] = { 3., 3., 3., 3. ,3.,3.,3.,3.,3.,3.,3.};

//! below necessary if neutral field is read from file
constexpr D_REAL norm_externNeutralDensity[] = {1.,1.,1.,1.};

constexpr D_REAL norm_externNeutralVelocity[] = {1.,1.,1.,1.};

//! file names for reading extern fields
constexpr char extern_NeutralField_name[][80] = { "../../../1236.dat",
			"lineY",
			"lineZ"};
							
////!


//! -------------- Reactions, collisions, charge exchange and recombination --------------------------------------------------------
							
//! Reaction Rates for Reactions of the type							
//! spec X + NeutralSpecies Y1 -> DestSpec Z + NeutralSpecies Y2
//! ReactionRate[X][][] := Source Ion Species
//! ReactionRate[][Z][] := Destination Ion Species
//! ReactionRate[][][Y] := Neutral Species
//! note that reactions include elastic collisions and charge exchange
//! array entries have to be in normalized units
//! since rates are usually given in cm^3/s, normalize by multiplying with	
constexpr D_REAL RateNorm = SI_t0*1e-6*SI_n0;  //! NOTE need dt??
constexpr D_REAL ReactionRate[][num_Particle_Species][NUM_NEUTRAL_SPECIES] =
  //each row below is a desitination species. for example:
  // S 0:
  // {
  //   { X, X }  //destination is S0 (if S0 is above, this would be elastic collsions)
  //   { Y, Y }  //destination is S1 (if S0 is above, this would be charge exchange (to S1)
  //   { Z, Z }  //destination is S2 (is S0 is above, this would be charge exchaneg (to S2)
  // }
  //
  // then would come blocks for S1, S2,... (source ion species)

  // -------------------------- Upstream Species: -------------------------- 
  //O+ (Current Sheet)
  {/* // Species 0    
    {
      {ateNorm*8.95e-10,RateNorm*6.64e-10,0.},                 // {O+ + CO2 --> O+,   O+ + O2 --> O+}   
        for S0=H2: // {RateNorm*29.e-10, RateNorm*22.5e-10},  // see S&N eqn 4.88 (check units if re-using eqn!)
      {0.,0.,0.},                                               // {O+ + CO2 --> CO2+, O+ + O2 --> CO2+}
      {0.,0.,0.},                                                // {O+ + CO2 --> O2+,  O+ + O2 --> O2+}
      {0.,0.,0.}
    },
    
    // ----------------------- Callisto Atmosphere: -------------------------- 
    CO2/
    {
      {0.,0.,0.},                                               // {CO2+ + CO2 --> O+,   CO2+ + O2 --> O+}
      {RateNorm*2.39e-10, 0.,0.},                                // {CO2+ + CO2 --> CO2+, CO2+ + O2 --> CO2+}
      {0.,0.,0.},                                                // {CO2+ + CO2 --> O2+,  CO2+ + O2 --> O2+}
      {0.,0.,0.}
    },

    //O2
    / Species 2 */
        
    {
      //O2+}
      {0.},                                  // {O2+ + CO2 --> O2+,  O2+ + O2 --> O2+}
      {0.},
    },
        
    {
      {0.},
      {0.},
    },
  };

//! apart from the constant rates in the array above, it is also possible to include a velocity-dependence of the reaction rate
//! the expressions for the velocity dependence have to be provided in 
constexpr bool use_velocity_dependent_reaction_rates = false;
								
//! activate recombination for respective ion species
//! recombination rate is calculated in parameters.cpp
//! in function void CBlock::calc_RecombinationAlphaField(void)
constexpr bool Recombination_for_Species[] = {0, 0, 0, 0, 0, 0,0,0,0};

//! check max value of reaction rates
constexpr bool check_max_reaction_rates = true;
//! check max value of p = dt/tau
constexpr bool check_max_reaction_probability = true;

//! --------------------------------------------------------------------------------------------------------------------------------


//!----------------- Photoionisation, eletron impact ionization and stellar wind production ----------------------------------------

//! two ways to describe ionization processes:
//! a) specify ionisation rate if neutral density is set
//	advantage: no global value for production rate has to be known
//! b) specify total ion production rate

//! Characteristic Frequency of Photoionisation in 1/s
//! normalization with SI_t0
// constexpr D_REAL PhotoionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	// //	   /*Neutralspec 0*/, /*Neutralspec 1*/ ----
	// /* Species 0*/ /*{*/    6.4e-11*SI_t0   ,  /* 0    },*/
	// /* Species 1*/ /*{*/    6.1e-10*SI_t0   ,  /* 0    },*/
	// /* Species 2*/ /*{*/    3.7e-09*SI_t0   ,  /* 0    },*/
	// /* Species 3*/ /*{*/     0		,  /* 0    },*/
	// /* Species 4*/ /*{*/    1.4e-10*SI_t0  /* ,   0    }*/
	// };
constexpr D_REAL PhotoionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	//	   /*Neutralspec 0*/
	/* Species 0*/ /*{*/    0.0     /* 0    },*/
	};

//! value comet Churyumov-Gerasimenko	
// constexpr D_REAL COM_nu = 5.88e-7;
	
//! values Enceladus solar maximum	
// constexpr D_REAL PhotoionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
// 	//	   /*Neutralspec 0*/, /*Neutralspec 1*/ ----
// 	/* Species 0*/ {    2.4e-10*SI_t0   ,   0    },
// 	/* Species 1*/ {    1.6e-09*SI_t0   ,   0    },
// 	/* Species 2*/ {    9.1e-09*SI_t0   ,   0    },
// 	/* Species 3*/ {     0		    ,   0    },
// 	/* Species 4*/ {    4.5e-10*SI_t0   ,   0    }
// 	};	
				
	
//! decide if electron ionisation rate depends on electron temperature
//! NOTE not yet implemented for electron pressure equation!
constexpr bool ElectronionisationRate_Te_dependent = false;
//! if false, normalization with SI_t0
//! if true, ElectronionisationRate[][] is multiplied with value of
//! field id_ElectronTemperature (Te in Kelvin)!
//! if expression for ElectronionisationRate is more complex (e.g. ce_sqrt(Te)*const)
//! this has to be changed in CBlk_IonProfile.cpp in insert_ions_from_neutral_profile
// constexpr D_REAL ElectronionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	// //	                           /*Neutralspec 0*/, /*Neutralspec 1*/ ----
	// /* Species 0*/ /*{*/ Ence_nu_e_total*0.007/1.26*SI_t0   ,  /* 0    },*/
	// /* Species 1*/ /*{*/ Ence_nu_e_total*0.222/1.26*SI_t0   ,   /*0    },*/
	// /* Species 2*/ /*{*/ Ence_nu_e_total*0.958/1.26*SI_t0   ,   /*0    },*/
	// /* Species 3*/ /*{*/ 		0 		    ,   /*0    },*/
	// /* Species 4*/ /*{*/ Ence_nu_e_total*0.0759/1.26*SI_t0  ,   /*0    }*/
	// };	
constexpr D_REAL ElectronionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	0.0
	};	
	
	
//! this replaces the old and ugly parameter obstacle_MP_weight_each_t0[]
//! entries in ions/s
constexpr D_REAL Global_IonProduction_Rate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	//	     /*Neutralspec 0*/, /*Neutralspec 1*/ ----
	/*Species 0*/ /*{*/    0.0    ,    /*0.0     },*/
	/*Species 1*/ /*{*/    0.0    ,    /*0.0     }*/
	};	


//! num macro Particle each TL
//! representing the heavy ions
//! NOTE: if ions are inserted by ionization of neutral profile,
//! the number of inserted ions is adjusted to this value6
constexpr INT32 obstacle_MP_num_each_TL[] = {     
    0,
    300,
    300,
    300
};

		

////! relevant if nonadiabatic_gradPE_COLLISION_TERM is true

  //! collision parameter for electron-neutral-interactions
  //! the product of the electron density (ion charge density in normalized units)
  //! and the neutral density
  //! and the en_coll_param is the factor
  //! n_e*m_e /(m_e+m_n) * nu_{e,n}
  //! -> see Bachelorthesis (Hendrik Ranocha)
  //! since rates are usually given in cm^3/s, normalize by multiplying with RateNorm               
  //! implementation requires addition factor m_e /(m_e+m_n)                                          
  constexpr D_REAL en_coll_param[] = {                                                                                                                          
      0.* 4.427e-11*RateNorm*mass_electron_norm/(mass_electron_norm+Neutral_Masses[0]), 
      0.,
      0.,
      0.,
      0.
    };

////! end relevant for nonadiabatic_gradPE_COLLISION_TERM  
  
  
						    
////! relevant if use_dust_species_as_field is true

	//! Number of neutral species:
	//! used as parameters in nonadiabatic_gradPE_TERM computation
	constexpr INT32 num_Dust_Species = 0;

	//! NOTE parameters should be renamed
	//! use analytical expression for dust density, but add this density to density from file
	constexpr bool analytical_dust_plume = false; 	
	//! use analytical expression and no file
	constexpr bool analytical_dust_plume_only = false; 

	//! below necessary if dust field is read from file

	constexpr D_REAL norm_externDustDensity[] = {1.,1,1};

	constexpr D_REAL norm_externDustVelocity[] = {1.e-03,1,1};

	//! file names for reading extern fields
	constexpr char extern_DustField_name[][80] = { "../../../StaubFeld10B.dat",
							"lineY",
							"lineZ"};
							
	//! Full dust density at TL
	constexpr INT32 TL_DUST_MAX = -1;//40000*(1.-0.5*LOWRES)*(1.- 49./50.*VOYAGER);
							
	//! if density corrections according to dusty plasma shall be considered, potential of 
	//! dust cloud may be included in extern rho file
	constexpr bool phi_cloud_from_extern = false;						    
							
	//! smooth dust field
	constexpr INT32 num_smooth_extern = 0;
	constexpr D_REAL smooth_extern[] = {0.0,0.,0.,0.4,0.8};

////! end relevant for use_dust_species_as_field


	
////! Relevant if use_ion_production_as_field is true
	
	//! ion production field is photoionization frequency * neutral density
	//! due to shadow or absorption of photons (chapman profile), this may be
	//! different from the neutral density profile
	//! the neutral density itself is, however, necessary for the electron pressure
	//! and therefore two different fields are necessary

	//! use ion_prod density profile from file
	//! If this is false, you MUST fill in the electron ionization array above: ElectronionisationRate
	//! If an external file is used, all ion production should be there (electron, photoionization, etc)
	//! Currently internal and external production cannot be combined without modification
	constexpr bool ion_prod_profile_from_file = true; 
	
	// -------- NOTE: --------
	//If using relative scale heights, UPDATE parameters.cpp!!!!!
	// ----- SEE ABOVE! ----

	//! file names for reading extern fields
	constexpr char extern_IonProdField_name[][80]= { "../data/atm_O2", "../data/atm_H2", "../data/atm_H2O"};

	//! calculate the ion production profile of the neutral field
	//! that results from photoabsorption (Chapman profile)	
	constexpr bool calc_ionProd_from_neutral_field = false;

	//! set analytical ion production profile to field
	//! e.g. analytical Chapman profile or modify neutral profile
	//! to consider shadow of moon or planet
	//! This will use set_ion_production_profile defined in parameters.cpp
	constexpr bool set_analytical_ionProd_profile = false;
	
	constexpr D_REAL norm_externIonProdRate[] = {1.,1.,1.};

	constexpr D_REAL norm_externIonProdVelocity[] = {1.,1.,1.};

////! end relevant for use_ion_production_as_field


//!--------------------------------------------------------
//!----- extern fields related  --------------------------
//!---------------------------------------------------------

//! This parameter applies to both, neutral and dust fields
constexpr bool serialize_reading_of_extern_field = true;

//! Extern Fields for future usages
//! NOTE extern fields are currently included as neutral fields, ion production fields 
//! and dust fields. However there may be future usage as fields for multi-scale simulations.
//! TODO Horia...

//! Indices of inflow species.
//! Exactly "num_externRhoVelocityFields" integers have to be provided.
//! Be carefull not to use the same species for inflow,obstacle and extern ions
//! at the same time.
constexpr INT32 index_externRho_Species[5] = {};

//! normalization of extern Density:
//! e.g. Mercury:
//! n0 = 32/cm^3
constexpr D_REAL norm_externRho[] = {12.,1,1};

//! normalization of extern Velocity:
//! e.g. Mercury:
//! B0 = 21nT
//! n0 = 32.e6 1/m^3
//! -> v0 = 80.9728 km/s = 8.09728*10^6 cm/s
constexpr D_REAL norm_externUi[] = {1.82501515*8.09728 *1.e6,1,1};;

//! file names for reading extern fields
constexpr char extern_Field_name[][RUN_NAME_SIZE] = { "/velocity_ionization_field.txt",
                                                    "lineY",
                                                    "lineZ"};




//!-------------------------------------------------------------//
//!-------------- 6) Obstacle Parameter ------------------------//
//!-------------------------------------------------------------//

//! decide if use resistivity inside obstacle
//! details are specified in CBlk_EtaProfiles.cpp
constexpr bool use_resistive_obstacle = true;

//! Obstacle Resistivity
	//! (if Diffusion shall compensate lack of convection within
	//!  the obstacle, ETA has to be larger than V_sw*2*R_Obstacle)
	//!  NOTE: either 
	//!        1) switch on "advance_obstacle_B"
	//!           when using high eta values (>>1.)
	//!           -> Crank Nicholson discretisation for BField
	//!     or 2) set ETA_TERM_BField to true
	//!     	 above for low eta values (~1.)
	//!           -> leap Frog discretisation for BField
constexpr D_REAL Eta_Obstacle = 1./(0.5*1./(mu_0*SW_v*R_Moon*SI_x0))*e*SI_n0/SI_B0*10;

//! define radius in which Eta_Obstacle shall be used
constexpr D_REAL R_Eta_Obstacle = 0.96*R_Obstacle;

//! Either use fermi profile or smoothing for transition from
//! Eta_Obstacle to background eta (Eta_sw above)
//! Turn this off if using a comet eta profile (below)
constexpr bool use_eta_fermi = true;
//! fermi profile for resistivity given by
//! Eta= Eta_Obstacle *1./(exp((r - R_Eta_Obstacle)*fermi_slope)+1.);
constexpr D_REAL fermi_slope = 30.;

//! Decide if a inverse fermi profile should be used
constexpr bool use_eta_comet = false;
//! Resistivity in the solar wind at infinity 
//! the value of Eta_comet_innerComa will be added
constexpr D_REAL Eta_comet_inf = 3.e-1;
//! Resistivity in the inner Coma
constexpr D_REAL Eta_comet_innerComa = 1.e-3;
//! Determine the radius of low resistivity area
constexpr D_REAL R_Eta_innerComa = 4.;
//! slope of change
constexpr D_REAL fermi_slope_comet = 0.025;

//! Optionally set an intermediate core to bridge the resistivity gap
//! between the Fermi plateau and the numerical core.
//! Extent of the numerical core determined by bridge_profile_start
//! and bridge_profile_end
constexpr bool use_intermediate_core = true;
//! This is where the intermediate profile will start:
constexpr D_REAL bridge_profile_start =  0.5;
//! This is where the intermediate profile will connect to the Fermi profile
constexpr D_REAL bridge_profile_end = 0.7;

//! within the core no magnetic field is advanced nor smoothed
//! (-> it will remain the initial field forever)
//! in % of obstacle radius
constexpr D_REAL obstacle_core_fraction = 0.3 * use_dipole;

//! for two body simulations: activate second obstacle
constexpr bool use_second_obstacle = false;
constexpr D_REAL R_SecondObstacle = 0;
constexpr D_REAL Position_SecondObstacle[3] = {10,0,0}; 


//! OBSTACLE MAGNETIC Field calculation.
//! OM1) Switch on/off calculation inside the Obstacle.
//! TODO: does not work with AMR up to now
constexpr bool advance_obstacle_B = true;

//! OM2) choose relaxations parameter for SOR
constexpr D_REAL B_SOR_omega_L[] = {1.1, 1.15, 1.2, 1.25};
//{1.2, 1.3, 1.35, 1.4};

//! OM3) choose break criteria for SOR
constexpr D_REAL B_SOR_max_error = 1.e-4;

//! OM4) B_SOR_num_cycle_base^L cycles will be performed
//!      in level L each iterations
constexpr INT32 B_SOR_num_cycle_base = 2;

//! OM5) choose maximal number of iterations for SOR
constexpr INT32 B_SOR_max_Iteration = 2000;

//! OM6) calc error every xx TL
constexpr INT32 B_SOR_calc_error_step = 3;


//! DIVERGENCE CLEANER (DC)
//! DC1) Switch on/off BField-divergence cleaning inside the Obstacle.
constexpr bool div_cleaner = true;

//! DC2) choose relaxations parameter for SOR
constexpr D_REAL DC_omega_L[] = {1.5, 1.6, 1.7, 1.8};

//! DC3) choose break criteria for SOR
constexpr D_REAL DC_max_error = 1.e-4;

//! DC4) B_SOR_num_cycle_base^L cycles will be performed
//!      in level L each iterations
constexpr INT32 DC_num_cycle_base = 2;

//! DC5) choose maximal number of iterations for SOR
constexpr INT32 DC_max_Iteration = 2000;

//! DC6) calc error every xx TL
constexpr INT32 DC_calc_error_step = 3;


//! The relaxation parameters B_SOR_omega_L[] and DC_omega_L[]
//! crucially determine the convergence speed of the solvers.
//! But they can be determined only empirically for each simulation setup.
//! Thus check which value yields fastest convergence (smallest error)
//! testing values from B_SOR_omega_L-steps/2 in steps of 0.01
//! set to zero to switch optimization off
constexpr INT32 num_optimize_SOR_steps = 0;
//! TL starting the optimization 
//! At TLstart_optimize_B_SOR_omega + num_optimize_SOR_steps, the
//! optimal value will be written in the logfile and this value
//! should be used for next simulation.
constexpr INT32 TLstart_optimize_B_SOR_omega = 10;




//! Switch on/off electric field calculation inside the Obstacle.
//! should usually be switched on. In case of inert moons end high
//! resistivities its more stable to switch it of.
constexpr bool advance_obstacle_E = true;


//! Obstacle's intrinsic Dipole Field
//! In order to fix it inside obstacle use
//! a obstacle_core_fraction>0.
//! MM is defined such that at the origin MM and BField are parallel.
constexpr D_REAL Magnetic_Moment[3] = {-2.910e18/SI_M0*use_dipole, 7.170e18/SI_M0*use_dipole, -1.310e20/SI_M0*use_dipole};
constexpr D_REAL MM[3] = {SI_M0*Magnetic_Moment[0], SI_M0*Magnetic_Moment[1], SI_M0*Magnetic_Moment[2]};

//! if Magnetic Moment is from obstacle, then the corresponding field
//! MUST NOT be set at ghostnodes! In case it is extern (such as Saturn's
//! dipole field for Enceladus) ghost nodes have to be initialized with 
//! the dipole field
constexpr bool is_intern_Magnetic_Moment = true;

//! Tilt angle in degree
constexpr D_REAL Magnetic_Moment_Angle[2] = {0., 0.};

//! Position Pffset in normalized units
// constexpr D_REAL Magnetic_Moment_Offset[3] = {-6., 5., 0.0};
// constexpr D_REAL Magnetic_Moment_Offset[3] = {0., 0., +10.06};
constexpr D_REAL Magnetic_Moment_Offset[3] = {0., 0., 0.};

//! specify angular velocity for obstacle
constexpr D_REAL omega_rotating_obs[3] = {0., 0., 0.};



//! densities for each ion species to eliminate strong electron pressure 
//! Gradients at Obstacle surface.
//! - inner densities are set INTERNALLY for gradPE calculation
//!   and NOT visible in visualization!
constexpr D_REAL obstacle_rho[] = { MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField};


//! For logging
constexpr D_REAL ram = SI_m0 * SI_n0 * SW_v*SW_v;
constexpr D_REAL thermal = SI_n0*kB*(Ion_Betas[0]/calcBeta+ Te);
constexpr D_REAL mag = SI_B0 * SI_B0 / (2*mu_0);


//! Abreviations below were originally defined in absolute_Globals.h.
//! In rare cases this leads to incorrect values, since the order
//! of value assinement is not well defined.
//! Defining these parameters at the end of this file always yields
//! correct assinement.
constexpr INT32 num_root_blocks    = RB_X*RB_Y*RB_Z;
constexpr INT32 num_nodes_in_block = BlkNds_X *BlkNds_Y *BlkNds_Z;

