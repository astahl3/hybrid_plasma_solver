#pragma once

#include "types.h"

//! perform aditional interrogations to check for rarely occuring errors
constexpr bool DEBUG_MODUS = false;

//! real Inprecission interception:
//! used in move and insertion procedure
//! should be 1.e-6 in case PARTICLE_REAL is float
//! should be 1.e-16 in case PARTICLE_REAL is double
//! (see comments in CBlock::move_particle for details)
constexpr double PART_REAL_PRECISION = 1.e-16;

#define BUILD_SUM	0
#define BUILD_MAX	1

#define UNIFORM		0
#define STAGGERED	1

//!-------------------------------------
//!-- PHYSICAL FIELDS FOR CALCULATION --
//!-- (-> allocate memory for each)   --
//!-------------------------------------
//! EM Field ids
#define id_BTotal	0
#define id_BEven	1
#define id_BOdd	    2
#define id_Bcfbg	3
#define id_EField   4


//! density Field ids
#define id_Lam		5
#define id_rho_n	6
#define id_rho_np1	7
#define id_rho_rez	8
// #define id_allRhoSpecies 	 9

//! ion mean velocity Field ids
#define id_Gam		9
#define id_UI_plus	10
#define id_UI_minus	11

//! misc field id's
#define id_Eta 			12
#define id_PhiDC		13
#define id_BDerivative	14

//! Field Groups
#define id_ALL_FIELDS			15
#define id_UIp_Gam_allRhoSpec	16


//!------------------------------------------
//!-- FIELDS EXCLUSIVELY FOR VISUALIZATION --
//!--   (-> share one "scratch" field )    --
//!------------------------------------------

//! derived physical fields
#define id_gradPI0	    17
#define id_rotB 		18
#define id_divB 		19
#define id_divE 		20
#define id_PMagnetic	21
#define id_PTotal		22

//! time 
#define id_FieldTime		23
#define id_ParticleTime	    24
#define id_scratch_vector	25
#define id_scratch_scalar	26

//just used as scratchpad by HK-Cod
#define id_ExB	25

//! refinement
#define id_gradB					27
#define id_Refine_Rating			28
#define id_ElectronTemperature	    29
#define id_rho_np1_recombined		30

#define id_Bcfbg_HIMesh			    31

//! species related
#define id_allRhoSpecies        32
#define id_allPISpecies         33
#define id_allGyroSpecies       34
#define id_allGyroELSpecies     35
#define id_allPESpecies         36
#define id_allUISpecies         37
#define id_allForcesSpecies     38

#define RUN_NAME_SIZE	40
#define NUM_CRITERIA	5

//! Box - Boundarie Abbriviations
#define Xmin_BB	 0
#define Xmax_BB	 1
#define Ymin_BB	 2
#define Ymax_BB	 3
#define Zmin_BB	 4
#define Zmax_BB	 5

#define _Im1_		0
#define _Ip1_		1
#define _Jm1_		2
#define _Jp1_		3
#define _Km1_		4
#define _Kp1_		5
#define _ANY_		6
#define _INIT_		6
#define _OBS_		7

//! RESTORE STATE CONVERSION:
#define LEAVE_UNCHANGED      0
#define CONVERT_TO_TRACK     1
#define CONVERT_TO_ORDINARY  2

//! Block redistribution
#define NUM_REDIST_OPTIONS			 11
#define AVERAGE_REF_VALUES			 0
#define TOTAL_TIME					 1
#define FIELD_TIME					 2
#define BLOCK_NUMBER				 3
#define PARTICLE_TIME				 4
#define PARTICLE_NUMBER				 5

#define TOTAL_TIME_INC_CHILDREN		  6
#define FIELD_TIME_INC_CHILDREN		  7
#define BLOCK_NUMBER_INC_CHILDREN	  8
#define PARTICLE_TIME_INC_CHILDREN	  9
#define PARTICLE_NUMBER_INC_CHILDREN  10

//! uniform_grid output
#define uniform_grid_TYPE_BFIELD			0
#define uniform_grid_TYPE_EFIELD			1
#define uniform_grid_TYPE_DENSITY_FIELD	    2
#define uniform_grid_TYPE_VELOCITY_FIELD	3

//! Define use_vectorclass if you want to use the optimized
//! vectorclass from http://www.agner.org/optimize/#vectorclass
//!
//! If use_vectorclass is true,
//! NOTE: num_nodes_in_block has to be a multiple of 4
//! because of CBlock::sub_squaredField_Multiply
#define use_vectorclass false

#if (use_vectorclass)

//! Use Intels Short Vector Math library (SVML)
//! see vectorclass/VectorClass.pdf, page 32 etc.
//! for details about it
//! TODO: Does not work correctly up to now...
#define VECTORMATH 2

#include "vectorclass/vectorclass.h"
#include "vectorclass/special/vectormath.h"

//! NOTE: If you change something in the other typedefs,
//! you have to make the same changes here!
typedef Vec4d VEC4_PARTICLE_REAL;
typedef Vec4f VEC4_F_REAL;
typedef Vec4d VEC4_D_REAL;
typedef Vec2d VEC2_D_REAL;

#endif


//! ----------------------------------------------------------------------


//! to run on NESS and HECToR below has to be defined
#define MPICH_IGNORE_CXX_SEEK


constexpr bool DO_PARENT_BUFFER_UPDATE = true;
constexpr bool NO_PARENT_BUFFER_UPDATE = false;

constexpr bool INCL_PHYS_FACE = true;
constexpr bool SKIP_PHYS_FACE = !INCL_PHYS_FACE;

//! The more Ion Species are allocated, the more
//! Fields must be allowed

#if (NONADIABATIC_GRADPE_TERM || USE_NEUTRAL_SPECIES_AS_FIELD) && USE_DUST_SPECIES_AS_FIELD && USE_ION_PRODUCTION_AS_FIELD

#define NUM_FIELDS 114
#define NUM_SCALAR_FIELDS_EACH_BLK (49 +13*NUM_CHARGED_SPECIES + NUM_PARTICLE_SPECIES+ 6*NUM_NEUTRAL_SPECIES+ 4*num_Dust_Species+5*num_ion_prod_fields+4*num_externRhoVelocityFields +num_scalar_average_fields)

#elif (NONADIABATIC_GRADPE_TERM || USE_NEUTRAL_SPECIES_AS_FIELD) && USE_DUST_SPECIES_AS_FIELD && !USE_ION_PRODUCTION_AS_FIELD

#define NUM_FIELDS  114
//! TODO NOTE I have absolutely no clue why it has to be 5*num_Dust_Species instead of 4*num_Dust_Species, but 4 produces a seg fault
#define NUM_SCALAR_FIELDS_EACH_BLK (49 +13*NUM_CHARGED_SPECIES + NUM_PARTICLE_SPECIES+6*NUM_NEUTRAL_SPECIES+ 5*num_Dust_Species+4*num_externRhoVelocityFields +num_scalar_average_fields)

#elif (NONADIABATIC_GRADPE_TERM || USE_NEUTRAL_SPECIES_AS_FIELD) && USE_ION_PRODUCTION_AS_FIELD  && !USE_DUST_SPECIES_AS_FIELD

#define NUM_FIELDS 200
//! TODO NOTE I have absolutely no clue why it has to be 5*num_ion_prod_fields instead of 4*num_ion_prod_fields, but 4 produces a seg fault
#define NUM_SCALAR_FIELDS_EACH_BLK (100 +13*NUM_CHARGED_SPECIES + 6*NUM_NEUTRAL_SPECIES+ 5*num_ion_prod_fields+4*num_externRhoVelocityFields +num_scalar_average_fields)

#elif (NONADIABATIC_GRADPE_TERM || USE_NEUTRAL_SPECIES_AS_FIELD) && !USE_DUST_SPECIES_AS_FIELD && !USE_ION_PRODUCTION_AS_FIELD

#define NUM_FIELDS 114
#define NUM_SCALAR_FIELDS_EACH_BLK (49 +13*NUM_CHARGED_SPECIES + NUM_PARTICLE_SPECIES+ 6*NUM_NEUTRAL_SPECIES+4*num_externRhoVelocityFields +num_scalar_average_fields)

#else

#define NUM_FIELDS 114
#define NUM_SCALAR_FIELDS_EACH_BLK (49 +13*NUM_CHARGED_SPECIES+ NUM_PARTICLE_SPECIES +4*num_externRhoVelocityFields +num_scalar_average_fields)

#endif


//! unfortunately the IBM compiler does not
//! support dynamic stack allocation of strings
//! so leave this value constant
#define INFO_ARRAY_SIZE 30

//! protocol particle that impact on the obstacle
constexpr bool USE_PROTOCOL_OBSTACLE_PARTICLE = false;
//! Comment below two defines out if above is false
// #define startTL_PROTOCOL_OBSTACLE_PARTICLE 3
// #define TL_PROTOCOL_ANY_PARTICLE 2


//! In case Neumann (outflow) Boundaries are used:
//! MARS  (fast plasma): set NM_BOUND_NODE to 1
//! DIONE (slow plasma): set NM_BOUND_NODE to 2
//! For some reason, Mars-Run got instable (NaN in B)
//! in case NM_BOUND_NODE was set to 2.
//! In case of slow plasma both cases do work, but in 
//! case of 1 strong wave reflection occured at the 
//! outflow boundary. This could be strongly decreased
//! by setting source to 2.
constexpr int NM_BOUND_NODE = 2;


constexpr int NUM_REQUESTS = 16;

constexpr int gathNEIB	= 6;


//! NOTE:
//! all numbers from 
//! id_rhoSpecies1 to id_rhoSpecies1 +5*NUM_CHARGED_SPECIES-1
//! are used for each species density, velocity and temperature !!!

//! E.g. in case of 2 ion species:
//! id 28: rho of species 1
//! id 29: rho of species 2
//! id 30: velocity of species 1
//! id 31: velocity of species 2
//! id 32: gyroradius of species 1
//! id 33: gyroradius of species 2
//! id 34: gyroradius of electrons species 1
//! id 35: gyroradius of electrons species 2
//! id 36: temperature of species 1
//! id 37: temperature of species 2

//!----------------------------------------
//!- DO NOT USE THESE NUMBERS ELSEWHERE!  -
//!----------------------------------------
#define id_rhoSpecies1		39
#define id_UI_Species1		(id_rhoSpecies1		+ NUM_CHARGED_SPECIES)
#define id_gyro_Species1	(id_UI_Species1		+ NUM_CHARGED_SPECIES)
#define id_gyro_el_Species1	(id_gyro_Species1		+ NUM_CHARGED_SPECIES)
#define id_recomb_Species1	(id_gyro_el_Species1	+ NUM_CHARGED_SPECIES)
#define id_ForceSpecies1	(id_recomb_Species1	+ NUM_CHARGED_SPECIES)
#define id_PISpecies1		(id_ForceSpecies1		+ NUM_CHARGED_SPECIES)
#define id_PESpecies1		(id_PISpecies1		+ NUM_CHARGED_SPECIES)
#define id_PEtotal			(id_PESpecies1		+ NUM_CHARGED_SPECIES)

//! neutral particles are used as parameters in
//! NONADIABATIC_GRADPE_TERM computations
//! NOTE: You may have to adjust NUM_FIELDS defined above
//!       if you use NONADIABATIC_GRADPE_TERM
// #define id_allRhoNeutralSpecies				(id_PEtotal + 1)
#define id_allRhoNeutralSpecies				    (id_PESpecies1 + NUM_CHARGED_SPECIES +1)
// #define id_numberdensity_neutralSpecies1	    (id_allRhoNeutralSpecies + 1)
#define id_numberdensity_neutralSpecies1		(id_PESpecies1 + NUM_CHARGED_SPECIES +2)
#define id_allUNeutralSpecies					(id_numberdensity_neutralSpecies1 + NUM_NEUTRAL_SPECIES)
#define id_velocity_neutralSpecies1			    (id_numberdensity_neutralSpecies1 + NUM_NEUTRAL_SPECIES + 1)
#define id_allPNeutralSpecies					(id_velocity_neutralSpecies1 + NUM_NEUTRAL_SPECIES)
#define id_pressure_neutralSpecies1			    (id_velocity_neutralSpecies1 + NUM_NEUTRAL_SPECIES + 1)
#define id_allnewBetaNeutralSpecies			    (id_pressure_neutralSpecies1 + NUM_NEUTRAL_SPECIES)
#define id_new_electron_beta_neutralSpecies1	(id_pressure_neutralSpecies1 + NUM_NEUTRAL_SPECIES + 1)

//! Use id_PESpecies1 as temporary scratch for PE_odd
#define id_PE_even	id_PEtotal
#define id_PE_odd		id_PESpecies1

#if (NONADIABATIC_GRADPE_TERM || USE_NEUTRAL_SPECIES_AS_FIELD)

    #define id_externRho1 (id_new_electron_beta_neutralSpecies1 + NUM_NEUTRAL_SPECIES)

#else

	#define id_externRho1 (id_PESpecies1 + NUM_CHARGED_SPECIES + 1)

#endif


#define id_extern_Ui1 (id_externRho1 + num_externRhoVelocityFields)


#if USE_DUST_SPECIES_AS_FIELD && USE_ION_PRODUCTION_AS_FIELD

	#define id_density_ionProdSpecies1	(id_extern_Ui1 + num_externRhoVelocityFields)
	#define id_velocity_ionProdSpecies1	(id_density_ionProdSpecies1 + num_ion_prod_fields)

	#define id_density_dustSpecies1       (id_velocity_ionProdSpecies1 + num_ion_prod_fields)
	#define id_velocity_dustSpecies1		(id_density_dustSpecies1 + num_Dust_Species)

	#define id_average_Field1 			(id_velocity_dustSpecies1 + num_Dust_Species)

#elif !USE_DUST_SPECIES_AS_FIELD && USE_ION_PRODUCTION_AS_FIELD

	#define id_density_ionProdSpecies1	(id_extern_Ui1 + num_externRhoVelocityFields)
	#define id_velocity_ionProdSpecies1	(id_density_ionProdSpecies1 + num_ion_prod_fields)

	#define id_density_dustSpecies1  -1
	#define id_velocity_dustSpecies1 -1

	#define id_average_Field1 			(id_velocity_ionProdSpecies1 + num_ion_prod_fields)

#elif USE_DUST_SPECIES_AS_FIELD && !USE_ION_PRODUCTION_AS_FIELD

	#define id_density_ionProdSpecies1 -1
	#define id_velocity_ionProdSpecies1 -1

	#define id_density_dustSpecies1       (id_extern_Ui1 + num_externRhoVelocityFields)
	#define id_velocity_dustSpecies1		(id_density_dustSpecies1 + num_Dust_Species)

	#define id_average_Field1 			(id_velocity_dustSpecies1 + num_Dust_Species)
#else

	#define id_density_ionProdSpecies1 -1
	#define id_velocity_ionProdSpecies1 -1

	#define id_density_dustSpecies1 -1
	#define id_velocity_dustSpecies1 -1

	#define id_average_Field1	(id_extern_Ui1  + num_externRhoVelocityFields)

#endif 



#define id_divrotB		(id_average_Field1 + num_average_fields +1)
#define id_divU			(id_divrotB + 1)
#define id_gradPE			(id_divU + 1)

#define id_BxgradB		(id_gradPE + 1)
#define id_RespProc		(id_BxgradB + 1)

//! refinement
#define id_Refine_Status	(id_RespProc + 1)

#define id_extern_PhiC	(id_Refine_Status + 1)

//! Note:
//! Dot not use numnbers higher 
//! than NUM_FIELDS defined above !!!
#define id_BOld		id_BOdd
#define id_BNew		id_BEven
#define id_KX			id_rotB
#define id_BMidStep	id_rotB
#define id_B_HI		id_rotB

#define id_rho_np1_HIMesh	id_PESpecies1
#define id_rho_n_HIMesh 	id_PISpecies1

//! define id for output of energy spectrum
//! (has to be larger than highest field id )
#define ENERGY_SPECTRUM (id_extern_PhiC + 1)