#include "parameters.h"
#include "absolute_Globals.h"
#include "CBlk.h"
#include "utils.h"

#include <math.h>



//! -------------------------------------------------------------------------------------------- !//
//! -------------------------------------------------------------------------------------------- !//
//! ------------------------------------FUNCTIONS/METHODS--------------------------------------- !//
//! -------------------------------------------------------------------------------------------- !//
//! -------------------------------------------------------------------------------------------- !//


//! -------------------------------------------------------------- !//
//! ------------------- Neutral Profile -------------------------- !//
//! -------------------------------------------------------------- !//
// Define the neutral profile to be used internally if one is not read from file(s)

// Ganymede neutral profile
void CBlock::neutral_profile(INT32 neutralSpecies, D_REAL *x, D_REAL *Nneutral, D_REAL *VNeutral)
{
  D_REAL z=vec_len(x); // (x,y,z) position in normalized units (i.e., position in meters divided by SI_x0)
  double R_Moon_km = R_Moon*SI_x0/1000; // un-normalize moon radius and divide by 1000 to get units of km
  double r_cut = 4.0*R_Moon; // maximum extent of atmosphere in normalized units

  // O2 at Ganymede
  if(neutralSpecies == 0) {
    if ((z > R_Moon) && (z < r_cut))
    {
      double exponent = (z/R_Moon - 1.0) * R_Moon_km; // (height in R_Moon / R_Moon - 1.0) * (moon radii in km)
      double nap0 = 6.0e+13; // surface density in 1/m^3
      double h1 = 250.0; // scale height in km
      Nneutral[0] = (nap0 * exp(-exponent/h1)) / SI_n0; // density in normalized units (norm'd to upstream density)
    }
    else {Nneutral[0] = 0.0;}      
  }

  // H2 at Ganymede
  else if(neutralSpecies == 1) {
    if ((z > R_Moon) && (z < r_cut))
    {
      double exponent = (z/R_Moon - 1.0)*R_Moon_km; // same as above
      double nap0 = 7.5e+12; // same as above
      double h1 = 1000.0; // same as above
      Nneutral[0] = (nap0 * exp(-exponent/h1)) / SI_n0; // same as above
    }
    else {Nneutral[0] = 0.0;}
  }

  // H2O at Ganymede
  else if (neutralSpecies == 2) {
    if ((z > R_Moon) && (z < r_cut))
    {
      // see ../data/asym_atmosphere_H2O.cpp for more detailed comments on the asymmetry
      // H2O is generally attenuated by a factor of cos()^gamma about the subsolar point

      double exponent = (z/R_Moon - 1.0) * R_Moon_km; // same as above
      double nap0 = 2.2e+14; // same as above
      double h1 = 200.0; // same as above
      double beta = 6.0; // attenuation exponent -> cos()^beta
     
      double sub_solar_angle = 132.0 * (M_PI/180); // see asym_atmosphere_H2O.cpp
      double x_ssp = cos(sub_solar_angle);
      double y_ssp = sin(sub_solar_angle);
      double z_ssp = 0.0;
      double r_ssp[3] = {R_Moon*x_ssp, R_Moon*y_ssp, R_Moon*z_ssp};

      double r_mag = vec_len(x);
      double r_ssp_mag = vec_len(r_ssp);
      double cos_alpha = (r_ssp[0]*x[0] + r_ssp[1]*x[1] + r_ssp[2]*x[2]) / (r_mag * r_ssp_mag);
      double asym_factor = pow(cos_alpha, beta);			  
      
      if (cos_alpha > 0.0) {
	      Nneutral[0] = (nap0 * asym_factor * exp(-exponent/h1)) / SI_n0;
      }
      else {
	      Nneutral[0] = 0.0;
      }
    }
    else {Nneutral[0] = 0.0;}
  }

  VNeutral[0]=0.;
  VNeutral[1]=0.;
  VNeutral[2]=0.;
}





//! -------------------------------------------------------------- !//
//! ------------------- Ionospheric Profile ---------------------- !//
//! -------------------------------------------------------------- !//
// Define the ionospheric production profile to be used internally if one is not read from file(s)
// This function is only used if set_analytical_ionProd_profile is True and
// ion_prod_profile_from_file is false.

// Enceladus
void CBlock::set_ion_production_profile(INT32 neutral_species)
{
	//! Saturn local time of all flybys
	D_REAL SLT[] = {22.6,17.0,17.1,23.2,22.7,22.5,22.5,10.9,3.7,9.3,3.9,3.6,8.8,8.8,23.8,23.8,23.8,0.6,0.6,0.5};
	//! Saturn's sub solar latitude of all flybys (in degrees)
	D_REAL SSL[] = {-22,-22,-21,-7.5,-5,-4.5,-4.25,1.275,1.565,3.954,4.132,5.543,7.132,7.435,11.39,11.62,12.05,13.8,13.92,14.28};

	D_REAL SLT_Degree, Sun_ENIS_angle, Sun_ENIS_vec[3];  
	  
	 SLT_Degree = SLT[flyby]/24.*360;
	 Sun_ENIS_angle = (270 - SLT_Degree)/360.*2.*M_PI;
	 
	 //! comment to SSL:
	 //! negative SSL is southern summer, that is a positive Sun_ENIS_vec[2]
	 //! thus pi + SSL with SSL<0 !
	 Sun_ENIS_vec[0] = cos(Sun_ENIS_angle)*sin(M_PI/2.+2.*M_PI*SSL[flyby]/360.);
	 Sun_ENIS_vec[1] = sin(Sun_ENIS_angle)*sin(M_PI/2.+2.*M_PI*SSL[flyby]/360.);
	 Sun_ENIS_vec[2] = cos(M_PI/2.+2.*M_PI*SSL[flyby]/360.);
	 
	bool Saturn_Shadow = false; 
	
	//! comment: length of geometric shadow in SLT:
	//! whole orbit is RS: 2 pi * 3.95
	//! divide by 24 to get 1 hour in RS or inverse to have 1 RS in hours
	//! ---> 1.934 local time hours
	//! no photionization from about 23 SLT to 1 SLT
	if(flyby==14 || flyby==15 || flyby==16)
	  Saturn_Shadow = true;

	
	INT32  ind[3];
  INT32 i_j_k;

	D_REAL *IonProd, *UnX, *UnY, *UnZ;
  D_REAL dist_to_shadow_axis, abs2_Sun_ENIS_vec, x_scalar_Sun_ENIS_vec;

	D_REAL x[3], Nneutral[1], Vneutral[3], cell_intern_r[3];
	D_REAL temp[3];

	memset(Nneutral,0,sizeof(D_REAL));
	memset(Vneutral,0,3*sizeof(D_REAL));
	memset(cell_intern_r,0,3*sizeof(D_REAL));

	//! Set pointer to intern (Block) field
	IonProd= Field_Type[id_density_ionProdSpecies1 + neutral_species];
	UnX = Field_Type[id_velocity_neutralSpecies1 + neutral_species];
	UnY = UnX +num_nodes_in_block;
	UnZ = UnY +num_nodes_in_block;

	
  for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
  {
    for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
    {
      for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
      {
        i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
          +ind[1]*BlkNds_Z
          +ind[2];

        if(Saturn_Shadow)
        {
          IonProd[i_j_k] = 0;
          continue;
        }
        
        intern2normedCoords(x,cell_intern_r,ind);
          
        x_scalar_Sun_ENIS_vec = x[0]*Sun_ENIS_vec[0]
                + x[1]*Sun_ENIS_vec[1]
                + x[2]*Sun_ENIS_vec[2];
        abs2_Sun_ENIS_vec = vec_len2(Sun_ENIS_vec);
        
        temp[0] = x[0]-
          (Sun_ENIS_vec[0]*x_scalar_Sun_ENIS_vec)/abs2_Sun_ENIS_vec;
                  
                  
        temp[1] = x[1]-
          (Sun_ENIS_vec[1]*x_scalar_Sun_ENIS_vec)/abs2_Sun_ENIS_vec;

            
        temp[2] = x[2]-
          (Sun_ENIS_vec[2]*x_scalar_Sun_ENIS_vec)/abs2_Sun_ENIS_vec;
      
        
        dist_to_shadow_axis = vec_len(temp);

        
        if(dist_to_shadow_axis < R_Obstacle && x_scalar_Sun_ENIS_vec >0)
        IonProd[i_j_k] = 0;

      }
    }
  }
}


//! -------------------------------------------------------------- !//
//! ------------------- Recombination Rates ---------------------- !//
//! -------------------------------------------------------------- !//
//! Calculate the recombination rates of the ionospheric species.
//! This is only used if any of the boolean switches in
//! Recombination_for_Species[] are set to true/1 in parameters.h

//!-------------------------------------------------------------//
//! Use Eq. 11 after Gombosi 1996
//!-------------------------------------------------------------//
void CBlock::calc_RecombinationAlphaField(void)
{	
	D_REAL *RecombinationRate;
	
	for(INT32 species=0; species<num_Particle_Species; species++)
	if(Recombination_for_Species[species])	
	{	
		RecombinationRate = Field_Type[id_recomb_Species1 +species];
		
		for (INT32 node=0; node < num_nodes_in_block; node++)
		{
// 			if(Te[node]<=200)
// 			{
// 			RecombinationRate[node] = sqrt(300/Te[node]);
// 			}
// 			else 
// 			{
// 			D_REAL exp=0.2553-0.1633*log(Te[node]);
// 			RecombinationRate[node] = 2.342*pow(Te[node],exp);
// 			}
		        //if(species==0)RecombinationRate[node]=3.5e-12*sqrt(300/Te[node])*SI_t0*1e-6*SI_n0;
			//O2
		  //		        if(species==1)RecombinationRate[node]=2.4e-7*pow(300./Te_ionosphere,0.7)*SI_t0*1e-6*SI_n0;  //! from NS
		  //	if(species==2)
		  //	  if(Te_ionosphere<800){
		  //	    RecombinationRate[node]=1.57e-5*pow(Te_ionosphere,-0.569)*SI_t0*1e-6*SI_n0;
		  //	  }else if(Te_ionosphere<4000){
		  //	      RecombinationRate[node]=4.73e-5*pow(Te_ionosphere,-0.74)*SI_t0*1e-6*SI_n0;
		  //	    }else RecombinationRate[node]=1.03e-3*pow(Te_ionosphere,-1.111)*SI_t0*1e-6*SI_n0;
			//not needed, only recombination for species 1 and 2
			//if(species==3)RecombinationRate[node]=7.e-7*sqrt(300/Te[node])*SI_t0*1e-6*SI_n0;
			//if(species==4)RecombinationRate[node]=1.9e-6*sqrt(300/Te[node])*SI_t0*1e-6*SI_n0;
			//if(species==5)RecombinationRate[node]=4.8e-7*sqrt(300/Te[node])*SI_t0*1e-6*SI_n0; //! krasnopolsky
			//if(species==6)RecombinationRate[node]=1.e-6*sqrt(300/Te[node])*SI_t0*1e-6*SI_n0;
			//if(species==7)RecombinationRate[node]=1e-6*sqrt(300/Te[node])*SI_t0*1e-6*SI_n0;
			//if(species==8)RecombinationRate[node]=1e-6*sqrt(300/Te[node])*SI_t0*1e-6*SI_n0;
			
		}
	}
}
