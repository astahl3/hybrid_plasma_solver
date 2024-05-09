

#include "CBlk.h"
#include "parameters.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <math.h>


using namespace std;




//!-------------------------------------------------------------//
//! init_Block: 									//
//!-------------------------------------------------------------//
void CBlock::init_Eta_Profile(void)
{
	// Create Eta profile according to variables and switches set in parameters.h
	
	INT32 u_v_w;
	INT32 cell_indices[3];

	D_REAL r_vec[3];
	D_REAL null_vec[3] = {0.,0.,0.};

	D_REAL dist_to_orig, eta_join, eta;
	D_REAL xhalfaxis, zhalfaxis, superellipseexponent;     
	D_REAL asym_of_box;

	D_REAL* Eta = Field_Type[id_Eta];

	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		u_v_w  = u*BlkNds_Y*BlkNds_Z 
			      +v*BlkNds_Z 
			      +w;
		

		cell_indices[0] = u;
		cell_indices[1] = v;
		cell_indices[2] = w;

		if(mesh_type == STAGGERED)
		intern2normedCoords_HIMesh(r_vec, null_vec, cell_indices);
		else
		intern2normedCoords(r_vec, null_vec, cell_indices);


		dist_to_orig = vec_len(r_vec);

		eta_join = Eta_Obstacle * 1./(exp((bridge_profile_end* R_Obstacle -R_Eta_Obstacle)*fermi_slope)+1.) + Eta_sw;


		//! set Eta_sw as background-eta value
                
		Eta[u_v_w] = Eta_sw;

		//! better use fermi function than define L0 & smooth:
		//! In case only defined on L0, Eta has to be interpolated
		//! to higher levels which lead to "rectangular shape".
	
		//! NOTE:
		//! In case "smooth ETA" version is used, do not forget to
		//! include field_from_parent(Eta) function in refine_Block !
		if (use_eta_fermi)
			Eta[u_v_w]= Eta_Obstacle *1./(exp((dist_to_orig -R_Eta_Obstacle)*fermi_slope)+1.) + Eta_sw;
		else if (!use_eta_fermi && dist_to_orig < R_Eta_Obstacle)
			Eta[u_v_w] = Eta_Obstacle;
		else if (use_eta_comet)
			Eta[u_v_w]= Eta_comet_inf *(1.-1/(exp((dist_to_orig -R_Eta_innerComa)*fermi_slope_comet)+1.)) + Eta_comet_innerComa;

		if (use_intermediate_core && dist_to_orig < bridge_profile_end*R_Obstacle)
		    Eta[u_v_w]= eta_join/((bridge_profile_end -bridge_profile_start)*R_Obstacle)*(dist_to_orig -bridge_profile_start*R_Obstacle);
	
		
                
		//! eta should be zero inside obstacle core
		if (use_intermediate_core && dist_to_orig < bridge_profile_start* R_Obstacle)
			Eta[u_v_w] = 0.;
		else if (!use_intermediate_core && dist_to_orig < obstacle_core_fraction* R_Obstacle)
		    Eta[u_v_w] = 0.;

	}

	//! set Eta boundaries
	if(eta_Alfven_Wing_boundary)
	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		u_v_w  = u*BlkNds_Y*BlkNds_Z 
			     +v*BlkNds_Z 
			     +w;
		


		cell_indices[0] = u;
		cell_indices[1] = v;
		cell_indices[2] = w;
		intern2normedCoords(r_vec, null_vec, cell_indices);
		
// 		//! upper boundary
// 		eta = resistive_bound_eta[4]/resistive_bound_dist[4] *( r_vec[2] - ( LZ - Box_Origin[2] - resistive_bound_dist[4]));
// 		
// 		if(eta>Eta_sw)
// 		Eta[u_v_w] = eta;
// 		
// 		//! lower boundary
// 		eta = resistive_bound_eta[5]/resistive_bound_dist[5] *( -r_vec[2] + ( -Box_Origin[2] + resistive_bound_dist[5]) );
// 		

		//! below are now in parameters.h
		// xhalfaxis = 6*R_Moon;  
		// zhalfaxis = 7*R_Moon;  
		// superellipseexponent = 6;  
		
		asym_of_box = Box_Origin[2]-0.5*LZ;
		
		eta = resistive_bound_eta[4]*(pow((r_vec[0]+ Box_Origin[0])/eta_Alfven_xhalfaxis,superellipseexponent) + pow((r_vec[2]+ asym_of_box)/eta_Alfven_zhalfaxis,superellipseexponent) - 1)
			/(pow((LX-Box_Origin[0]+ Box_Origin[0])/eta_Alfven_xhalfaxis,superellipseexponent) + pow((LZ-Box_Origin[2]+ asym_of_box)/eta_Alfven_zhalfaxis,superellipseexponent) - 1);
	
		if(eta>Eta_sw && vec_len(r_vec)>5.*R_Moon)
		Eta[u_v_w] = eta;
	
	    }


	//! set Eta boundaries
	//!----------- -X boundary ------------------------
	if(is_box_boundary[_Im1_] && resistive_bound_cells[_Im1_])
	 for(INT32 u=0; u<resistive_bound_cells[_Im1_]; u++)
	  for(INT32 v=0; v<BlkNds_Y; v++)
	   for(INT32 w=0; w<BlkNds_Z; w++)
	   {
		u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Im1_];

	    }



	//!----------- +X boundary ------------------------
	if(is_box_boundary[_Ip1_] && resistive_bound_cells[_Ip1_])
	 for(INT32 u=BlkNds_X-resistive_bound_cells[_Ip1_]; u<BlkNds_X; u++)
	  for(INT32 v=0;					     v<BlkNds_Y; v++)
	   for(INT32 w=0;				      w<BlkNds_Z; w++)
	   {
		u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Ip1_];

	   }


	//!----------- -Y boundary ------------------------
	if(is_box_boundary[_Jm1_]  && resistive_bound_cells[_Jm1_])
	 for(INT32 u=0; u<BlkNds_X; u++)
	  for(INT32 v=0; v<resistive_bound_cells[_Jm1_]; v++)
	   for(INT32 w=0; w<BlkNds_Z; w++)
	   {
		u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Jm1_];

	    }

	//!----------- +Y boundary ------------------------
	if(is_box_boundary[_Jp1_]  && resistive_bound_cells[_Jp1_])
	 for(INT32 u=0;					   u<BlkNds_X; u++)
	  for(INT32 v=BlkNds_Y-resistive_bound_cells[_Jp1_]; v<BlkNds_Y; v++)
	   for(INT32 w=0;				      w<BlkNds_Z; w++)
	   {

		u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Jp1_];

	   }

	//!----------- -Z boundary ------------------------
	if(is_box_boundary[_Km1_]  && resistive_bound_cells[_Km1_])
	 for(INT32 u=0; u<BlkNds_X; u++)
	  for(INT32 v=0; v<BlkNds_Y; v++)
	   for(INT32 w=0; w<resistive_bound_cells[_Km1_]; w++)
	   {

		u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Km1_];

	  }


	//!----------- +Z boundary ------------------------
	if(is_box_boundary[_Kp1_]  && resistive_bound_cells[_Kp1_])
	 for(INT32 u=0;					    u<BlkNds_X; u++)
	  for(INT32 v=0;					     v<BlkNds_Y; v++)
	   for(INT32 w=BlkNds_Z -resistive_bound_cells[_Kp1_]; w<BlkNds_Z; w++)
	   {

		u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Kp1_];

	   }





}
