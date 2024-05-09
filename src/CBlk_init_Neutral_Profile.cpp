
#include "CBlk.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>

using namespace std;

extern D_REAL **delta_of_L;
extern D_REAL *CellVol_of_L;
extern D_REAL **Blk_Length_of;
extern INT32 i_j_k;


//!-------------------------------------------------------------//
//! init_Neutral_Profile: 					//
//!-------------------------------------------------------------//

//! set neutral profile simply by using a new function in parameters.cpp
//! neutral_profile(NeutralSpecies, x[3], Nneutral, VNeutral[3])

//!-------------------------------------------------------------//
//! for neutral density					//
//!-------------------------------------------------------------//
void CBlock::set_analytical_neutral_profile(INT32 neutralSpecies)
{
	INT32  ind[3];

	D_REAL *NeutralRho, *UnX, *UnY, *UnZ, *p, *beta_new;

	D_REAL x[3], Nneutral[1], Vneutral[3], cell_intern_r[3];
	
	memset(Nneutral,0,sizeof(D_REAL));
	memset(Vneutral,0,3*sizeof(D_REAL));
	memset(cell_intern_r,0,3*sizeof(D_REAL));

	//! Set pointer to intern (Block) field
	NeutralRho= Field_Type[id_numberdensity_neutralSpecies1 + neutralSpecies];
	UnX = Field_Type[id_velocity_neutralSpecies1 + neutralSpecies];
	UnY = UnX +num_nodes_in_block;
	UnZ = UnY +num_nodes_in_block;
	p  = Field_Type[id_pressure_neutralSpecies1 + neutralSpecies];
	beta_new = Field_Type[id_new_electron_beta_neutralSpecies1 + neutralSpecies];
	

	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
	
	
		i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
			
		intern2normedCoords(x,cell_intern_r,ind);
		
		neutral_profile(neutralSpecies,x,Nneutral,Vneutral);
			
		NeutralRho[i_j_k] = Nneutral[0];	
		
		UnX[i_j_k] = Vneutral[0];
		UnY[i_j_k] = Vneutral[1];
		UnZ[i_j_k] = Vneutral[2];
		
		
		//!pressure is plasma beta at background density * normed density
		p[i_j_k]  = NeutralRho[i_j_k]*Neutral_Betas[neutralSpecies]; 
                                
		beta_new[i_j_k] = NewElectron_Betas[neutralSpecies];
	  }
	  
}


//!-------------------------------------------------------------//
//! set_RhoUi_extern: 								//
//!-------------------------------------------------------------//
void CBlock::set_RhoUi_extern(INT32 id_densityfield, INT32 id_velocityfield, INT32* num_Nodes, D_REAL* Origin, D_REAL* Length, D_REAL* rho, INT32& num_values_not_in_extern_box)
{



	INT32  ind[3];
	INT32  num_extern_box_nodes;
	INT32 a,b,c, a_b_c;
	INT32 ap1_b_c, a_bp1_c, a_b_cp1;
	INT32 ap1_bp1_c, ap1_b_cp1, a_bp1_cp1, ap1_bp1_cp1;

	D_REAL *RHO, *UiX, *UiY, *UiZ;
	D_REAL *extern_UiX, *extern_UiY, *extern_UiZ;
	D_REAL r[3], extern_delta[3], shape_func[8];

	PARTICLE_REAL x_BlockNode[3], x_extern[3], cell_intern_r[3];
	memset(cell_intern_r,0,3*sizeof(PARTICLE_REAL));


	for(INT32 comp=0; comp<3; comp++)
	extern_delta[comp] = Length[comp]/(num_Nodes[comp]-1);


	//! Set pointer to extern field
	num_extern_box_nodes =   num_Nodes[0]
				*num_Nodes[1]
				*num_Nodes[2];

	extern_UiX = rho +num_extern_box_nodes;
	extern_UiY = extern_UiX +num_extern_box_nodes;
	extern_UiZ = extern_UiY +num_extern_box_nodes;
	
	
	//! Set pointer to intern (Block) field
	RHO = Field_Type[id_densityfield];

	UiX = Field_Type[id_velocityfield];
	UiY = UiX +num_nodes_in_block;
	UiZ = UiY +num_nodes_in_block;


	//! use a,b,c for extern mesh 
	//! use i,j,k for intern mesh
	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
	
	
		i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
		//! get Coordinate of intern Block Node 
		intern2normedCoords(x_BlockNode, cell_intern_r, ind);
		
		//! find lower next neighbour node in read Mesh
		x_extern[0] = (x_BlockNode[0] +Origin[0])/extern_delta[0];
		x_extern[1] = (x_BlockNode[1] +Origin[1])/extern_delta[1];
		x_extern[2] = (x_BlockNode[2] +Origin[2])/extern_delta[2];

		 a  = int(x_extern[0]);
		 b  = int(x_extern[1]);
		 c  = int(x_extern[2]);


		if(   (x_extern[0]<0) || a>num_Nodes[0]-2
		    ||  x_extern[1]<0 || b>num_Nodes[1]-2
		    ||  x_extern[2]<0 || c>num_Nodes[2]-2)
		  {

			num_values_not_in_extern_box++;


		  }
		  else
		  {
		  
			r[0] = x_extern[0]-a;
			r[1] = x_extern[1]-b;
			r[2] = x_extern[2]-c;
			
			shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2]);
			shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2]);
			
			shape_func[4] = (   r[0])*(   r[1])*(1.-r[2]);
			shape_func[5] = (   r[0])*(1.-r[1])*(   r[2]);
			shape_func[6] = (1.-r[0])*(   r[1])*(   r[2]);
			shape_func[7] = (   r[0])*(   r[1])*(   r[2]);
			
			//! -----------------------------------------------
			a_b_c   =      a*num_Nodes[1]*num_Nodes[2]    +b*num_Nodes[2]      +c;
			
			//! -----------------------------------------------
			ap1_b_c = (a+1)*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2]     +c;
			a_bp1_c =     a*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2]     +c;
			a_b_cp1 =     a*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2] +(c+1);
			
			//! ------------------------------------------------
			ap1_bp1_c = (a+1)*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2]     +c;
			ap1_b_cp1 = (a+1)*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2] +(c+1);
			a_bp1_cp1 =     a*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2] +(c+1);
			
			ap1_bp1_cp1 = (a+1)*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2] +(c+1);

			RHO[i_j_k]  =  rho[  a_b_c] * shape_func[0]
				      +rho[ap1_b_c] * shape_func[1]
				      +rho[a_bp1_c] * shape_func[2]
				      +rho[a_b_cp1] * shape_func[3]
			
				      +rho[  ap1_bp1_c] * shape_func[4]
				      +rho[  ap1_b_cp1] * shape_func[5]
				      +rho[  a_bp1_cp1] * shape_func[6]
				      +rho[ap1_bp1_cp1] * shape_func[7];

			//! neutral density must not be smaller than zero	      
			if(RHO[i_j_k]  < 0.  )
			{
	
				log_file << endl << "ERROR:" << endl
				<< " Negative Values not allowed." << endl;

				log_file <<  "ind[0]: " << ind[0] << endl;
				log_file <<  "ind[1]: " << ind[1] << endl;
				log_file <<  "ind[2]: " << ind[2] << endl;

				log_file <<  " x_extern[0]: " << x_extern[0] << endl;
				log_file <<  " x_extern[1]: " << x_extern[1] << endl;
				log_file <<  " x_extern[2]: " << x_extern[2] << endl;

				log_file <<  " a: " << a << endl;
				log_file <<  " b: " << b << endl;
				log_file <<  " c: " << c << endl;

				for(int f=0; f<8; f++)
				log_file <<  "shape_func["<<f<<"]: " << shape_func[f] << endl;


				log_file <<  "RHO[i_j_k]: " << RHO[i_j_k]<< endl;
				
// 				finalize_MPI();
				exit(1);
			

			}

			UiX[i_j_k]  =  extern_UiX[  a_b_c] * shape_func[0]
				      +extern_UiX[ap1_b_c] * shape_func[1]
				      +extern_UiX[a_bp1_c] * shape_func[2]
				      +extern_UiX[a_b_cp1] * shape_func[3]
			
				      +extern_UiX[  ap1_bp1_c] * shape_func[4]
				      +extern_UiX[  ap1_b_cp1] * shape_func[5]
				      +extern_UiX[  a_bp1_cp1] * shape_func[6]
				      +extern_UiX[ap1_bp1_cp1] * shape_func[7];

			UiY[i_j_k]  =  extern_UiY[  a_b_c] * shape_func[0]
				      +extern_UiY[ap1_b_c] * shape_func[1]
				      +extern_UiY[a_bp1_c] * shape_func[2]
				      +extern_UiY[a_b_cp1] * shape_func[3]
			
				      +extern_UiY[  ap1_bp1_c] * shape_func[4]
				      +extern_UiY[  ap1_b_cp1] * shape_func[5]
				      +extern_UiY[  a_bp1_cp1] * shape_func[6]
				      +extern_UiY[ap1_bp1_cp1] * shape_func[7];

			UiZ[i_j_k]  =  extern_UiZ[  a_b_c] * shape_func[0]
				      +extern_UiZ[ap1_b_c] * shape_func[1]
				      +extern_UiZ[a_bp1_c] * shape_func[2]
				      +extern_UiZ[a_b_cp1] * shape_func[3]
			
				      +extern_UiZ[  ap1_bp1_c] * shape_func[4]
				      +extern_UiZ[  ap1_b_cp1] * shape_func[5]
				      +extern_UiZ[  a_bp1_cp1] * shape_func[6]
				      +extern_UiZ[ap1_bp1_cp1] * shape_func[7];
			
	      
 
			if(Flag[i_j_k])
			{
				RHO[i_j_k] = 0;
				UiX[i_j_k] = 0;
				UiY[i_j_k] = 0;
				UiZ[i_j_k] = 0;
			}	

		}
	  }

}


// TODO: Separate testing suite

void CBlock::init_neutral_profile_test(INT32 neutralSpecies, D_REAL* x, D_REAL* Nneutral, D_REAL* VNeutral)
{
  //! Set density
  Nneutral[0]  = 2.00e8;
  
  //! Set velocity
  VNeutral[0] = V_sw[0];
  VNeutral[1] = V_sw[1];
  VNeutral[2] = V_sw[2];
}
// nrho[u_v_w]     = 2.00e9;
// ux[u_v_w]       = V_sw[0];
// uy[u_v_w]       = V_sw[1];
// uz[u_v_w]       = V_sw[2];
// p[u_v_w]        = 2.776e6;
// beta_new[u_v_w] = 1.38e-3;

