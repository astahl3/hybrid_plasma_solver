#include <cmath>
#include <iostream>
#include <ctime>
#include <fstream>
#include <iostream>

#define R_C 1.
#define dt 0.014
#define boxdim 8. //Don't forget the '.'

const int grid=600; //For grid>300 when compiling, use 'g++ -mcmodel=large ...'
const int num_pts_mesh=grid*grid*grid;

double Hrel = 10;

//! 0:            test/defaut case
//! flyby number: see below
int flyby = 0;


double LT;
double SLT;
double SSL; 

double asym;

double origin[3]={0.5001*boxdim,0.5001*boxdim,0.5001*boxdim};
double x[grid],y[grid],z[grid];
double pos_array[3*num_pts_mesh];
double suv[grid][grid][grid];
double suv_array[num_pts_mesh];
double r;
double vec_x[3];
double sun[3];

//! NOTE: make sure to redifine sun vector based on SLT in switch below
//! (should already be included in the switch statement below)

double sigma_CO2_abs[37]={1.550,4.616,9.089,14.361,16.505,19.016,17.518,21.492,21.594,23.574,25.269,24.871,28.271,29.526,30.254,31.491,33.202,34.200,34.913,35.303,34.300,34.447,33.699,23.518,32.832,93.839,61.939,26.493,39.831,13.980,44.673,52.081,42.869,50.311,15.100,14.200,18.241};

double sigma_CO2_ion[37]={0.447,2.083,4.960,8.515,11.113,13.004,11.906,14.390,14.414,15.954,18.271,17.982,21.082,24.378,27.163,30.138,31.451,32.382,33.482,34.318,33.795,34.003,32.287,20.856,27.490,86.317,51.765,21.676,34.094,10.930,7.135,0.,0.,0.,0.,0.,0.};


double sigma_O2_abs[37]={1.316,3.806,7.509,10.900,13.370,15.790,14.387,16.800,16.810,17.438,18.320,18.118,20.310,21.910,23.101,24.606,26.040,22.720,26.610,28.070,32.060,26.017,21.919,27.440,28.535,20.800,18.910,26.668,22.145,16.631,8.562,12.817,18.730,21.108,1.630,1.050,1.346};

double sigma_O2_ion[37]={1.316,2.346,4.139,6.619,8.460,9.890,9.056,10.860,10.880,12.229,13.760,13.418,15.490,16.970,17.754,19.469,21.600,18.840,22.789,24.540,30.070,23.974,21.116,23.750,23.805,11.720,8.470,10.191,10.597,6.413,5.494,9.374,15.540,13.940,1.050,0.,0.259};



double dlambda[37]={50, 50, 50, 50, 1, 1, 50, 1, 1, 50, 1, 50, 50, 1, 50, 50, 1, 1, 50,1, 1, 50, 50, 1, 50, 1, 1, 1, 50, 50, 50, 50, 1, 50, 1, 1, 50};


//! I_Inf_i * 10^13 for SI
//
double I_inf_i[37]={1.2,0.45,4.8,3.1,0.46,0.21,1.679,0.8,6.9,0.965,0.650,0.314,0.383,0.29,0.285,0.452,0.72,1.27,0.357,0.53,1.59,0.342,0.32,0.36,0.141,0.17,0.26,0.702,0.758,1.625,3.537,3,4.4,1.475,3.5,2.1,2.467};
//
double A[37]={1.0017e-2,7.1250e-3,1.3375e-2,1.9450e-2,2.7750e-3,1.3768e-1,2.6467e-2,2.5000e-2,3.3333e-3,2.2450e-2,6.5917e-3,3.6542e-2,7.4082e-3,7.4917e-3,2.0225e-2,8.7583e-3,3.2667e-3,5.1583e-3,3.6583e-3,1.6175e-2,3.3250e-3,1.18e-2,4.2667e-3,3.0417e-3,4.7500e-3,3.8500e-3,1.2808e-2,3.2750e-3,4.7667e-3,4.8167e-3,5.6750e-3,4.9833e-3,3.9417e-3,4.4168e-3,5.1833e-3,5.2833e-3,4.3750e-3};

double F107, F107A, flyby_dist_sun_jupiter;


double R_cut=5.; //! R_c is not included here
double max_rate=0;
double pos_max[3]={0.,0.,0.};


int count_neg_value=0;
int count_beneath_surface=0;
int count_shadow=0;
int count_point_ionosphere=0;
int count_progress=0;

using namespace std;

//double cos_chi(double *a, double *b)
//{
// double hlp=(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/sqrt((a[0]*a[0]+a[1]*a[1]+a[2]*a[2])*(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]));
// if(hlp>0.05)return hlp;
//  else return 0.05;
//}

void calc_solar_flux_flyby()
{
  cout << "Flyby number: " << flyby << endl;
  switch(flyby)
    {//! Quick standard secnario
      case 0: 
	F107=150.; //! daily flux
	F107A=150.; //! averaged flux +-1 month
	flyby_dist_sun_jupiter=5. ;
	LT = 18.;
	SSL = M_PI/2;
	SLT = -M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
        sun[1]=sin(SSL)*sin(SLT);
        sun[2]=cos(SSL);
	break;
      //! C3 on date
      case 3: 
	F107=67.8; //! daily flux
	F107A=72.6; //! averaged flux +-1 month
	flyby_dist_sun_jupiter=5.15;  //!AE
	LT=7.33;
	SSL=-(-1.145*M_PI)/180.+M_PI/2.;
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);
	break;
      //! C9 on date
      case 9:
	F107=74.1; //! daily flux
	F107A=74.7; //! averaged flux +-1 month
	flyby_dist_sun_jupiter=5.07;  //!AE
	LT=5.53;
	SSL=-(-0.147*M_PI)/180.+M_PI/2.;
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);
	break;
	//! C10 on date
      case 10: 
	F107=94.0; //! daily flux
	F107A=88.4; //! averaged flux +-1 month  NOTE this is approximated
	flyby_dist_sun_jupiter=5.05;  //!AE
	LT=5.03;
	SSL=-(0.362*M_PI)/180.+M_PI/2.; //SSL in deg *pi/180 + pi/2 --> SSL in rad (needed)
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);	
	break;
    }
    
  double P=(F107+F107A)/2.;
  for(int i=0;i<37;i++)
  {
    I_inf_i[i]*=(1+A[i]*(P-80))/(flyby_dist_sun_jupiter*flyby_dist_sun_jupiter);   
  }
  
  
  
}

//! Neutral CO2 profile
double CO2(double z, bool &in_atmosphere, double *z_vec)
{
  asym = sqrt(0.5*(1-(z_vec[0])/(sqrt(z_vec[0]*z_vec[0]+z_vec[1]*z_vec[1]+z_vec[2]*z_vec[2]))));

  if((abs(z_vec[0])<=R_cut && abs(z_vec[1])<=R_cut && abs(z_vec[2])<=R_cut))
    //Comment out asym term for a symmetric atmosphere!
    return asym/*asym_sun*/*(4.e+14/Hrel)*exp(-2410*(z-R_C)/(23*Hrel)); //! actual scale height is 23 km, 4e14 in m^-3
  else
  {
    in_atmosphere=false;
    return asym/*asym_sun*/*(4.e+14/Hrel)*exp(-2410*(z-R_C)/(23*Hrel));
  }
}


//! Neutral H2 profile
double H2(double z, bool &in_atmosphere, double *z_vec)
{
  double R_Moon=R_C;
  double nap0=3.75e+12;
  double h0=1000.0;
  double nL,n;  
  double exponent;
  double rmag;
  double rmax = 4.0;

  // only populate neutral atmosphere to certain distance from moon
  rmag = sqrt(z_vec[0]*z_vec[0]+z_vec[1]*z_vec[1]+z_vec[2]*z_vec[2]);
  if (rmag <= rmax) {
    exponent=(z-R_Moon)*2634.1;
    nL=nap0*exp(-exponent/h0);
  }
  else {
    nL=0.0; 
  }
  n=nL;
  return (n);
}

double calc_photo_rate_trapezrule(double *z_vec, double z0, int i, int j, int k, double Hrel, double RC)
{
  double not_in_shadow=1.;
  long double Prod=0.;
  bool do_integrate=true;
  double H2x0=H2(z0,do_integrate,z_vec);
  double temp=0.;
  double temp1=0.;
  double temp2=0.;
  double a=z0;
  double theta;
  double theta_SS = 132.0; // sub-solar point in xy plane in degrees CCW from x=0
  double theta_LB,theta_UB; // LB = theta - 90 and UB = theta + 90
  double e_rate = 6.0e-8; // rate in 1/s
  double p_rate = 3.13e-9; // rate in 1/s
  double p_rate_temp = 0.0;

  z_vec[0]+=dt*sun[0];
  z_vec[1]+=dt*sun[1];
  z_vec[2]+=dt*sun[2];

  double b=sqrt(z_vec[0]*z_vec[0]+z_vec[1]*z_vec[1]+z_vec[2]*z_vec[2]);
  
  theta = atan2(z_vec[1],z_vec[0]) * 180.0 / M_PI;
  if (theta < 0) {
    theta = theta + 360.0; // atan2 gives theta in [-pi,pi], so add 360 to negative to get [0, 2pi]
  }  

  theta_LB = theta_SS - 90.0;
  theta_UB = theta_SS + 90.0;
  
  if (theta_LB < 0) {
    theta_LB = theta_LB + 360;
  }
  theta_UB = fmod(theta_UB,360); 

  if ((theta >= theta_LB) && (theta <= theta_UB)) {
    p_rate_temp = p_rate; 
  }

  else { 
    p_rate_temp = 0.0; 
  }

  //  return e_rate*H2x0 + p_rate_temp*H2x0;
  return e_rate*H2x0;
}


int main()
{
  clock_t start,finish;
  double time;
  start = clock();
  double zvec[3];
  //!Fill grid with coordinates, Box with 20^3 R_C, Origin in center, This is the global mesh
  for(int i=0;i<grid;i++)
  {
    x[i]=-origin[0]+i*boxdim/(grid-1);
    y[i]=-origin[1]+i*boxdim/(grid-1);
    z[i]=-origin[2]+i*boxdim/(grid-1);
  }
  int i_j_k;
  
  //calc_solar_flux_flyby();
  int percentage=0;
  //! Calculate SUV rate at every grid point
  for(int i=0;i<grid;i++){
    if (((int((double(i)/grid)*100.)%10)==0)&&(percentage!=int(double(i)/grid*100.))){
      printf("%d Percent\n",int(double(i)/grid*100.));
      percentage=int(double(i)/grid*100.);
    }
    for(int j=0;j<grid;j++)
      for(int k=0;k<grid;k++)
      {
	count_progress++;

	zvec[0]=x[i];zvec[1]=y[j];zvec[2]=z[k]; //! this vector has global cartesic coordinates
	r=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);  //! R_C not included meaning this is NOT altitude above ground


	if(r>=R_C)//! This guarantees that our starting point is at least on the surface
	{
	  i_j_k=i*grid*grid+j*grid+k;	  
	  //! SUV-rate function call
	  suv[i][j][k]=calc_photo_rate_trapezrule(zvec,r,i,j,k,Hrel,R_C); 
// 	  pos_array[0*num_pts_mesh+i_j_k]=x[i];
// 	  pos_array[1*num_pts_mesh+i_j_k]=y[j];
// 	  pos_array[2*num_pts_mesh+i_j_k]=z[k];
	  if(suv[i][j][k]<0)count_neg_value++; //!negative values are wrong, check the code
	  suv_array[i_j_k]=suv[i][j][k];
	  //!Global Maxium determination
	  if(suv[i][j][k]>max_rate)
	  {
	    max_rate=suv[i][j][k];
	    pos_max[0]=x[i];
	    pos_max[1]=y[j];
	    pos_max[2]=z[k];
	  }

	}
	else //! obviously we are in the planet with our point here, no ionization rate here
	{
	  count_beneath_surface++;
	  suv[i][j][k]=0.; //! Value of 1. for better Visualization
	}
		
	//if(count_progress%(num_pts_mesh/100)==0)
	//  cout<<"\rDone with "<<100*count_progress/num_pts_mesh<<" % of total mesh points"<<flush;

      }
  }
  count_point_ionosphere=grid*grid*grid-count_beneath_surface;
  //cout<<"\nNumber points beneath Callisto surface is: "<<count_beneath_surface<<endl;
  //cout<<"Number points in Callisto shadow is: "<<count_shadow<<endl;
  //cout<<"Resolution: "<<boxdim*1560.8/grid<<"km"<<endl;
  //cout<<"Max value at grid: "<<max_rate<<endl;
  //cout<<"Position of max value: R: "<<sqrt(pos_max[0]*pos_max[0]+pos_max[1]*pos_max[1]+pos_max[2]*pos_max[2])<<", theta: "<<acos(pos_max[2]/sqrt(pos_max[0]*pos_max[0]+pos_max[1]*pos_max[1]+pos_max[2]*pos_max[2]))<<", phi: "<<atan(pos_max[1]/pos_max[0])<<endl;
  //cout<<"neg values : "<<count_neg_value<<endl;
 
//!Output stuff, if grid is changed, output must also be changed for same data
  //!Output z=0 ebene 
//   ofstream out;
//   out.open("N2-3d.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//       for(int k=0;k<grid;k++)
// 	out<<x[i]<<" "<<y[j]<<" "<<z[k]<<" "<<suv[i][j][k]<<endl;
//   out.close();
  
  //! 3D output in txt, length renormed to meters
//   ofstream outSI;
//   outSI.open("H2-3dSI.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//       for(int k=0;k<grid;k++)
// 	outSI<<x[i]*2575000<<" "<<y[j]*2575000<<" "<<z[k]*2575000<<" "<<suv[i][j][k]<<endl;
//   outSI.close();
//   
//   
// //! 2d output
//   ofstream out_2d;
//   out_2d.open("H2-2dSI.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
// 	out_2d<<x[i]*2575000<<" "<<y[j]*2575000<<" "<<suv[i][j][grid/2]<<endl;
//   out_2d.close();
  
//!Binary output, for ptracer usage, The writing procedure must exactly match the one to be read in CHybrid_IonProfiles read_extern_IonProfile_pTracer()
//! BE CAREFUL CHANGING
//! Write header for 
  char Tracer_Info[42]="test"; //! this name doesn't really matter...
  short num_extern_Nodes_short[3]={grid,grid,grid};
  int TL=0;
  float extern_Origin_FILEREAL[3]={origin[0],origin[1],origin[2]};
  float extern_Length_FILEREAL[3]={boxdim,boxdim,boxdim};
  float Radius=1;
  float SI_Quantities[4]={1,1,1,1};//!
  
  short num_comps=1;
  int ftype=0;
  char Field_Name[50]="ion_prod_rate";
//   float field_array_box=

  
  ofstream out_bin;
//
//! These output file names can be changed
//! They'll be the name of the output file that is read into AIKEF
//
  out_bin.open("atm_H2",ofstream::binary);
  out_bin.write(reinterpret_cast<char*> (Tracer_Info),42*sizeof(char));
  out_bin.write(reinterpret_cast<char*> (num_extern_Nodes_short), 3*sizeof(short));
  out_bin.write(reinterpret_cast<char*> (&TL), sizeof(int));
  out_bin.write(reinterpret_cast<char*> (extern_Origin_FILEREAL), 3*sizeof(float));
  out_bin.write(reinterpret_cast<char*> (extern_Length_FILEREAL), 3*sizeof(float));
  out_bin.write(reinterpret_cast<char*> (&Radius), sizeof(float));
  out_bin.write(reinterpret_cast<char*> (SI_Quantities), 4*sizeof(float));
  
  out_bin.write( reinterpret_cast<char*> (&num_comps),sizeof(short));  //!Number of components of the field (1 or 3)
  out_bin.write( reinterpret_cast<char*> (&ftype),sizeof(int));  	//!field type 9,1,2,3 isnt even used after read?
  out_bin.write( reinterpret_cast<char*> (Field_Name), 50*sizeof(char));  //! also not used after read
  out_bin.write( reinterpret_cast<char*> (suv_array), num_comps*num_pts_mesh*sizeof(double));
  
  
  out_bin.close();

//! Output along single line (parallel to grid).
//! These are textfiles, can write a quick program to visualize 
//! ionization rates along a single line.
   ofstream lineoutx;
   lineoutx.open("O2_036be_quick-lineoutx.txt");
   for(int i=0;i<grid;i++)
     //if(sqrt(x[i]*x[i]+y[j]*y[j])>=1.9 && sqrt(x[i]*x[i]+y[j]*y[j])<=2.1)
     lineoutx<<x[i]<<" "<<suv[i][grid/2][grid/2]<<endl;
   lineoutx.close();
   
   ofstream lineouty;
   lineouty.open("O2_036be_quick-lineouty.txt");
   for(int i=0;i<grid;i++)
     //if(sqrt(x[i]*x[i]+y[j]*y[j])>=1.9 && sqrt(x[i]*x[i]+y[j]*y[j])<=2.1)
     lineouty<<y[i]<<" "<<suv[grid/2][i][grid/2]<<endl;
   lineouty.close();
  
//! Output cos_chi for checking
//   ofstream out_chi;
//   out_chi.open("cos_chi.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//     {
//       vec_x[0]=x[i];vec_x[1]=y[j];vec_x[2]=0.;
//       out_chi<<x[i]<<" "<<y[j]<<" "<<cos_chi(vec_x,sun)<<endl;
//     }
//   out_chi.close(); 
  
//   ofstream tau_out;
//   tau_out.open("tau.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//     {
//       vec_x[0]=x[i];vec_x[1]=y[j];vec_x[2]=0.;
//       tau_out<<x[i]<<" "<<y[j]<<" "<<tau[i][j][150]<<endl;
//     }
//   tau_out.close(); 
//   (x[i]*sun[1]+y[j]*sun[2]+a[3]*b[3])/sqrt((a[1]*a[1]+a[2]*a[2]+a[3]*a[3])*(b[1]*b[1]+b[2]*b[2]+b[3]*b[3]))
//   ofstream test;
//   test.open("dat4.txt");
// //       if(sqrt(x[i]*x[i]+y[j]*y[j])>=1.9 && sqrt(x[i]*x[i]+y[j]*y[j])<=2.1)
//       test<<x[40+16]<<" "<<z[40+16]<<" "<<suv2[40+16][0][40+16]<<endl;
//       test<<x[40+18]<<" "<<z[40+14]<<" "<<suv2[40+18][0][40+14]<<endl;
//   test.close();
  
  //cout<<"done"<<endl;   
  //finish = clock();
  //time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  //cout<< " Time: " << time/3600 << "hr." << endl << endl;   
  
  
  
}
