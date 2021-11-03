#include "header.h"

struct globals gb;


int main(int argc,char *argv[])
{	
	extern struct globals gb;
 	
 	initialize();

	long fish_dim;
	long int i,j,l, q, s,t;

	fish_dim = gb.Npars;	

	double pk_zmax = 14.;
	double pk_kmax = 200.;
	double inc     = 0.01;
	double M_min   = 1.e9;
	long mode_mf   = ST;

	/*
	* Define which lines to compute, and how many points to have for z-interpolation
	*/
	int nlines     = 7;
        int lines[7]   = {CO10, CO21, CO32, CO43, CO54, CO65, CII};
        int JJ[7]      = {1,2,3,4,5,6,0};
        size_t ninterp = 150;

	double gb_pars[] = {gb.logAs,gb.ns,gb.h,gb.Omega_b,gb.Omega_cdm,gb.sigFOG0};

	struct Cosmology Cx_ref;
	for(i=0;i<gb.Npars;i++){
		Cx_ref.cosmo_pars[i] = gb_pars[i];
	}

	/*
	* Initialize the cosmology structure, which includes CLASS cosmology and Line structures
	*/
	clock_t tic_r = clock();
	Cosmology_init(&Cx_ref, pk_kmax, pk_zmax, nlines, lines, ninterp, M_min, mode_mf);
	printf("Reference Cosmology initialized\n");
	clock_t toc_r = clock();
	printf("Elapsed: %f seconds for ref cosmology\n", (double)(toc_r - tic_r)/CLOCKS_PER_SEC);

	
	/*
	* -----------------------------
	* Depending on what quantities ou want to compute, you need to modify this part of the main() function.
	* -----------------------------
	*/

	/* 
	 * Set the k and z arrays
	 */
      	int nk = 200;
  	int nz = 50;
	double *k  = loginit_1Darray(nk, 1.e-3, 8.);
 	double *z = init_1Darray(nz,0.0,11.);
	
	/* 
	* Compute the line clustering signal using halo model 
	*/
 	int line_id = 0;
 	for(i=0;i<nlines;i++){
 		line_id = i;
 		for(j=0;j<nz;j++){
			for(l=0;l<nk;l++){
				PS_line_HM(&Cx_ref, k[l], z[j], M_min, mode_mf, lines[i], line_id);
			}
		}	
	}
	
	/* 
	* Compute the shot noise (beyond the Poisson limit), using halo model
	*/
	double Omegam   = (Cx_ref.cosmo_pars[3] + Cx_ref.cosmo_pars[4]);
	double rhom_bar = Omegam * rhoc(&Cx_ref,0.); 
	double lbias_arr[2];
	double input[5];
	double b1   = 0.;
	double b2   = 0.;
	double Tbar = 0.;

	FILE *fp1;
  	char filename1[FILENAME_MAX];
 	for(i=0;i<nlines;i++){
 		for(j=0;j<nz;j++){
 			sprintf(filename1,"%s/line/Tbar_J%d_z%d.txt", gb.output_dir, JJ[i],(int)z[j]);
	 		fp1    = fopen(filename1, "ab");

	 		Tbar = Tbar_line(&Cx_ref, i, z[j]);
			line_bias(Cx_ref.Lines[i],z[j],lbias_arr);
			b1 = lbias_arr[0];
			b2 = lbias_arr[1];
			fprintf(fp1,"%d %12.6e %12.6e %12.6e %12.6e \n",i, z[j], b1, b2, Tbar); 

			input[0] = Tbar_line(&Cx_ref, i,z[j]); 
			input[1] = b1;
			input[2] = 0.25 * pow(b2,2.) * pow(input[0],2.) * b22_ls(&Cx_ref,z[j]);
			input[3] = PS_shot(&Cx_ref, z[j], i);
			input[4] = rhom_bar; 

			//If accounting for the nfw profile, the shot noise will be scale-dependent. So loop over k-values
			// for(l=0;l<nk;l++){
			// 	PS_shot_HM(&Cx_ref, k[l], z[j], M_min, input, mode_mf, lines[i]);
			// }

			// When neglecting the nfw profile for comparison with sims, the shot noise is scale-independent. 
			// So uncomment the following line and comment out the above for loop to just compute it for a single k
			PS_shot_HM(&Cx_ref, k[0], z[j], M_min, input, mode_mf, lines[i]);
			
		}	
	}	
	free(k);
 	free(z);

	/* 
	* -----------------------------
	* Modify up to here 
	* -----------------------------
	*/

	cleanup(&Cx_ref);
	
	return 0;
}


void initialize()
{
	extern struct globals gb;

	getcwd(gb.project_home , sizeof(gb.project_home));  //This gives the path to the source directory
	chdir("..");  ///Change the path to the parent directory
	getcwd(gb.project_home , sizeof(gb.project_home));

	sprintf(gb.output_dir,"%s/Output", gb.project_home);
	sprintf(gb.data_dir,"%s/Input", gb.project_home);
	sprintf(gb.data_priors,"%s/Priors", gb.data_dir);

	sprintf(gb.SFR_filename,"%s/release-sfh_z0_z8_052913/sfr/sfr_release.dat",gb.data_dir);  

	logSFR_alloc_init();

	///// cosmological parameters corresponding to initial conditions of HV simulations
	gb.c 		 = 2.99792458e5;  /// In units of km/s
	gb.h 		 = 0.677;
	gb.Omega_cdm     = 0.11923/pow(gb.h,2.);  ////omega_cdm = Omega_cdm h^2 ;
	gb.Omega_b 	 = 0.02247/pow(gb.h,2.);   ///omega_b = Omega_b h^2;
	gb.Omega_r 	 = 0.0000910958;  ////radiation = photons + neutrinos
	gb.Omega_g 	 = 5.3956421715871286e-05;   ////photons, input for Class
	gb.Omega_nu  = 0.00;    ////neutrinos
	gb.ns        = 0.96824;
	gb.As 	     = 2.1085e-9;
	gb.logAs     = log(gb.As*pow(10.,10.));  ///3.0665
	gb.kp 	     = 0.05;  //in unit of Mpc^-1

	gb.sigFOG0 = 250.; 	

	gb.PS_kmin = 1.e-4;
	gb.PS_kmax = 200.;

	gb.Npars   = 6;
	
	return;
}


void cleanup(struct Cosmology *Cx)
{
	extern struct globals gb;	

	Cosmology_free(Cx);
	SFR_Behroozi_free();

	return;
}
	
 	
	
