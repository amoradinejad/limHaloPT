

/** @file main.c Documented main module, including functions to initilize and cleanup the cosmology structure
 *  and examples of calls to functions in other modules to compute the line clustering and shot power spectrum. 
 *  
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 *  In order to call any function from the package, the function calls should be placed in the marked section of main() function.
 */


#include "header.h"
struct globals gb;


int main(int argc,char *argv[])
{	
 	
 	/**
	 * -----------------------------
	 * If you want to change the values of cosmological parmaeters compared to the 
	 * default values of the code, you should do so inside the initialize() function.
	 * For instance change the value of gb.h to set the value of hubble parmaeter in your cosmology.
	 * -----------------------------
	 */

 	/**
 	 * Initialize various ingredients for the rest of the code. More specifically, set the relavent directory paths, values of cosmological parameters and in general any parameter
 	 * that has to be passed to CLASS Boltzmann code. Furthermore initialize some of the interpolators needed
 	 * for the package.
 	 */
 	initialize();



 	/**
	 * -----------------------------
	 * Change if you want to assume a different mass function or minimum halo mass
	 * -----------------------------
	 */
 	
 	/**
 	 * Set the parameters for computing quantities related to spectral lines. 
 	 * More specifically, the mimum mass for halos that can host line-emitting galaxies M_min
 	 * and the the theoretical mass function to be used, ST: Sheth-Tormen, PSC: Press-Schecter, TK:: Tinker
 	 */
	double M_min   = 1.e9;
	long mode_mf   = ST;



 	/**
	 * -----------------------------
	 * Change if you want to limit the lines to include in your computation.
	 * -----------------------------
	 */
	/**
	 * Define which lines to compute, and how many points to have for z-interpolation
	 * 
	 * Note: The main line properties used in the calculation of the line power spectrum, which includes
	 * first and second moments of line luminosity (weighted by halo mass function) and 
	 * luminosity-weighted linear and quadratic line biases are computed in function Line_alloc_init().
	 * In that function, these quantities are calculate over an array of redshift values and then an interpolator 
	 * for these quantities is initialized. The parameter "ninterp" below, sets the number of interpolation points  
	 * for this interpolator. If you set it too low, the redshift-dependance of line properties would not 
	 * be accurate. In the tests done during code developement, ninterp = 150 give results at the precision 
	 * needed. However large values of ninterp results in longer run time.
	 * 
	 * For testing installation and running of the code, you can set lower value of ninterp, to make the code 
	 * run faster. But remember in your actuall runs, set ninterp to larger numbers. 
	 * 
	 * If you do not need to compute the line properties for all the 7 lines below, 
	 * you can set nlines to the number of lines you are interested, change the arrays of line[nlines] and JJ[nlines]
	 * to correspond to the lines you need to calculate. 
	 */
	size_t ninterp = 150;
	int nlines     = 7;
      int lines[7]   = {CO10, CO21, CO32, CO43, CO54, CO65, CII};
      int JJ[7]      = {1,2,3,4,5,6,0};




	/**
	 * -----------------------------
	 * In general, you should not need to change this portion of the code. But in order to call most functions of 
	 * limHaloPT, you furst need to initilize the Cosmology structure as done below. 
	 * -----------------------------
	 */
      
      /** 
       * Set values of cosmological parameter and finger-of-god paramater (peculiar velocity of line emitters) 
       * and initialize cosmo_pars array in the cosmology structure Cx_ref
       */
	double gb_pars[] = {gb.logAs,gb.ns,gb.h,gb.Omega_b,gb.Omega_cdm,gb.sigFOG0};

	struct Cosmology Cx_ref;
	for(int i=0;i<gb.Npars;i++){
		Cx_ref.cosmo_pars[i] = gb_pars[i];
	}

	/**
	 * Initialize the cosmology structure, which includes CLASS cosmology and Line structures
	 */
	clock_t tic_r = clock();
	Cosmology_init(&Cx_ref, gb.PS_kmax, gb.PS_zmax, nlines, lines, ninterp, M_min, mode_mf);
	printf("Reference Cosmology initialized\n");
	clock_t toc_r = clock();
	printf("Elapsed: %f seconds for ref cosmology\n", (double)(toc_r - tic_r)/CLOCKS_PER_SEC);




	/**
	 * -----------------------------
	 * Depending on what quantities ou want to compute, you need to modify this part of the main() function.
	 * -----------------------------
	 */

	/** 
	 * Set the k and z arrays
	 * 
	 * For testing the code, you can reduce the range of k, and z and also the number of points for each, nk, and nz.
	 */
 	int nk     = 200;
 	int nz     = 50;
	double *k  = loginit_1Darray(nk, 1.e-3, 2.);
 	double *z  = init_1Darray(nz,0.0,11.);
 	
 	double ps_clust_hm = 0.;
 	double ps_shot_hm  = 0.;
	/** 
	 * Compute the line clustering signal using halo model 
	 */
 	printf("Calculating halo-model line power spectrum\n");
 	int line_id = 0;
 	for(int i=0;i<nlines;i++){
 		line_id = i;
 		for(int j=0;j<nz;j++){
			for(int l=0;l<nk;l++){
				ps_clust_hm = PS_line_HM(&Cx_ref, k[l], z[j], M_min, mode_mf, lines[i], line_id);
				printf("%d %12.6e %12.6e %12.6e \n", i, z[j], k[l], ps_clust_hm);
			}
		}	
	}
	
	/**
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

 	int nd = 2;  //This is the number of digits to show in th evalue of z for sprintf()


  	printf("Calculating line mean bright temprature and halo-model shot nosie\n");
 	for(int i=0;i<nlines;i++){
 		for(int j=0;j<nz;j++){
 			//Set the path for the output of Tbar of the lines. By default the files are saved
 			//in "Output" directory. You can add subdirectories to Output directory by changing the path below if needed. 
 			sprintf(filename1,"%s/Tbar_J%d_z%.*f.txt", gb.output_dir, JJ[i], nd, z[j]);
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
			for(int l=0;l<nk;l++){
				ps_shot_hm = PS_shot_HM(&Cx_ref, k[l], z[j], M_min, input, mode_mf, lines[i]);
				printf("%d %12.6e %12.6e %12.6e \n", i, z[j], k[l], ps_shot_hm);

			}

			// When neglecting the nfw profile for comparison with sims, the shot noise is scale-independent. 
			// So uncomment the following line and comment out the above for loop to just compute it for a single k
			// PS_shot_HM(&Cx_ref, k[0], z[j], M_min, input, mode_mf, lines[i]);
			
		}	
	}	
	free(k);
 	free(z);

	/** 
	 * -----------------------------
	 * Modify up to here 
	 * -----------------------------
	 */

	cleanup(&Cx_ref);
	
	return 0;
}



	
 	
	
