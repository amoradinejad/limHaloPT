

/** @file main.c Documented main module, including functions to initilize and cleanup the cosmology structure
 *  and calls to fthree main functions to compute the line mean brightness temprature, linear and quadratic biases 
 * and clustering and shot contributions to line power spectrum. 
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
 	 * Initialize various ingredients for the rest of the code. 
 	 * 
 	 * More specifically, initilize the filename and interpolator of SFR.
 	 * Read in the values of the parameters from .ini file and set the values of the elements of the glocal structure
 	 * Furthermore, the initialize() function returns a structure containing parameters such as as mimumum halo mass that are read from the .ini file.
 	 */

	struct initialize_struct *init_struct  = (struct initialize_struct *)malloc(sizeof(struct initialize_struct));

	initialize(argv,init_struct);

 	/**
 	 * Extract from the return struct of initialize() function, values of the parameters for computing quantities related to spectral lines. 
 	 * Note that these values were read from the .ini file within the initialize() function.
 	 * 
 	 * M_min : mimum mass for halos that can host line-emitting galaxies in unit of M_sun
 	 * mode_mf: switch to set the theoretical mass function to be used, ST: Sheth-Tormen, PSC: Press-Schecter, TK:: Tinker
 	 * ninterp: how many points to have for z-interpolation of line bias and first and second luminosity moments
 	 * nlines: the number of emission lines you want to include in the analysis
 	 * 
 	 * Note about ninterp: The main line properties used in the calculation of the line power spectrum, which includes
	 * first and second moments of line luminosity (weighted by halo mass function) and 
	 * luminosity-weighted linear and quadratic line biases are computed in function Line_alloc_init().
	 * In that function, these quantities are calculate over an array of redshift values ranging from 0 to line_zmax (a paramter of .ini file). 
	 * Subsequently, an interpolator for these quantities is initialized. The parameter "ninterp" below, sets the number of interpolation points  
	 * for this interpolator. If you set ninterp too low and the redshift range is wide, the redshift-dependance of line properties would not 
	 * be accurate. In the tests done during code developement, ninterp = 150 give results at the precision 
	 * needed when comparing with MithraLIMSims. However large values of ninterp results in longer run time. For testing installation and running the code using test.c module, 
	 * lower value of ninterp is used, to have the code run fast. But remember in your actuall runs, set ninterp to larger numbers. 
 	 *
 	 * Note about nlines: If you do not need to compute the line properties for all the 7 lines implemented in the code, 
	 * you can set nlines in the .ini file to the number of lines you are interested. Note that at the moment 
	 * you can only call the lines in the order below. 
	 *
 	 */

	double M_min   = init_struct->Mh_min;
	long mode_mf   = init_struct->mode_mf;
 	size_t ninterp = init_struct->ninterp;
 	int nlines     = init_struct->nlines;

	/**
	 * Given the number of lines you chose with nline, the line_type[] array and JJ[] arrays will be initialized here. 
	 */

      int J_val;
      int JJ[nlines];
      int line_type[nlines];

 	for(int i=0;i<nlines;i++){
 		line_type[i] = init_struct->lines[i];
 		printf("line type %d %d \n", init_struct->lines[i], line_type[i]);
 	}


      /**
       * set the value of J for the lines. Note that value of J is onlyy important for determining the frequency of CO rotational lines
       * and for CII, J=0 has no physical implication. Here we set the values of J just for labeling of the output files
       * */
      for(int i=0;i<nlines;i++){
      	if(line_type[i] == CO10) 
	      	J_val = 1;
	    	else if(line_type[i] == CO21) 
	      	J_val = 2;
	    	else if(line_type[i] == CO32) 
	      	J_val = 3;
	    	else if(line_type[i] == CO43)
	      	J_val = 4; 
	    	else if(line_type[i] == CO54)
	      	J_val = 5; 
	    	else if(line_type[i] == CO65)
	      	J_val = 6; 
	      else if(line_type[i] == CII)
	      	J_val = 0;
      	JJ[i]    = J_val;
      	printf("check that lines are read and assigned correctly %d %d \n", line_type[i], JJ[i]);;
      }

      
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
	Cosmology_init(&Cx_ref, gb.PS_kmax, gb.PS_zmax, nlines, line_type, ninterp, M_min, mode_mf);
	printf("Reference Cosmology initialized\n");
	clock_t toc_r = clock();
	printf("Elapsed: %f seconds for ref cosmology\n", (double)(toc_r - tic_r)/CLOCKS_PER_SEC);


	/**
	 * -----------------------------
	 * At the moment, limHaloPT has three possible outputs; the mean brightness temprature, the clustering power or the shot power.
	 * The latter two are computed within halo model. The clustering power is real-space (RSD not included).
	 * 
	 * If you want to call any other function of limHaloPT, you should add them to to this section of the code 
	 * -----------------------------
	 */

 	if(gb.switch_Tbar == MEAN){
 		printf("You asked to compute the line mean bright temprature \n");

 		int nz_mean     = gb.mean_nz;
 		double *z_mean  = init_1Darray(nz_mean,0.0,gb.mean_zmax);

		double lbias_arr[2];
		double b1   = 0.;
		double b2   = 0.;
		double Tbar = 0.;

 		FILE *fp1;
  		char filename1[FILENAME_MAX];

	 	for(int i=0;i<nlines;i++){
	 		sprintf(filename1,"%s/Tbar_J%d.txt", gb.output_dir, JJ[i]);
			fp1    = fopen(filename1, "ab");
			fprintf(fp1,"%s %s %s %s %s \n","#", "z", "b1", "b2", "Tbar");
	 		for(int j=0;j<nz_mean;j++){
		 		Tbar = Tbar_line(&Cx_ref, i, z_mean[j]);
				line_bias(Cx_ref.Lines[i],z_mean[j],lbias_arr);
				b1   = lbias_arr[0];
				b2   = lbias_arr[1];
				fprintf(fp1,"%12.6e %12.6e %12.6e %12.6e \n",z_mean[j], b1, b2, Tbar); 
			}
			fclose(fp1);
		}
		free(z_mean);
	}

	if(gb.switch_HMshot == HMclust || gb.switch_HMshot == HMshot){

	 	int nz_power     = gb.power_nz;
 		double *z_power  = init_1Darray(nz_power,0.0,gb.power_zmax);

		if(gb.switch_HMclust == HMclust){
 			printf("You asked to compute the line clustering power within halo model \n");
	 	
		 	FILE *fp2, *fp3, *fp4;
	  		char filename2[FILENAME_MAX];
	 		int nk_power     = gb.power_nk;
			double *k_power  = loginit_1Darray(nk_power, gb.power_kmin, gb.power_kmax);

			double ps_clust_hm = 0.;
	 		int nd = 2;  //This is the number of digits to show in th evalue of z for sprintf()

		 	for(int i=0;i<nlines;i++){
		 		for(int j=0;j<nz_power;j++){
		 			sprintf(filename2,"%s/clust_J%d_z%.*f.txt", gb.output_dir, JJ[i], nd, z_power[j]);
		 			fp2    = fopen(filename2, "ab");
					fprintf(fp2,"%s %s %s %s %s \n","#", "line_id", "z", "k", "ps_clust-[Mpc^3]");
					
					for(int l=0;l<nk_power;l++){
						ps_clust_hm = PS_line_HM(&Cx_ref, k_power[l], z_power[j], M_min, mode_mf, line_type[i], i);
						fprintf(fp2,"%d %12.6e %12.6e %12.6e \n", i, z_power[j], k_power[l], ps_clust_hm);
					}
					fclose(fp2);
				}	
			}
			free(k_power);
		}
		else if (gb.switch_HMshot == HMshot){
			printf("You asked to compute the non-Poisson shot noise of the line \n");

			FILE *fp5, *fp6;
	  		char filename5[FILENAME_MAX];

			double Omegam   	= (Cx_ref.cosmo_pars[3] + Cx_ref.cosmo_pars[4]);
			double rhom_bar 	= Omegam * rhoc(&Cx_ref,0.); 
			double input[3];
		 	double ps_shot_hm = 0.;

		 	double k_shot = 0.01;
	 		int nd = 2;  //This is the number of digits to show in th evalue of z for sprintf()

		 	for(int i=0;i<nlines;i++){
		 		sprintf(filename5,"%s/shot_J%d.txt", gb.output_dir, JJ[i]);
				fp5    = fopen(filename5, "ab");
				fprintf(fp5,"%s %s %s %s %s %s \n","#", "line_id", "z", "ps_shot [Mpc^3]");

		 		for(int j=0;j<nz_power;j++){
					input[0]   = Tbar_line(&Cx_ref, i, z_power[j]); 
					input[1]   = rhom_bar; 
					input[2]   = PS_shot(&Cx_ref, z_power[j], i);
					ps_shot_hm = PS_shot_HM(&Cx_ref, k_shot, z_power[j], M_min, input, mode_mf, line_type[i], i);
					fprintf(fp5,"%d %12.6e %12.6e \n", i, z_power[j], ps_shot_hm);
				}
				fclose(fp5);
				fclose(fp6);

			}	
		}
		free(z_power);
	}

	cleanup(&Cx_ref);
	
	return 0;
}



	
 	
	
