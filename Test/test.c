
/** @file test.c This module includes two test functions that check the output of teh clustering and shot contributions
 * 
 * Azadeh Moradinezhad Dizgah, April 3rd 2022
 *
 */


#include "test_header.h"
struct globals gb;


int main(int argc,char *argv[] )
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

      
      double M_min   = init_struct->Mh_min;
      long mode_mf   = init_struct->mode_mf;
      size_t ninterp = init_struct->ninterp;
      int nlines     = init_struct->nlines;

      /**
       * Given the number of lines you chose with nline, set line_type[]
       */

      int line_type[nlines];
      for(int i=0;i<nlines;i++){
            line_type[i] = init_struct->lines[i];
            printf("line type %d %d \n", init_struct->lines[i], line_type[i]);
      }


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

      double k    = 0.01;
      double z    = 0.1;
      int line_id = 0;

      /**
       * Compute the line biases and mean brightness temprature for a fixed value of redshift and only for CO10  
       */
      double lbias_arr[2];
      line_bias(Cx_ref.Lines[line_id],z,lbias_arr);
      double b1   = lbias_arr[0];
      double b2   = lbias_arr[1];
      double Tbar = Tbar_line(&Cx_ref, line_id, z);
      printf("SUCESS Tbar and biases: b1 = %12.6e b2 = %12.6e Tbar = %12.6e \n", b1, b2, Tbar); 

      /**
       * Compute the clustering contribution of line power spectrum within halo model
       */
      double ps_clust_hm = PS_line_HM(&Cx_ref, k, z, M_min, mode_mf, line_type[0], 0);
      printf("SUCESS Power_clust:  ps_clust_hm = %12.6e \n", ps_clust_hm);


      /**
       * Compute the shot contribution of line power spectrum within halo model
       */
      double Omegam   = (Cx_ref.cosmo_pars[3] + Cx_ref.cosmo_pars[4]);
      double rhom_bar = Omegam * rhoc(&Cx_ref,0.); 
      
      double input[3];
      input[0] = Tbar;
      input[1] = rhom_bar; 
      input[2] = PS_shot(&Cx_ref, z, line_id);
     
      double ps_shot_hm = PS_shot_HM(&Cx_ref, k, z, M_min, input, mode_mf, line_type[0], 0);
      printf("SUCESS POwer_shot: ps_shot_Poisson= %12.6e ps_shot_hm = %12.6e \n", input[2], ps_shot_hm);

      /**
       *  Free the cosmology structure
       */
      cleanup(&Cx_ref);
      

      return 0;

}

