
/** @file test.c This module includes two test functions that check the output of teh clustering and shot contributions
 * 
 * Azadeh Moradinezhad Dizgah, April 3rd 2022
 *
 */


#include "test_header.h"
struct globals gb;


int main()
{

      initialize();

      gb.line_zmax = 1.;  //To carryy on the test faster, re-write the value of maximum z for line calculatoin 
                          // so that you can run the code for less number of interpolation points ninterp.     
      
      double M_min   = 1.e9;
      long mode_mf   = ST;

      size_t ninterp = 10;
      int nlines     = 1;
      int lines[1]   = {CO10};

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


      double k    = 0.01;
      double z    = 0.2;
      int line_id = 0;

      double ps_clust_hm = PS_line_HM(&Cx_ref, k, z, M_min, mode_mf, lines[0], line_id);
      printf("ps_clust_hm = %12.6e \n", ps_clust_hm);


      double lbias_arr[2];
      line_bias(Cx_ref.Lines[0],z,lbias_arr);
      double b1   = lbias_arr[0];
      double b2   = lbias_arr[1];
      double Tbar = Tbar_line(&Cx_ref, line_id, z);
      printf("b1 = %12.6e b2 = %12.6e Tbar = %12.6e \n", b1, b2, Tbar); 


      double Omegam   = (Cx_ref.cosmo_pars[3] + Cx_ref.cosmo_pars[4]);
      double rhom_bar = Omegam * rhoc(&Cx_ref,0.); 
      
      double input[5];
      input[0] = Tbar;
      input[1] = b1;
      input[2] = 0.25 * pow(b2,2.) * pow(input[0],2.) * b22_ls(&Cx_ref,z);
      input[3] = PS_shot(&Cx_ref, z, line_id);
      input[4] = rhom_bar; 

      double ps_shot_hm = PS_shot_HM(&Cx_ref, k, z, M_min, input, mode_mf, lines[0]);
      printf("ps_shot_hm = %12.6e \n", ps_shot_hm);


      cleanup(&Cx_ref);
      

      return 0;

}

