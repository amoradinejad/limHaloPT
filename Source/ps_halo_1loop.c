

/** @file ps_halo_1loop.c Documented real-space, direct integration computation of 1loop contributions of the halo/galaxy power spectrum
 * See arXiv:2010.14523 for explicit expressions 
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 * This module computes the 1loop halo/galaxy power sprtcurm in real-space via direct numerical integration.
 * IR-resummation and EFT counter terms are included. In addition to loops due to gravitational loops, terms arising only in the presence of local PNG are 
 * also included. The explicit expressions of all the loops are given in 2010.14523. 
 *
 * In summary, the following functions can be called from other modules:
 * -# PS_hh_G()
 * -# PS_hh_PNG()
 * -# Compute_Gloops()
 * -# Compute_PNGloops()
 * -# F2_s()
 * -# F3_s()
 * -# S2_s()
 * -# F2()
 * -# S2()
 */

#include "header.h"
struct globals gb;

/**
 * Compute the contributions up to 1loop to halo power spectrum for Gaussian initial conditions
 *
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber 
 * @param z            Input: redshift of interest
 * @param M            Input: halo mass, used in computing the halo bias
 * @param mode_pt      Input: switch to decide whether to compute tree-level halo power spectrum or the 1loop 
 * @param IR_switch    Input: switch to decide whether to perform IR resummation or no
 * @param SPLIT        Input: switch to set the method to perform the wiggle-nowiggle split of matter power spectrum
 * @param mode_mf      Input: switch to set the theoretical model of the mass function used to compute the halo biases
 * @return G loop contributions of P_h
 */

double PS_hh_G(struct Cosmology *Cx, double k,  double z, double M, long mode_pt, long IR_switch, long SPLIT, long mode_mf)
{
      
      int cleanup = 0;

      double pm_lin = 0., pm_lin_IR = 0., pm_1loop_IR = 0., pm_22 = 0., pm_13 = 0., pm_1loop =0., pm_ct = 0., ph_tot = 0.;

      double *bias_arr = make_1Darray(4);
      halo_bias(Cx, M, z, mode_mf, bias_arr);
      double b1  = bias_arr[0];
      double b2  = bias_arr[1];
      double bG2 = bias_arr[2];
      double btd = bias_arr[3];

      printf("Bias Vals: %12.6e %12.6e %12.6e %12.6e %12.6e  %12.6e \n",k,M,b1,b2,bG2,btd);

      if(mode_pt == LOOP){
            double *ps_hloops, *ps_mloops;
            ps_hloops = make_1Darray(6);
            ps_mloops = make_1Darray(2);

            if(IR_switch == WIR){
                Compute_G_loops(Cx, k, z, WIR, HALO, SPLIT, ps_hloops);
                Compute_G_loops(Cx, k, z, WIR, MATTER, SPLIT,ps_mloops);
            }
            else if(IR_switch == NOIR){
                Compute_G_loops(Cx, k, z, NOIR, HALO, SPLIT, ps_hloops);
                Compute_G_loops(Cx, k, z, NOIR, MATTER, SPLIT,ps_mloops);
            }  

            double pb1b2    = b1 * b2  * ps_hloops[0];
            double pb1bg2   = 2. * b1 * bG2 * ps_hloops[1];
            double pb22     = 0.25 * pow(b2, 2.) * ps_hloops[2];
            double pbg22    = pow(bG2, 2.)  * ps_hloops[3];
            double pb2bg2   = b2 * bG2 * ps_hloops[4];
            double pb1b3nl  = 2. * b1 * (bG2 + 2./5. * btd) * ps_hloops[5];
            double ph_loops =  pb1b2 + pb1bg2 + pb22+ pbg22 + pb2bg2 + pb1b3nl;

            double cs2      = 0.2;
            double khat     = 1. * gb.h;
            
            if(IR_switch == NOIR){
                pm_lin   = Pk_dlnPk(Cx, k, z, LPOWER);
                pm_22    = ps_mloops[0];
                pm_13    = ps_mloops[1]; 
                pm_1loop = pm_lin + pm_22 + pm_13;
                pm_ct    = - 2. * cs2 * pow(k, 2.) * pm_lin;
                ph_tot   = (pow(b1, 2.) * (pm_1loop + pm_ct) + ph_loops); 
            }
            else if(IR_switch == WIR){
                pm_lin      = Pk_dlnPk(Cx,k, z, LPOWER);
                pm_22       = ps_mloops[0];
                pm_13       = ps_mloops[1]; 
                pm_lin_IR   = pm_IR_LO(Cx,k, z, SPLIT);
                pm_1loop_IR = pm_IR_NLO(Cx,k, z, SPLIT);
                pm_ct       = - 2. * cs2 * pow(k, 2.) * pow(k, 2.)/(1.+pow(k/khat,2.))* pm_lin_IR;
                ph_tot      = (pow(b1, 2.) * (pm_1loop_IR + pm_ct) + ph_loops);
            }  

            FILE *fp;
            char filename1[FILENAME_MAX];
            sprintf(filename1,"%s/ph_IR_z%d_comps.txt", gb.output_dir,(int)z);

            fp = fopen(filename1, "ab");
            fprintf(fp, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n"\
                    ,k, pow(b1, 2.) * pm_lin, pow(b1,2.)* pm_1loop_IR, pow(b1,2.) *pm_ct, pow(b1,2.) * pm_22, pow(b1,2.) *pm_13, pb1b2, pb1bg2, pb22, pbg22, pb2bg2, pb1b3nl,ph_tot);
            fclose(fp);
            
            free(ps_hloops);
            free(ps_mloops);         
      }
      else if(mode_pt == TREE){
            ph_tot  = pow(b1, 2.) * Pk_dlnPk(Cx, k, z, LPOWER);
      }

      free(bias_arr);

      return ph_tot; 
}


/**
 * Compute contributions up to 1loop to halo power spectrum arising from non-Gaussian initial conditions of local shape
 *
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber 
 * @param z            Input: redshift of interest
 * @param M            Input: halo mass, used in computing the halo bias
 * @param mode_pt      Input: switch to decide whether to compute tree-level halo power spectrum or the 1loop 
 * @param IR_switch    Input: switch to decide whether to perform IR resummation or no
 * @param SPLIT        Input: switch to set the method to perform the wiggle-nowiggle split of matter power spectrum
 * @param mode_mf      Input: switch to set the theoretical model of the mass function used to compute the halo biases
 * @return PNG loop contributions of P_h
 */

double PS_hh_PNG(struct Cosmology *Cx, double k, double z, double M, long mode_pt, long IR_switch, long SPLIT, long mode_mf)
{

      double f          = 0.;
      double deltac     = 1.686;
      double pm_lin = 0., pm = 0., ph_loops = 0., ph_tot = 0.;

      static int cleanup = 0;
  
      double *bias_arr = make_1Darray(4);
      halo_bias(Cx, M, z, mode_mf, bias_arr);
      double b1  = bias_arr[0];
      double b2  = bias_arr[1];
      double bG2 = bias_arr[2];
      double btd = bias_arr[3];

      double bz  = 2.* deltac * (b1-1.);
      double bzd = 2.* ((1.-b1) + deltac * (b2 - 8./21.*(b1-1.))) + bz;
      double fnl = Cx->cosmo_pars[5L];

      double tk_inv  = 1./Mk_dlnMk(Cx, k, z, TRANS);

      if(mode_pt == LOOP){
            double *ps_hloops, *ps_mloops;
            ps_hloops  = make_1Darray(8);
            ps_mloops  = make_1Darray(2);

            if(IR_switch == WIR){
                Compute_PNG_loops(Cx, k, z, WIR, SPLIT, ps_hloops);
                Compute_G_loops(Cx, k, z, WIR, MATTER, SPLIT,ps_mloops);
            }
            else if(IR_switch == NOIR){
                Compute_PNG_loops(Cx, k, z, NOIR, SPLIT, ps_hloops);
                Compute_G_loops(Cx, k, z, NOIR, MATTER, SPLIT,ps_mloops);
            } 

            double pm_22   = ps_mloops[0];   
            double pm_13   = ps_mloops[1];   

            //////Comment out if computing the averaged model from the chains
            double pb1bz   = 2. * b1  * bz  * ps_hloops[0];
            double pb1bzd  = 2. * b1  * bzd * ps_hloops[1];
            double pb2bz   = b2 * bz  * ps_hloops[2];
            double pbg2bz  = 2. * bG2 * bz  * ps_hloops[3];
            double pb2bzd  = b2 * bzd * ps_hloops[4];
            double pbg2bzd = 2. * bG2 * bzd * ps_hloops[5];
            double pb3nlbz = 2. * bz  * (bG2 + 2./5. * btd) * pm_IR_LO(Cx, k, z, SPLIT) * tk_inv * ps_hloops[6]; 
	          double pb1bz_fnl2 = 2. * b1 * bz * pm_IR_LO(Cx, k, z, SPLIT) * pow(tk_inv,2.) * ps_hloops[7];

            if(IR_switch == NOIR)
                pm = Pk_dlnPk(Cx, k, z, LPOWER) + pm_13;
            else if(IR_switch == WIR)
                pm = pm_IR_NLO(Cx, k, z, SPLIT) - pm_22;  
                
            ph_tot = 2. * fnl * b1 * bz * tk_inv * pm + pow(fnl * bz * tk_inv, 2.) * pm_IR_LO(Cx, k, z, SPLIT) + ph_loops;            

            // FILE *fp;
            // // fp = fopen("results/ps_bis_th/z1/NGps_M1_p4_z11_comps.txt", "a");
            // fprintf(fp, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n"\
            //         ,k, ph_tot, 2. * b1bz * tk_inv * pm, bz2*pow(tk_inv, 2.)*pm_IR_LO(k, z, Cx, SPLIT)\
            //         ,pb1bz, pb1bzd, pb2bz, pbg2bz , pb2bzd, pbg2bzd, pb3nlbz, pb1bz_fnl2);
            // fclose(fp);

            free(ps_hloops);
            free(ps_mloops);
           
      }
      else if(mode_pt == TREE){
            pm_lin = Pk_dlnPk(Cx, k, z, LPOWER);
            ph_tot = 2. * fnl * b1 * bz * tk_inv * pm_lin + pow(fnl * bz * tk_inv,2.) * pm_lin;
      }
      
      free(bias_arr);

      // printf("ok, non-Gaussian piece computed\n");

      return ph_tot; 
}




/**
 * Compute the loop contributions dure to nonlinear evolution of matter fluctuations and nonlinear halo bias, present for Gaussian initial conditions
 * The function G_loop_integrands() defines the integrand and Compute_G_loops() computes the integrals 
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber 
 * @param z            Input: redshift of interest
 * @param M            Input: halo mass, used in computing the halo bias
 * @param IR_switch    Input: switch to decide whether to perform IR resummation or no
 * @param hm_switch    Input: switch to decide whether to compute the 1loop terms due to matter or bias
 * @param SPLIT        Input: switch to set the method to perform the wiggle-nowiggle split of matter power spectrum
 * @param result       Output: an output array containing the results of the 1loop terms, 
 *                             has 2 elements for hm_switch=MATTER, and 6 elements for hm_switch=HALO
 * @return void
 */

void Compute_G_loops(struct Cosmology *Cx, double k, double z, long IR_switch, long hm_switch,long SPLIT, double *result)
{
      struct integrand_parameters2 par;

      int ncomp = 0;

      if(hm_switch == HALO)
        ncomp = 6;
      else if(hm_switch == MATTER)
        ncomp = 2;

      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];

      double plin_IR_k   = pm_IR_LO(Cx, k, z, SPLIT);

      par.p1  = Cx;
      par.p4  = k;
      par.p5  = z;
      par.p6  = log(1.e-4);
      par.p7  = log(100.);
      par.p13 = IR_switch;
      par.p14 = hm_switch;
      par.p15 = SPLIT;
      par.p8  = plin_IR_k;  /* Note: since cuhre integrator is parallelized, if evaluating all the 
                             * IR resummed power spectra in the integrand, it will build the interpolators
                             * For each thread and creats a mess. Here we evaluate the p(k) (which builds the interpolator)
                             * The first time it is called, and then call the p(q) and p(kmq) and p(kpq) inside the integrand
                             */

      // int ndim = 2, nvec = 1, verbose = 0, last = 4, seed = 0, 
      //     mineval = 0, maxeval = 2e6, key1 =50, key2 = 50, key3 = 1, maxpass = 200, 
      //       border =0, maxchisq = 10, mindeviation=0.25,
      //       ngiven = 0, ldxgiven = ndim, nextra = 0; 
      // int nregions, neval;   

      // Divonne(ndim, ncomp, G_loop_integrands, &par, nvec,
      //           RelErr, AbsErr, verbose, seed,
      //          mineval, maxeval,  key1, key2, key3, maxpass,
      //         border, maxchisq, mindeviation,
      //         ngiven, ldxgiven, NULL, nextra, NULL,
      //        NULL, NULL, &nregions, &neval, fail, result, error, prob);

      int key;
      if(hm_switch == MATTER)
          key = 7;
      else if(hm_switch == HALO)
          key = 13;

      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 1e8;
      int nregions, neval;
      
      Cuhre(ndim,ncomp, G_loop_integrands, &par, nvec,
              RelErr, AbsErr, verbose | last,
               mineval, maxeval, key,
               NULL, NULL, &nregions, &neval, fail, result, error, prob);
      
      // for(int i =0; i<ncomp; i++)
      //     printf("Gloops integral : %d %12.6e %12.6e %12.6e %12.6e %d \n", i, k, result[i], error[i], prob[i], fail[i]);  


      return;
}

static int G_loop_integrands(const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)    
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);      

      struct Cosmology *Cx = pij.p1;
      double k            = pij.p4;
      double z            = pij.p5;
      double logqmin      = pij.p6;
      double logqmax      = pij.p7;  
      long IR_switch      = pij.p13;
      long hm_switch      = pij.p14;
      long SPLIT          = pij.p15;
      double plin_IR_k   = pij.p8;


      double logq  = (x[0] * (logqmax - logqmin) + logqmin); 
      double cos   = (2. * x[1] - 1.); 
      double q     = exp(logq);         
      double kmq   = pow(fabs(pow(q, 2.) + pow(k, 2.) - 2. * q * k * cos), 1./2.);
      double kpq   = pow(fabs(pow(q, 2.) + pow(k, 2.) + 2. * q * k * cos), 1./2.);

      static int cleanup = 0;

      /// Model used in 1907.06666, the integrals are given in the appendix, Eq. A1, note that my S2_s = sigma^2(q,k-1) and F2_s = F2(q,k-q) in their notation.
      /// Factor of 2. * (logqmax - logqmin) is due to change of variable from 0 to logarithmic k, and a factor of 2*PI is due to integration over azimuthal angle. 
      /// Note that to compare the theoretical predictions against Emiliano's measurement, since he is using a different notation for Fourier transform, I need to 
      /// devide each 0 power spectrum by a factor of 1/pow(2.*M_PI,3.), which I do in my pk_lin() function. If using another notation for Fourier transform 
      /// (the one that I usually use, which has a factor of 1/pow(2*M_PI,3) in the definition), you need to multiply these integrands by a factor of 1/pow(2*M_PI,3).
      
      /// The integrands below correspond to the follwing bias combinaions: 
      /*  ff[0] : b1b2
          ff[1] : b1bG2
          ff[2] : b2b2
          ff[3] : bG2bG2
          ff[4] : b2bG2
          ff[5] : b1b3nl, b3nl = 2*(bG2+2/5*btd)

          ff[0] : P22_m
          ff[1] : P13_m
      */

      ////Defining regularized version of p22 integrand as in Senatore et al. We need to first define heavisde step function

      double theta_kmq = 0, theta_kpq = 0.;
      if(kmq>=q)
        theta_kmq = 1.;
      else
        theta_kmq = 0.;

      if(kpq>=q)
        theta_kpq = 1.;
      else
        theta_kpq = 0.;

      int nn =0;
      if(hm_switch == HALO)
         nn = 6;
      else if(hm_switch == MATTER)
         nn = 2; 
       
      if(kmq <= exp(logqmax) && kmq>=exp(logqmin) && kpq <= exp(logqmax) && kpq>=exp(logqmin)){
        if(IR_switch == NOIR){
            double plin_k   =  Pk_dlnPk(Cx, k, z, LPOWER);
            double plin_q   =  Pk_dlnPk(Cx, q, z, LPOWER);
            double plin_kmq =  Pk_dlnPk(Cx, kmq, z, LPOWER);
            double plin_kpq =  Pk_dlnPk(Cx, kpq, z, LPOWER);

            if(hm_switch == HALO)
            {
                ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * F2_s(q, k, cos); 
                ff[1] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * F2_s(q, k, cos) * S2_s(q, k, cos);  
                ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * (plin_kmq - plin_q);  
                ff[3] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * pow(S2_s(q, k, cos), 2.);  
                ff[4] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * S2_s(q,k,cos); 
                ff[5] = 1./pow(2.*M_PI,3.) * 4. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_k * plin_q   * S2_s(q, k, cos) * F2(q, k, -cos);  
            }
            else if(hm_switch == MATTER)
            {                
                ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * pow(F2_s(q, k, cos), 2.);
                ff[1] = 1./pow(2.*M_PI,3.) * 6. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_k * plin_q * F3_s(k, q, cos);
            }
            for(int i =0;i<nn;i++){
              if (isnan(ff[i]))
                 printf("%d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",i, k, q, kmq, plin_q, plin_kmq, ff[i]);
            }

        }      
        else if(IR_switch == WIR){
            double plin_IR_q   = pm_IR_LO(Cx, q, z, SPLIT);
            double plin_IR_kmq = pm_IR_LO(Cx, kmq, z, SPLIT);
            // double plin_IR_kpq = pm_IR_LO(Cx, kpq, z, SPLIT);

            if(hm_switch == HALO)
            {
                ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * F2_s(q, k, cos); 
                ff[1] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * F2_s(q, k, cos) * S2_s(q, k, cos);  
                // ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * (plin_IR_kmq - plin_IR_q);  
                ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * (plin_IR_kmq);  
                ff[3] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow(S2_s(q, k, cos), 2.);  
                ff[4] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * S2_s(q, k, cos); 
                ff[5] = 1./pow(2.*M_PI,3.) * 4. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q * S2_s(q, k, cos) * F2(q, k, -cos);  
            }
            else if(hm_switch == MATTER)
            {            
                // ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow(F2_s(q, k, cos), 2.) * theta_kmq\
                //       + 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kpq * pow(F2_s(q, k, -cos), 2.) * theta_kpq;
                ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow(F2_s(q, k, cos), 2.);
                ff[1] = 1./pow(2.*M_PI,3.) * 6. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q * F3_s(k, q, cos);
            }


            for(int i =0;i<nn;i++){
              if (isnan(ff[i]))
                 printf("%d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",i, k, q, kmq, plin_IR_q, plin_IR_kmq, ff[i]);

            } 
        }
    }           
    else{
      for(int i =0;i <nn; i++)
        ff[i] = 0.;
    }

      return 0;
}


/**
 * Compute the loop contributions dure to nonlinear evolution of matter fluctuations and nonlinear halo bias, rising from non-Gaussian initial conditions of local shape
 * The function PNG_loop_integrands() defines the integrand and Compute_PNG_loops() computes the integrals 
 * 
 * @param Cx           Input: pointer to cosmology structure 
 * @param k            Input: wavenumber 
 * @param z            Input: redshift of interest
 * @param IR_switch    Input: switch to decide whether to perform IR resummation or no
 * @param SPLIT        Input: switch to set the method to perform the wiggle-nowiggle split of matter power spectrum
 * @param result       Output: an output array containing the results of the 1loop terms, has 8 elements for hm_switch=HALO
 * @return void
 */

void Compute_PNG_loops(struct Cosmology *Cx, double k, double z,  long IR_switch, long SPLIT, double *result)
{
      struct integrand_parameters2 par;

      int    ncomp = 8;                               
      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];

      par.p1 = Cx;
      par.p4 = k;
      par.p5 = z;
      par.p6 = log(1.e-4);
      par.p7 = log(100.); 
      par.p13 = IR_switch;
      par.p14 = SPLIT;

      int ndim = 2, nvec = 1, verbose = 0, last = 4, seed = 0, 
          mineval = 0, maxeval = 2e6, key1 =50, key2 = 50, key3 = 1, maxpass = 200, 
            border =0, maxchisq = 10, mindeviation=0.25,
            ngiven = 0, ldxgiven = ndim, nextra = 0; 
      int nregions, neval;   

      Divonne(ndim, ncomp, PNG_loop_integrands, &par, nvec,
                 RelErr, AbsErr, verbose, seed,
                 mineval, maxeval,  key1, key2, key3, maxpass,
                 border, maxchisq, mindeviation,
                 ngiven, ldxgiven, NULL, nextra, NULL,
                 NULL, NULL,
                 &nregions, &neval, fail, result, error, prob);
             
      // int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 2.e8, key = 13;
      // int nregions, neval;
      
      // Cuhre(ndim,ncomp, G_loop_integrands, &par, nvec,
      //          RelErr, AbsErr, verbose | last,
      //          mineval, maxeval, key,
      //          NULL, NULL,
      //          &nregions, &neval, fail, result, error, prob);
      
      return;
}

static int PNG_loop_integrands(const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)    
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);      

      struct Cosmology *Cx = pij.p1;
      double k            = pij.p4;
      double z            = pij.p5;
      double logqmin      = pij.p6;
      double logqmax      = pij.p7;
      long IR_switch      = pij.p13;
      long SPLIT          = pij.p14;
      static int cleanup  = 0;
       
      double logq     = (x[0] * (logqmax - logqmin) + logqmin); 
      double cos      = (2. * x[1] - 1.); 
      double q        = exp(logq); 
      double kmq      = pow(fabs(pow(q, 2.) + pow(k, 2.) - 2. * q * k * cos), 1./2.);
      double A        = (q * k * cos - pow(q, 2.))/pow(kmq, 2.);
      double tk_q_inv = 1./Mk_dlnMk(Cx, q, z, TRANS);
      double tk_kmq = Mk_dlnMk(Cx, kmq, z, TRANS);
      double tk_q   = Mk_dlnMk(Cx, q, z, TRANS);
      
      /// Factor of 2. * (logqmax - logqmin) is due to change of variable from 0 to logarithmic k, and a factor of 2*PI is due to integration over azimuthal angle. 
      /// Note that to compare the theoretical predictions against Emiliano's measurement, since he is using a different notation for Fourier transform, I need to 
      /// devide each 0 power spectrum by a factor of 1/pow(2.*M_PI,3.), which I do in my pk_lin() function. If using another notation for Fourier transform 
      /// (the one that I usually use, which has a factor of 1/pow(2*M_PI,3) in the definition), you need to multiply these integrands by a factor of 1/pow(2*M_PI,3).
      
      /// The integrands below correspond to the follwing bias combinaions: 
      /*  ff[0] : b1 b_zeta
          ff[1] : b1 b_zetadelta
          ff[2] : b2 b_zeta
          ff[3] : bG2 b_zeta  (I-term)
          ff[4] : b2 b_zetadelta
          ff[5] : bG2 b_zetadelta
          ff[6] : b3nl b_zeta,  b3nl = 2*(bG2+2/5*btd) 
	  ff[7] : fnl^2 term from matter bispectrum
      */

    if(kmq>=exp(logqmin)){
  
      if(IR_switch == NOIR){
          double plin_k   = Pk_dlnPk(Cx, k, z, LPOWER);
          double plin_q   = Pk_dlnPk(Cx, q, z, LPOWER);
          double plin_kmq = Pk_dlnPk(Cx, kmq, z, LPOWER);

          ff[0] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq  * F2_s(q, k, cos) * A * tk_q_inv; 
          ff[1] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq  * F2_s(q, k, cos) * tk_q_inv;  
          ff[2] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * (plin_kmq * A + plin_q) * tk_q_inv;  
          ff[3] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq  * S2_s(q, k, cos) * A * tk_q_inv;  
          ff[4] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * (plin_kmq - plin_q) * tk_q_inv;  
          ff[5] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq  * S2_s(q, k, cos) * tk_q_inv;  
          ff[6] = 4. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * S2_s(q, k, cos) * F2(q, k, -cos);
	        ff[7] = 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * F2_s(q,k,cos) * (plin_q * tk_kmq/tk_q + plin_kmq * tk_q/tk_kmq);  
      }          
      else if(IR_switch == WIR){
          double plin_IR_k   = pm_IR_LO(Cx, k, z, SPLIT);
          double plin_IR_q   = pm_IR_LO(Cx, q, z, SPLIT);
          double plin_IR_kmq = pm_IR_LO(Cx, kmq, z, SPLIT);

          ff[0] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq  * F2_s(q, k, cos) * A * tk_q_inv; 
          ff[1] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq  * F2_s(q, k, cos) * tk_q_inv;  
          ff[2] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * (plin_IR_kmq * A + plin_IR_q) * tk_q_inv;  
          ff[3] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq  * S2_s(q, k, cos) * A * tk_q_inv;  
          ff[4] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * (plin_IR_kmq - plin_IR_q) * tk_q_inv;  
          ff[5] = 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq  * S2_s(q, k, cos) * tk_q_inv;  
          ff[6] = 4. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * S2_s(q, k, cos) * F2(q, k, -cos);  
	        ff[7] = 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * F2_s(q,k,cos) * (plin_IR_q * tk_kmq/tk_q + plin_IR_kmq * tk_q/tk_kmq);
	       
      }
    }
    else{
      for(int i =0; i<8; i++)
        ff[i] = 0.;
    }  

      // for(int i=0;i<8;i++)
      //   printf("integ NG %d %12.6e %12.6e %12.6e %12.6e\n",i, k, q, kmq, ff[i] );
   
      return 0;
}    





//*********************************************************************************////
//****** Second-order kernels *******//// 
//*********************************************************************************////

////This is F2(k1,k2-k1) since only this configuration apears in the loops ////
double F2_s(double k1,double k2,double mu)
{
      double f = 0.;
      // if(k2<0.01)
	     //   f = 3.*pow(k2/k1,2.) - 5./7. * pow(k2/k1,2.)* pow(mu,2.) + 13./14.* pow(k2/k1,3.) * mu - 10./7. * pow(k2/k1,3.)*pow(mu,3.);
      // else		
	       f = pow(k2,2.)*(3.*k1 + 7.*k2*mu - 10.*k1*pow(mu,2.))/(14.*k1*(pow(k1,2.) + pow(k2,2.) - 2.*k1*k2*mu));
      
      return f;
}


////This is S2(k1,k2-k1) since only this configuration apears in the loops ////
double S2_s(double k1,double k2,double mu)
{     
      double f = 0.;
      f = pow(-k1 + k2*mu,2.)/(pow(k1,2.) + pow(k2,2.) - 2.*k1*k2*mu) - 1.;
      
      return f;   
}

////This is the symmetrized F3(q,-q,kmq) given in 1603.04405 
double F3_s(double k,double q, double mu)
{

    double f = 0.;
    double kmq2   = fabs(pow(q,2.)+pow(k,2.)-2.*q*k*mu);
    double kdotq = k*q*mu;

    f = 1./kmq2*(5.*pow(k,2.)/63. - (11.*kdotq)/54. - pow(k,2.)*pow(kdotq,2.)/(6.*pow(q,4.))\
       +  19.*pow(kdotq,3.)/(63.*pow(q,4.))  - 23.*pow(k,2.)*kdotq/(378.*pow(q,2.))\
       - 23.*pow(kdotq,2.)/(378.*pow(q,2.)) + pow(kdotq,3.)/(9.*pow(k*q,2.)));

    return f;
}

////This is S2(k1,k2) ////
double S2(double mu)
{     
      double f = 0.;
      f = pow(mu,2.) - 1.;

      return f;   
}

double F2(double k1,double k2,double mu)
{
      double f=0.;
      
      f = 5./7. + 1./2. * (k1/k2 + k2/k1) * mu + 2./7. * pow(mu,2.);

      return f;
}

