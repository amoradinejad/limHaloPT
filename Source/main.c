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

	// int nlines     = 7;
      // int lines[7]   = {CO10, CO21, CO32, CO43, CO54, CO65, CII};
      // int JJ[7]      = {1,2,3,4,5,6,0};
      // size_t ninterp = 150;

	// int nlines     = 4;
 //      int lines[4]   = {CO21, CO32, CO43, CII};
 //      int JJ[4]      = {2, 3, 4, 0};
 //      size_t ninterp = 150;



	int nlines     = 1;
      int lines[1]   = {CII};
      int JJ[1]      = {0};
      size_t ninterp = 150;


	double gb_pars[] = {gb.logAs,gb.ns,gb.h,gb.Omega_b,gb.Omega_cdm,gb.sigFOG0};

	double min_pars[gb.Npars];
	double max_pars[gb.Npars];

	struct Cosmology Cx_ref;
	struct Cosmology Cx_min[gb.Npars];
      struct Cosmology Cx_max[gb.Npars]; 

	for(i=0;i<gb.Npars;i++){ 
		min_pars[i] = (1.-inc)*gb_pars[i];
		max_pars[i] = (1.+inc)*gb_pars[i];
	}
	

	for(i=0;i<gb.Npars;i++){
		Cx_ref.cosmo_pars[i] = gb_pars[i];
	}

	for(i=0;i<gb.Npars;i++){
		for(j=0;j<gb.Npars;j++){
			if(i==j){
				Cx_min[i].cosmo_pars[j] = min_pars[i];
				Cx_max[i].cosmo_pars[j] = max_pars[i];
			}
			else if (i!=j){	
				Cx_min[i].cosmo_pars[j] = gb_pars[j];
				Cx_max[i].cosmo_pars[j] = gb_pars[j];
			}
			//printf("i=%ld \t  j=%ld \t min = %12.6e \t max= %12.6e \n",i,j,Cx_min[i].cosmo_pars[j],Cx_max[i].cosmo_pars[j]);
		}
		//printf("\n \n");
	}

      // double scatter    = 0.37;
      // double scatter_1h = p_sig_shot(scatter);
      // double scatter_2h = p_sig_Tbar(scatter);
      // printf("%12.5e %12.5e \n",scatter_2h,scatter_1h);

	///**--------------------------*///
	///* compare the theoretical and measured HMFs
 	///**--------------------------*///

	// FILE *files[1];	
	// double *M, *dndM_PSC,*dndM_ST,*dndM_TR;
	// double n=1000;
	// long z[] = {1} ;
	// int  nz  = 1;

	// for(j=0;j<nz;j++){
	// 	char filename[20];
	// 	sprintf(filename,"hmf-th-z%ld.00.dat",z[j]);
	// 	files[j] = fopen(filename,"w+");

	// 	M = loginit_1Darray(n,pow(10,1),pow(10,15));   ///in unit of M_sun
	// 	dndM_PSC = make_1Darray(n);
	// 	dndM_ST = make_1Darray(n);
	// 	dndM_TR= make_1Darray(n);
	// 	for(i=0;i<n;i++){
	// 		dndM_PSC[i] = mass_func(&Cx_ref, M[i], z[j], PSC);
	// 		dndM_ST[i]  = mass_func(&Cx_ref, M[i], z[j], ST);
	// 		dndM_TR[i]  = mass_func(&Cx_ref, M[i], z[j], TR);
	// 		fprintf(files[j], "%12.6e %12.6e %12.6e %12.6e \n",M[i], M[i] * dndM_PSC[i],M[i] * dndM_ST[i],M[i] * dndM_TR[i]); 

	// 	}	
	// 	free(dndM_PSC);
	// 	free(dndM_ST);
	// 	free(dndM_TR);
	// 	free(M);

	// 	fclose(files[j]);
	// }



      // int nk          = 50;
	// double *k       = loginit_1Darray(nk, 1.e-3, 1.);
	// double z[2]     = {1.,4.5}; 
	// int nz	    = 2;

	// for(j=0;j<nz;j++){
	// 	for(i=0;i<nk;i++){
	// 		printf("%12.6e %12.6e %12.6e %12.6e %12.6e  \n", z[j], k[i], Pk_dlnPk_HV(&Cx_ref, k[i],z[j], LPOWER), Pk_dlnPk(&Cx_ref, k[i], z[j], LPOWER), pow(growth_D(&Cx_ref, z[j]),2.) *pk_Gfilter_nw(&Cx_ref, k[i], k[0]));
	// 	}
	// }

	// printf("%12.6e %12.6e \n" , mass_func_sims(&Cx_ref, 1.e10, 2., ST), mass_func(&Cx_ref, 1.e10, 2., ST)); 
	// printf("%12.6e \n",Hubble(&Cx_ref,2));	
 
	
	clock_t tic_r = clock();
	Cosmology_init(&Cx_ref, pk_kmax, pk_zmax, nlines, lines, ninterp, M_min, mode_mf);
	printf("Reference Cosmology initialized\n");
	clock_t toc_r = clock();
	printf("Elapsed: %f seconds for ref cosmology\n", (double)(toc_r - tic_r)/CLOCKS_PER_SEC);


	// printf("z = 2: %12.6e %12.6e \n", PS_shot(&Cx_ref, 2.,1)*pow(gb.h,3.), PS_shot(&Cx_ref, 0, 2.)*pow(gb.h,3.));
	// printf("z = 4: %12.6e %12.6e \n", PS_shot(&Cx_ref, 4.,1)*pow(gb.h,3.), PS_shot(&Cx_ref, 0, 4.)*pow(gb.h,3.));
	// printf("z = 6: %12.6e %12.6e \n", PS_shot(&Cx_ref, 6.,1)*pow(gb.h,3.), PS_shot(&Cx_ref, 0, 6.)*pow(gb.h,3.));

	// double shot = PS_shot(&Cx_ref, 2., 0);
	// double clust_pt = 0.;

	// int nk     = 400;
	// double *k  = loginit_1Darray(nk, 1.e-3, 20.);
	// for(i=0;i<nk;i++){
	// 	clust_pt = PS_line(&Cx_ref, k[i], 2, 0);
	// 	printf("%12.6e %12.6e %12.6e\n", k[i], shot, clust_pt);
	// }


	/**
	 * compute multipoles of power spectrum 
 	 */
	// double mono = 0., quad = 0., hexa = 0.;
	// int nk = 400, nz = 3, nl = 2;
	// double *k  = loginit_1Darray(nk, 1.e-3, 20.);
      // double z[] = {2.,4.};

	// FILE *fp1;
      // char filename1[FILENAME_MAX];

	// for(i=0;i<nl;i++){
	// 	for(j=0;j<nz;j++)
	// 		for(l=0;l<nl;l++){
	// 			sprintf(filename1,"%s/line/ps/RSD_pk_pt_J%d_z%d.txt", gb.output_dir, JJ[l],(int)z[j]);
	// 			fp1    = fopen(filename1, "ab");
	// 			mono = ps_line_multipoles(&Cx_ref, &Cx_ref, k[i], z[j], lines[l], 0);
	// 			quad = ps_line_multipoles(&Cx_ref, &Cx_ref, k[i], z[j], lines[l], 2);
	// 			hexa = ps_line_multipoles(&Cx_ref, &Cx_ref, k[i], z[j], lines[l], 4);
	// 			fprintf(fp1,"%12.6e %12.6e %12.6e %12.6e \n", k[i], mono, quad, hexa);
	// 		}
	// }
	// fclose(fp1);
	// free(k);
	




	/**
	 * compare my class matter power spectrum 
	 * with the input of hiddenvalley sims
 	 */
	
 	//CL_Cosmology_initilize(&Cx_ref,pk_kmax,pk_zmax);

 	// int nk           = 3000;
	// double *k        = loginit_1Darray(nk, 1.e-5, 400.);
	// double *pk_class = make_1Darray(nk);
	// double *pk_HV    = make_1Darray(nk);
	// double *pk_nw    = make_1Darray(nk);

 	// for(i=0;i<nk;i++){
 	// 	 pk_class[i] = Pk_dlnPk(&Cx_ref, k[i], 0., LPOWER);
 	// 	 pk_HV[i]    = Pk_dlnPk_HV(&Cx_ref, k[i], 0., LPOWER);
 	// 	 pk_nw[i]    = pm_nowiggle(&Cx_ref, k[i], 0., k[0], 0, GFILTER);
 	// 	printf("%12.6e %12.6e %12.6e %12.6e\n",k[i]/gb.h, pk_nw[i]*pow(gb.h,3.), pk_class[i]*pow(gb.h,3.), pk_HV[i]*pow(gb.h,3.));
 	//}
      // free(pk_class);
      // free(pk_HV);
      // free(pk_nw);
      // free(k);

	// double M[50], L[50];
	// FILE  *ifp = fopen("/Volumes/Data/Documents/Git/LIM_PS_HM/Output/line/lum/L_Line_ModelCII_A.dat","r");

 	// for(i=0L;i<50;i++)
 	//      	fscanf(ifp,"%lg %lg \n", &M[i],&L[i]);

 	// for(i=0L;i<50;i++)
 	// 		printf("%12.6e %12.6e \n",log10(M[i]/gb.h),logSFR_Behroozi(log10(M[i]/gb.h), 2., 0));
	
	// for(i=0L;i<50;i++)
	// 	printf("%12.9e %12.9e \n",log10(M[i]/gb.h),luminosity(M[i]/gb.h, 2., CII));


	///**--------------------------*///


	/**
  	 * compare the theoretical and measured HMFs
 	 */

	// FILE *files[3];	
	// double *M, *dndM_PSC,*dndM_ST,*dndM_TR;
	// double n=1000;
	// long z[] ={2,4,6} ;

	// for(j=0;j<3;j++){
	// 	char filename[20];
	// 	sprintf(filename,"hmf-th-z%ld.00.dat",z[j]);
	// 	files[j] = fopen(filename,"w+");

	// 	M = loginit_1Darray(n,pow(10,1),pow(10,15));   ///in unit of M_sun
	// 	dndM_PSC = make_1Darray(n);
	// 	dndM_ST = make_1Darray(n);
	// 	dndM_TR= make_1Darray(n);
	// 	for(i=0;i<n;i++){
	// 		dndM_PSC[i] = mass_func(&Cx_ref, M[i], z[j], PSC);
	// 		dndM_ST[i]  = mass_func(&Cx_ref, M[i], z[j], ST);
	// 		dndM_TR[i]  = mass_func(&Cx_ref, M[i], z[j], TR);
	// 		fprintf(files[j], "%12.6e %12.6e %12.6e %12.6e \n",M[i], M[i] * dndM_PSC[i],M[i] * dndM_ST[i],M[i] * dndM_TR[i]); 

	// 	}	
	// 	free(dndM_PSC);
	// 	free(dndM_ST);
	// 	free(dndM_TR);
	// 	free(M);

	// 	fclose(files[j]);
	// }


	/**
	 * Compute the halo biases as a function of mass for several redshift using Sheth-Tormen mass function
 	 */
	
	// double b1, b2, bg2, btd;
	// int nz = 4, nM=50;
 	// double *MM = loginit_1Darray(nM,1.e9,1.e16);
	// double z[4] = {0.,2.,4.,6.};
	// double hbias_arr[4];

 	// FILE *fp1;
 	// char filename1[FILENAME_MAX];

 	// for(i=0;i<nz;i++){
 	// 	sprintf(filename1,"%s/bias_halo_z%d.txt", gb.output_dir, (int)z[i]);
 	// 	fp1 = fopen(filename1, "ab");
 	// 	for(j=0;j<nM;j++){
	// 		halo_bias(&Cx_ref, MM[j], z[i], ST, hbias_arr);
	// 		b1  = hbias_arr[0];
	// 		b2  = hbias_arr[1];
	// 		bg2 = hbias_arr[2];
 	//          btd = hbias_arr[3];  
	// 	    	fprintf(fp1,"%12.6e %12.6e %12.6e %12.6e %12.6e \n",MM[j],b1,b2,bg2,btd);
	// 	}
	//  	fclose(fp1);
	// }      



 	/**
       * Compute the line biases as a function of redshift with and without the nfw profile	
       */
	
	// int nz = 50;
 	// double *zz = init_1Darray(nz,0.0,11.);
 	// double b1_l[nlines][nz], b2_l[nlines][nz];

      // FILE *fp2;
      // char filename2[FILENAME_MAX];
      // sprintf(filename2,"%s/bias_line.txt", gb.output_dir);
      // fp2 = fopen(filename2, "ab");
	// double *lbias_arr = make_1Darray(2);

	// double *result = make_1Darray(3);
 	// double b1, b2, bg2, btd;

  	//fpr(j=0;j<nlines;j++){
	 	// for(i=0;i<nz;i++){
		  	// line_bias(Cx_ref.Lines[j], zz[i], lbias_arr);
		  	// b1_l[j][i]  = lbias_arr[0];
		  	// b2_l[j][i]  = lbias_arr[1];
		     	//printf("%12.6e %12.6e %12.6e \n",zz[i],b1_l[j][i],b2_l[j][i]);
		// }		
		// free(lbias_arr);
  	//}

	// for(i=0;i<nz;i++)
	//       fprintf(fp2," %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",zz[i],b1_l[0][i],b2_l[0][i],b1_l[1][i],b2_l[1][i],b1_l[2][i],b2_l[2][i],b1_l[3][i],b2_l[3][i],b1_l[4][i],b2_l[4][i],b1_l[5][i],b2_l[5][i],b1_l[6][i],b2_l[6][i]);

	// fclose(fp2);
	// free(zz);
	

	/**
	 * Compute the line biases, accounting for nfw profile 
	 * for a fixed redshift as a function of k
	 */
	// int nk          = 400;
      // int nz 	    = 3;
	// double *k       = loginit_1Darray(nk, 1.e-3, 20.);
	// double z[3]     = {2.,4.,6.}; 

      // double result1[3], result2[3];
	// double mom1_z2, mom1_z4, b1_z2, b2_z2, bg2_z2, btd_z2, b1_z4, b2_z4, bg2_z4, btd_z4;
	// mom1_z2 = mass_moment1(&Cx_ref, z[0], 1.e9, ST, CII);
	// mom1_z4 = mass_moment1(&Cx_ref, z[1], 1.e9, ST, CII);

	// for(l=0;l<nk;l++){
	// 	HM_1h2h(&Cx_ref, k[0], z[0], 1.e9, ST, CII, LINE, result1);
	// 	HM_1h2h(&Cx_ref, k[0], z[1], 1.e9, ST, CII, LINE, result2);

	// 	b1_z2  = result1[1]/mom1_z2;
	// 	b2_z2  = result1[2]/mom1_z2;
	// 	bg2_z2 =  -2./7.  * (b1_z2-1.);
      // 	btd_z2 =  23./42. * (b1_z2-1.); 

      // 	b1_z4  = result2[1]/mom1_z4;
	// 	b2_z4  = result2[2]/mom1_z4;
	// 	bg2_z4 =  -2./7.  * (b1_z4-1.);
      // 	btd_z4 =  23./42. * (b1_z4-1.); 

	// 	printf("%12.6e %12.6e %12.6e %12.6e %12.6e \n", k[0], b1_z2,b2_z2,bg2_z2,btd_z2);
	// 	printf("%12.6e %12.6e %12.6e %12.6e %12.6e \n", k[1], b1_z4,b2_z4,bg2_z4,btd_z4);


 	// }	
 	// free(k);


	/* 
	 * Set the k and z arrays
	 */

      int nk = 200;
      double *k  = loginit_1Darray(nk, 1.e-3, 8.);
      
      // int nz;
      // double *z;
      // for(i=0;i<nlines;i++){
      //       if(JJ[i] == 0){
      //             nz   = 2;
      //             z    = make_1Darray(nz);
      //             z[0] = 1.;
      //             z[1] = 4.5;
      //       }
      //       else if(JJ[i] == 3){
      //             nz   = 2;
      //             z    = make_1Darray(nz);
      //             z[0] = 1.;
      //             z[1] = 2.;
      //       }
      //       else if(JJ[i] == 2 || JJ[i] == 4){
      //             nz   = 1;
      //             z    = make_1Darray(nz);
      //             z[0] = 1.;
      //       }
      // }


      /**
	 * Compute the line clustering power spectrum
	 * PS_line: PT, PS_HM: HM
	 */

 	// double b2, Tbar, b22_shot, lbias_arr[2];
	// for(j=0;j<nz;j++){	
	// 	Tbar   = Tbar_line(&Cx_ref, 0, z[j]); 
	// 	line_bias(Cx_ref.Lines[0],z[j],lbias_arr);
	// 	b2 = lbias_arr[1];
	// 	b22_shot = 0.25 * pow(b2,2.) * pow(Tbar,2.) * b22_ls(&Cx_ref,z[j]);

	// 	printf("b22_ls %12.6e %12.6e %12.6e %12.6e %12.6e\n",z[j], Tbar,  b2, b22_ls(&Cx_ref,z[j]), b22_shot);
	// }

 //      int nk = 200;
 //      double *k  = loginit_1Darray(nk, 1.e-3, 5.);

 // 	double clust_hm;
 // 	int line_id =0;

 // 	for(i=0;i<nlines;i++){
 // 		line_id = i;
 // 		for(j=0;j<nz;j++){
	// 		for(l=0;l<nk;l++){
	// 			clust_hm = PS_line_HM(&Cx_ref, k[l], z[j], M_min, mode_mf, lines[i], line_id);
	// 			printf("%d %d %d %12.6e %12.6e %12.6e \n", i, j, l, z[j], k[l]/gb.h, pow(gb.h,3.)*clust_hm);
	// 		}
	// 	}	
	// }

	
	
	double shot_p   = 0.;
	double Omegam   = (Cx_ref.cosmo_pars[3] + Cx_ref.cosmo_pars[4]);
	double rhom_bar = Omegam * rhoc(&Cx_ref,0.); 

	double lbias_arr[2];
	double input[5];

	double b1   = 0.;
	double b2   = 0.;
	double Tbar = 0.;

	///uncomment if saving Tbar and bias to a file.
	// FILE *fp1;
 // 	char filename1[FILENAME_MAX];
 	int nz = 50;
 	double *z = init_1Darray(nz,0.0,11.);

 	for(i=0;i<nlines;i++){
 		for(j=0;j<nz;j++){
 			// sprintf(filename1,"%s/line/Tbar_J%d_z%d.txt", gb.output_dir, JJ[i],(int)z[j]);
	 		// fp1    = fopen(filename1, "ab");

	 		Tbar = Tbar_line(&Cx_ref, i, z[j]);
			line_bias(Cx_ref.Lines[i],z[j],lbias_arr);
			b1 = lbias_arr[0];
			b2 = lbias_arr[1];
			// fprintf(fp1,"%d %12.6e %12.6e %12.6e %12.6e \n",i, z[j], b1, b2, Tbar); 
			// printf("%d %12.6e %12.6e %12.6e %12.6e \n",i, z[j], b1, b2, Tbar); 

			input[0] = Tbar_line(&Cx_ref, i,z[j]); 
			input[1] = b1;
			input[2] = 0.25 * pow(b2,2.) * pow(input[0],2.) * b22_ls(&Cx_ref,z[j]);
			input[3] = PS_shot(&Cx_ref, z[j], i);
			input[4] = rhom_bar; 

			//If accounting for the nfw profile, the shot noise will be scale-dependent. So loop over k-values
			// for(l=0;l<191;l++){
			// 	PS_shot_HM(&Cx_ref, k[l], z[j], M_min, input, mode_mf, lines[i]);
			// }

			// When neglecting the nfw profile for comparison with sims, the shot noise is scale-independent. 
			// So uncomment the following line and comment out the above for loop to just compute it for a single k
			PS_shot_HM(&Cx_ref, k[0], z[j], M_min, input, mode_mf, lines[i]);
			
		}	
	}	
	// free(k);
 	free(z);


	/**
	 * Compute the 1loop matter power spectrum. The p13, p22, ptot are outputed inside pm_IR_NLO function 
	 */

      // double * ps_hloops = make_1Darray(6);
      // double pm_lin_IR,pm_1loop_IR;

      // int nk 	    = 100;
	// double *k       = loginit_1Darray(nk, 1.e-3, 20.);
	// double z[2]     = {4.,2.}; 
	// int nz	    = 2;

	// for(j=0;j<nz;j++){
      //     for(i=0;i<nk;i++){
      //      	// pm_lin_IR   = pm_IR_LO(&Cx_ref, k[i], z, GFILTER);
      //      	// pm_1loop_IR = pm_IR_NLO(&Cx_ref, k[i], z, GFILTER);
      //      	// printf("%12.6e %12.5e %12.6e %12.7e \n", k[i],Pk_dlnPk(&Cx_ref, k[i], z, LPOWER),pm_lin_IR,pm_1loop_IR);

      //      	PS_hh_G(&Cx_ref, k[i], z[j], 1.e12, LOOP, WIR, GFILTER, ST);
      //      }
      // }	
      // free(ps_hloops);
	// free(k);


	/**
	 * Test the halo model prediction of matter power spectrum, with and without 1loop
	 */
	// long J 	    = 1;
	// int nk 	    = 100;
	// double z     = 1.;
	// double Mmin  = 1.e9;
	// long mode_mf = ST, mode_lum = CO10;
	// double *k    = loginit_1Darray(nk, 1.e-3, 20.);

	// double cs2  = 2.;
 	// double khat = 1.; //In unit of h/Mpc 
      // double Ms   = 1.e6;  //Lower mass limit for the halo model integration
      // double matter[2], hm_corrs[2];
      // double nfw, p1h_c, p2h_c;
      // double plin, pm_ct, pm_lin_IR, pm_1loop_IR;
      // double p1h, p2h, p2h_1loop, p2h_ct;

      // double Omegam = (Cx_ref.cosmo_pars[3] + Cx_ref.cosmo_pars[4]);
  	// double rhom_bar = Omegam * rhoc(0,&Cx_ref); 
      // mhmc(z,&Cx_ref,mode_mf,hm_corrs);

	// for(i=0;i<nk;i++){
	// 	nfw   = nfw_profile(k[i], Ms, z, &Cx_ref);
	// 	p1h_c = pow(nfw,2.)*Ms/rhom_bar * hm_corrs[0];
 	//    p2h_c = nfw * hm_corrs[1];
	// 	HM_1h2h(k[i], z, Mmin, &Cx_ref, mode_mf, mode_lum, MATTER, matter);  /// in unit of M_sun/Mpc^3
		
	// 	plin        = Pk_dlnPk(k[i], z, &Cx_ref, LPOWER) * pow(gb.h,3.);
	// 	pm_lin_IR   = pm_IR_LO(k[i], z, k[0], &Cx_ref, GFILTER) * pow(gb.h,3.);
      //    pm_1loop_IR = pm_IR_NLO(k[i], z, k[0], &Cx_ref, GFILTER) * pow(gb.h,3.) ;
      //    pm_ct       = - 2. * cs2 * pow(k[i]/gb.h, 2.)/(1.+pow((k[i]/gb.h)/khat,2.)) * pm_lin_IR;  //Eq. 2.22 of 2004.09515, Pade approximation to make p_ct finite at large k

	// 	p1h         = (pow(rhom_bar,-2.) * matter[0] + p1h_c) * pow(gb.h,3.);
	// 	p2h         = pow(matter[1]/rhom_bar + p2h_c, 2.) * plin;
	// 	p2h_1loop   = pow(matter[1]/rhom_bar + p2h_c, 2.) * pm_1loop_IR ;
	// 	p2h_ct      = pow(matter[1]/rhom_bar + p2h_c, 2.) * pm_ct;

	// 	printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n", k[i]/gb.h, p1h, p2h, p2h_1loop, p2h_ct, pm_1loop_IR, pm_ct, plin);
	// }	
	//free(k);


	cleanup(&Cx_ref);
	
	return 0;
}


void initialize()
{
	extern struct globals gb;



//printf("\n*************************************************************************\n");
	getcwd(gb.project_home , sizeof(gb.project_home));  //This gives the path to the source directory
	chdir("..");  ///Change the path to the parent directory
	getcwd(gb.project_home , sizeof(gb.project_home));

	sprintf(gb.output_dir,"%s/Output", gb.project_home);
	sprintf(gb.data_dir,"%s/Input", gb.project_home);
	sprintf(gb.data_priors,"%s/Priors", gb.data_dir);

	sprintf(gb.SFR_filename,"%s/release-sfh_z0_z8_052913/sfr/sfr_release.dat",gb.data_dir);  
	sprintf(gb.Planck_Fisher_filename,"%s/Fisher_Planck_LCDM_marginalized.dat",gb.data_priors);  

	logSFR_alloc_init();

	// //printf("*************************************************************************\n\n");

	///// consistent with base plikHM TTTEEE lowTEB lensing post BAO H070p6 JLA”
	/////https://wiki.cosmos.esa.int/planckpla2015/images/f/f7/Baseline_params_table_2015_limit68.pdf/////
	gb.c 		 = 2.99792458e5;  /// In units of km/s
	gb.h 		 = 0.677;
	gb.Omega_cdm = 0.11923/pow(gb.h,2.);  ////omega_cdm = Omega_cdm h^2 ;
	gb.Omega_b 	 = 0.02247/pow(gb.h,2.);   ///omega_b = Omega_b h^2;
	gb.Omega_r 	 = 0.0000910958;  ////radiation = photons + neutrinos
	gb.Omega_g 	 = 5.3956421715871286e-05;   ////photons, input for Class
	gb.Omega_nu  = 0.00;    ////neutrinos
	gb.ns        = 0.96824;
	// gb.As 	 = 2.10732e-9;
	gb.As 	 = 2.1085e-9;
	gb.logAs 	 = log(gb.As*pow(10.,10.));  ///3.0665
	gb.kp 	 = 0.05;  //in unit of Mpc^-1

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
	
  
	

	///*************************************************************************************/////
 			// Functions to put in main() for DEBUGGING purposes 
	///*************************************************************************************/////

	// CL_Cosmology_initilize(&Cx_ref,pk_kmax,pk_zmax);



	///**--------------------------*///
	///* compare my class matter power spectrum 
	///* with the input of hiddenvalley sims
 	///**--------------------------*///
	
  	// CL_Cosmology_initilize(&Cx_ref,pk_kmax,pk_zmax);

 	// int nk     = 3000;
	// double *k  = loginit_1Darray(nk, 4.e-5, 120.);
	// double *pk = make_1Darray(nk);
 	// for(i=0;i<nk;i++){
 	//  	pk[i] = Pk_dlnPk(&Cx_ref, k[i], 0., LPOWER);
 	// 	printf("%12.6e %12.6e \n",k[i]/gb.h,pk[i]*pow(gb.h,3.));
 	// }



	///**--------------------------*///
	///* compare the theoretical and measured HMFs
 	///**--------------------------*///

	// FILE *files[3];	
	// double *M, *dndM_PSC,*dndM_ST,*dndM_TR;
	// double n=1000;
	// long z[] ={1} ;

	// for(j=0;j<3;j++){
	// 	char filename[20];
	// 	sprintf(filename,"hmf-th-z%ld.00.dat",z[j]);
	// 	files[j] = fopen(filename,"w+");

	// 	M = loginit_1Darray(n,pow(10,1),pow(10,15));   ///in unit of M_sun
	// 	dndM_PSC = make_1Darray(n);
	// 	dndM_ST = make_1Darray(n);
	// 	dndM_TR= make_1Darray(n);
	// 	for(i=0;i<n;i++){
	// 		dndM_PSC[i] = mass_func(&Cx_ref, M[i], z[j], PSC);
	// 		dndM_ST[i]  = mass_func(&Cx_ref, M[i], z[j], ST);
	// 		dndM_TR[i]  = mass_func(&Cx_ref, M[i], z[j], TR);
	// 		fprintf(files[j], "%12.6e %12.6e %12.6e %12.6e \n",M[i], M[i] * dndM_PSC[i],M[i] * dndM_ST[i],M[i] * dndM_TR[i]); 

	// 	}	
	// 	free(dndM_PSC);
	// 	free(dndM_ST);
	// 	free(dndM_TR);
	// 	free(M);

	// 	fclose(files[j]);
	// }


	///**--------------------------*///
	///* Compute the halo biases as a function of mass for several redshift using Sheth-Tormen mass function*///
 	///**--------------------------*///
	
	// double b1, b2, bg2, btd;
	// int nz = 4, nM=50;
 	// double *MM = loginit_1Darray(nM,1.e9,1.e16);
	// double z[4] = {0.,2.,4.,6.};
	// double hbias_arr[4];

 	// FILE *fp1;
 	// char filename1[FILENAME_MAX];

 	// for(i=0;i<nz;i++){
 	// 	sprintf(filename1,"%s/bias_halo_z%d.txt", gb.output_dir, (int)z[i]);
 	// 	fp1 = fopen(filename1, "ab");
 	// 	for(j=0;j<nM;j++){
	// 		halo_bias(&Cx_ref, MM[j], z[i], ST, hbias_arr);
	// 		b1  = hbias_arr[0];
	// 		b2  = hbias_arr[1];
	// 		bg2 = hbias_arr[2];
 	//          btd = hbias_arr[3];  
	// 	    	fprintf(fp1,"%12.6e %12.6e %12.6e %12.6e %12.6e \n",MM[j],b1,b2,bg2,btd);
	// 	}
	//  	fclose(fp1);
	// }      



 	///**--------------------------*///
      ////Compute the line biases as a function of redshift with and without the nfw profile	
      ///**--------------------------*///
	
	// int nz = 50;
 	// double *zz = init_1Darray(nz,0.0,11.);
 	// double b1_l[nlines][nz], b2_l[nlines][nz];

      // FILE *fp2;
      // char filename2[FILENAME_MAX];
      // sprintf(filename2,"%s/bias_line.txt", gb.output_dir);
      // fp2 = fopen(filename2, "ab");
	// double *lbias_arr = make_1Darray(2);

	// double *result = make_1Darray(3);
 	// double b1, b2, bg2, btd;

  	//fpr(j=0;j<nlines;j++){
	 	// for(i=0;i<nz;i++){
		  	// line_bias(Cx_ref.Lines[j], zz[i], lbias_arr);
		  	// b1_l[j][i]  = lbias_arr[0];
		  	// b2_l[j][i]  = lbias_arr[1];
		     	//printf("%12.6e %12.6e %12.6e \n",zz[i],b1_l[j][i],b2_l[j][i]);
		// }		
		// free(lbias_arr);
  	//}

	// for(i=0;i<nz;i++)
	//       fprintf(fp2," %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",zz[i],b1_l[0][i],b2_l[0][i],b1_l[1][i],b2_l[1][i],b1_l[2][i],b2_l[2][i],b1_l[3][i],b2_l[3][i],b1_l[4][i],b2_l[4][i],b1_l[5][i],b2_l[5][i],b1_l[6][i],b2_l[6][i]);

	// fclose(fp2);
	// free(zz);
	

	///**--------------------------*///
	///Compute the line biases, accounting for nfw profile 
	///for a fixed redshift as a function of k
	///**--------------------------*///
	// int nk          = 400;
 // 	int nz 	    = 3;
	// double *k       = loginit_1Darray(nk, 1.e-1, 40.);
	// double z[3]     = {2.,4.,6.}; 
 // 	double clust_pt = 0.;
 // 	double clust_hm = 0.;

 // 	double result[3];
	// double mom1, b1, b2, bg2, btd;
	// mom1 = mass_moment1(&Cx_ref, z[0], 1.e9, ST, CO10);

	// for(l=0;l<nk;l++){
	// 	HM_1h2h(&Cx_ref, k[l], z[0], 1.e9, ST, CO10, LINE, result);
	// 	b1 = result[1]/mom1;
	// 	b2 = result[2]/mom1;
	// 	bg2 =  -2./7.  * (b1-1.);
 // 	      btd =  23./42. * (b1-1.);  
	// 	printf("%12.6e %12.6e %12.6e %12.6e %12.6e \n", k[l], b1,b2,bg2,btd);
 // 	}	



 	///**--------------------------*///
	///Compute the line clustering power spectrum
	/// PS_line: PT, PS_HM: HM
	///**--------------------------*///

 	// FILE *fp1;
 	// char filename1[FILENAME_MAX];

 	// for(i=0;i<nlines;i++){
 	// 		for(j=0;j<nz;j++){
	//  		sprintf(filename1,"%s/line/ps/pk_pt_J%d_z%d.txt", gb.output_dir, JJ[i],(int)z[j]);
	//  		fp1    = fopen(filename1, "ab");

	// 		for(l=0;l<nk;l++){
	// 			clust_hm = PS_line_HM(&Cx_ref, k[l], z[j], M_min, mode_mf, lines[i]);
	// 			clust_pt = PS_line(&Cx_ref, k[l], z[j], i);
	// 			fprintf(fp1,"%12.6e %12.6e %12.6e \n",k[l], clust_pt, clust_hm);
	// 			// printf("%d %12.6e %12.6e %12.6e\n",j,k[i],clust_pt,clust_hm);

	// 		}
	// 	}	
	// }	
 	// fclose(fp1);
	// free(k);


	///**--------------------------*///
	///* Compute the 1loop matter power spectrum. The p13, p22, ptot are outputed inside pm_IR_NLO function *///
	///**--------------------------*///

      // double * ps_hloops = make_1Darray(6);
      // double pm_lin_IR,pm_1loop_IR;

      // int nk 	    = 100;
	// double *k    = loginit_1Darray(nk, 1.e-3, 20.);

      // for(i=0;i<nk;i++){
      // 	pm_lin_IR   = pm_IR_LO(&Cx_ref, k[i], z, GFILTER);
      // 	pm_1loop_IR = pm_IR_NLO(&Cx_ref, k[i], z, GFILTER);
      // 	printf("%12.6e %12.5e %12.6e %12.7e \n", k[i],Pk_dlnPk(&Cx_ref, k[i], z, LPOWER),pm_lin_IR,pm_1loop_IR);

      // 	Gloops_halo_interp(&Cx_ref, 0.1, z, GFILTER, 0, ps_hloops);

      // }
      // Gloops_halo_interp(&Cx_ref, 0.1, z, GFILTER, 0, ps_hloops);
      // free(ps_hloops);

	//free(k);


	///**--------------------------*///
	///* Test the halo model prediction of matter power spectrum, with and without 1loop*///
	///**--------------------------*///
	// long J 	    = 1;
	// int nk 	    = 100;
	// double z     = 1.;
	// double Mmin  = 1.e9;
	// long mode_mf = ST, mode_lum = CO10;
	// double *k    = loginit_1Darray(nk, 1.e-3, 20.);

	// double cs2  = 2.;
 	// double khat = 1.; //In unit of h/Mpc 
      // double Ms   = 1.e6;  //Lower mass limit for the halo model integration
      // double matter[2], hm_corrs[2];
      // double nfw, p1h_c, p2h_c;
      // double plin, pm_ct, pm_lin_IR, pm_1loop_IR;
      // double p1h, p2h, p2h_1loop, p2h_ct;

      // double Omegam = (Cx_ref.cosmo_pars[3] + Cx_ref.cosmo_pars[4]);
  	// double rhom_bar = Omegam * rhoc(0,&Cx_ref); 
      // mhmc(z,&Cx_ref,mode_mf,hm_corrs);

	// for(i=0;i<nk;i++){
	// 	nfw   = nfw_profile(k[i], Ms, z, &Cx_ref);
	// 	p1h_c = pow(nfw,2.)*Ms/rhom_bar * hm_corrs[0];
 	//    p2h_c = nfw * hm_corrs[1];
	// 	HM_1h2h(k[i], z, Mmin, &Cx_ref, mode_mf, mode_lum, MATTER, matter);  /// in unit of M_sun/Mpc^3
		
	// 	plin        = Pk_dlnPk(k[i], z, &Cx_ref, LPOWER) * pow(gb.h,3.);
	// 	pm_lin_IR   = pm_IR_LO(k[i], z, k[0], &Cx_ref, GFILTER) * pow(gb.h,3.);
      //    pm_1loop_IR = pm_IR_NLO(k[i], z, k[0], &Cx_ref, GFILTER) * pow(gb.h,3.) ;
      //    pm_ct       = - 2. * cs2 * pow(k[i]/gb.h, 2.)/(1.+pow((k[i]/gb.h)/khat,2.)) * pm_lin_IR;  //Eq. 2.22 of 2004.09515, Pade approximation to make p_ct finite at large k

	// 	p1h         = (pow(rhom_bar,-2.) * matter[0] + p1h_c) * pow(gb.h,3.);
	// 	p2h         = pow(matter[1]/rhom_bar + p2h_c, 2.) * plin;
	// 	p2h_1loop   = pow(matter[1]/rhom_bar + p2h_c, 2.) * pm_1loop_IR ;
	// 	p2h_ct      = pow(matter[1]/rhom_bar + p2h_c, 2.) * pm_ct;

	// 	printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n", k[i]/gb.h, p1h, p2h, p2h_1loop, p2h_ct, pm_1loop_IR, pm_ct, plin);
	// }	
	//free(k);
					




	///*************************************************************************************/////
						//// END OF DEBUGGING ////
	///*************************************************************************************/////

		
	
