
/** @file setup_teardown.c Documented setup_teardown module, including the initialize and clean functions that are called first and last thing in limHaloPT package
 *
 * Azadeh Moradinezhad Dizgah & Alberto Vallinotto, April 4th 2022
 * 
 */

#include "header.h"
#include "read_input.h"


/**
 * Initialize the path to the required directories,and the interpolator of the SFR(M,z) from tabulated data provided in gb.SFR_filename.  
 * Read in the .ini file and set the values of parameters.
 * 
 * The global structure "gb" have several elements to hold the paths to project source directory, input, and output folders, 
 * and values of cosmological parmaeters. 
 * 
 * @return void 
 */
void initialize(char *argv[], struct initialize_struct *init_struct)
{
      getcwd(gb.project_home , sizeof(gb.project_home));  //This gives the path to the source directory      
      sprintf(gb.output_dir,"%s/Output", gb.project_home);
      sprintf(gb.data_dir,"%s/Input", gb.project_home);
      sprintf(gb.SFR_filename,"%s/release-sfh_z0_z8_052913/sfr/sfr_release.dat",gb.data_dir);  

      logSFR_alloc_init();


      FILE* fp;
      unsigned int n_parameters;
      unsigned int p_id = 0;
      InputParam* p_array;

      printf("argv[0] = '%s'\n", argv[0]);
      printf("argv[1] = '%s'\n", argv[1]);

      fp = fopen(argv[1], "r");

      n_parameters = lines_with_parameters(fp);

      printf("n_parameters = %d\n", n_parameters);

      p_array     = (InputParam*)malloc(sizeof(InputParam) * n_parameters);
    


      while(!feof(fp)){
            char line[1024];
            char* pname = NULL;
            char* pvalue = NULL;

            fgets(line, 1024, fp);
            parse_line(line, &pname, &pvalue);

            if(pname != NULL && pvalue != NULL){
                  p_array[p_id].name = pname;
                  p_array[p_id].value = pvalue;
                  p_array[p_id].nvalue = (test_float(p_array[p_id].value) ? atof(p_array[p_id].value) : NAN);
                  p_id++;
            } // End of if

      } // End of while

      /**
       * Commment out to see the names, values and numerical values of the parameters read from .ini file.
       */
      // for(unsigned c=0; c<n_parameters; c++){
      //       printf("\n------------------------------\n");
      //       printf("p_array[%d].name   = '%s'\n", c, p_array[c].name);
      //       printf("p_array[%d].value  = '%s'\n", c, p_array[c].value);
      //       printf("p_array[%d].nvalue = %12.6e\n", c, p_array[c].nvalue);
      // }

      // printf("\n");


      char *lines_input_switches[] = {"line1", "line2", "line3", "line4", "line5", "line6", "line7"};
      char *lines_input_names[]    = {"CO10", "CO21", "CO32", "CO43", "CO54", "CO65", "CII"};
      int line_vals[7]             = {CO10, CO21, CO32, CO43, CO54, CO65, CII};
      static int line_read = 0;
      /**
       * Save the values of input parmaeters to global structure or init_struct. The latter is returned as the output of the initialize() function.
       */ 
      for(unsigned c=0; c<n_parameters; c++)
      {
            if (strcmp(p_array[c].name, "h") == 0){
                  gb.h = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.h);
            }
            if (strcmp(p_array[c].name, "Ocdm") == 0){
                  gb.Omega_cdm = p_array[c].nvalue/pow(gb.h,2.);
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.Omega_cdm);
            }
            if (strcmp(p_array[c].name, "Ob") == 0){
                  gb.Omega_b = p_array[c].nvalue/pow(gb.h,2.);
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.Omega_b);
            }
            if (strcmp(p_array[c].name, "As") == 0){
                  gb.As = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.As);
            }
            if (strcmp(p_array[c].name, "ns") == 0){
                  gb.ns = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.ns);
            }
            if (strcmp(p_array[c].name, "Npars") == 0){
                  gb.Npars = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.Npars);
            }
            if (strcmp(p_array[c].name, "kp") == 0){
                  gb.kp = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.kp);
            }
            if (strcmp(p_array[c].name, "sigFOG0") == 0){
                  gb.sigFOG0 = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.sigFOG0);
            }
            if (strcmp(p_array[c].name, "PS_zmax") == 0){
                  gb.PS_zmax = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.PS_zmax);
            }
            if (strcmp(p_array[c].name, "line_zmax") == 0){
                  gb.line_zmax = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.line_zmax);
            }
            if (strcmp(p_array[c].name, "mean_zmax") == 0){
                  gb.mean_zmax = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.mean_zmax);
            }
            if (strcmp(p_array[c].name, "mean_nz") == 0){
                  gb.mean_nz = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.mean_nz);
            }
            if (strcmp(p_array[c].name, "power_zmax") == 0){
                  gb.power_zmax = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.power_zmax);
            }
            if (strcmp(p_array[c].name, "power_nz") == 0){
                  gb.power_nz = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.power_nz);
            }
            if (strcmp(p_array[c].name, "power_kmax") == 0){
                  gb.power_kmax = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.power_kmax);
            }
            if (strcmp(p_array[c].name, "power_kmin") == 0){
                  gb.power_kmin = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.power_kmin);
            }
            if (strcmp(p_array[c].name, "power_nk") == 0){
                  gb.power_nk = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, gb.power_nk);
            }
            if (strcmp(p_array[c].name, "mode_mf") == 0){ 
                  if(strcmp(p_array[c].value, "PSC") == 0)
                       init_struct->mode_mf = PSC;
                  else if(strcmp(p_array[c].value, "ST") == 0)
                       init_struct->mode_mf = ST;
                  else if(strcmp(p_array[c].value, "TK") == 0)
                       init_struct->mode_mf = TK;
                  printf("p_array[%d].name   = '%s' %d\n",c, p_array[c].name, init_struct->mode_mf);
            }
            if (strcmp(p_array[c].name, "Mh_min") == 0){
                  init_struct->Mh_min = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %12.6e\n",c, p_array[c].name, init_struct->Mh_min);
            }
            if (strcmp(p_array[c].name, "ninterp") == 0){
                  init_struct->ninterp = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %d\n",c, p_array[c].name,init_struct->ninterp);
            }
            if (strcmp(p_array[c].name, "nlines") == 0){
                  init_struct->nlines = p_array[c].nvalue;
                  printf("p_array[%d].name   = '%s' %d \n",c, p_array[c].name,init_struct->nlines);
                  line_read = 1;
            }
            if(line_read == 1){
                  int l = 0;
                  for(int i=0;i<7;i++){
                        if (strcmp(p_array[c].name, lines_input_switches[i]) == 0){
                              for(int j=0;j<7;j++){
                                    if(strcmp(p_array[c].value, lines_input_names[j]) == 0){    
                                          init_struct->lines[l] = line_vals[j];
                                          printf("p_array[%d].name   = '%s' %d\n",c, p_array[c].name, init_struct->lines[i]);
                                    }      
                              }
                        }
                        l += 1;            
                  }
            }
            if (strcmp(p_array[c].name, "switch_Tbar") == 0){
                  if(strcmp(p_array[c].value, "Tbar") == 0)
                        gb.switch_Tbar = MEAN;
                  else 
                        gb.switch_Tbar = 0.;  
                  printf("p_array[%d].name   = '%s' %d\n",c, p_array[c].name, gb.switch_Tbar);
            }
            if (strcmp(p_array[c].name, "switch_HMshot") == 0){
                  if(strcmp(p_array[c].value, "HMshot") == 0)
                        gb.switch_HMshot = HMshot;
                  else 
                        gb.switch_HMshot = 0.;  
                  printf("p_array[%d].name   = '%s' %d\n",c, p_array[c].name, gb.switch_HMshot);
            }
            if (strcmp(p_array[c].name, "switch_HMclust") == 0){
                  if(strcmp(p_array[c].value, "HMclust") == 0)
                        gb.switch_HMclust = HMclust;
                  else 
                        gb.switch_HMclust = 0.;  
                  printf("p_array[%d].name   = '%s' %d\n",c, p_array[c].name, gb.switch_HMclust);
            }

            // printf("p_array[%d].name   = '%s' %12.6e %s \n \n",c, p_array[c].name, p_array[c].value, p_array[c].nvalue);

      }      

      /**
       * Cosmological parameters corresponding to initial conditions of HiddenValley (HV) simulations
       */
      gb.c         = 2.99792458e5;             // In units of km/s
      gb.Omega_r   = 0.0000910958;             // radiation = photons + neutrinos
      gb.Omega_g   = 5.3956421715871286e-05;   // photons, input for Class
      gb.Omega_nu  = 0.00;    ////neutrinos
      gb.logAs     = log(gb.As*pow(10.,10.));  ///3.0665

      gb.PS_kmin = 1.e-4;
      gb.PS_kmax = 200.;
      gb.kmax_CLASS = PS_KMAX;
      

      return;
}


/**
 * Free the memory allocated to cosmology structure and SFR interpolator
 * 
 * @return void 
 */
void cleanup(struct Cosmology *Cx)
{
      Cosmology_free(Cx);
      SFR_Behroozi_free();

      return;
}



