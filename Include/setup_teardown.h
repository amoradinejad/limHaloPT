
#ifndef _SETUP_TEARDOWN_H_
#define _SETUP_TEARDOWN_H_

#define CO10            131L
#define CO21            132L
#define CO32            133L
#define CO43            134L
#define CO54            135L
#define CO65            136L
#define CII             137L

typedef struct initialize_struct{
      double Mh_min;
      long   mode_mf;
      size_t ninterp;
      int    nlines;
      int    lines[7];
}initialize_struct;


/** 
 * Function declarations of main.c module
 */
void  initialize(char *argv[], struct initialize_struct *init_struct);
void  cleanup();


#endif