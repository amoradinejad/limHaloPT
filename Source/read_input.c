

/** @file setup_teardown.c Documented read_input module. This module provides utilityy functions to read the .ini file 
 *
 * Alberto Vallinotto, April 4th 2022
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "read_input.h"


/**
 * See
 * https://stackoverflow.com/questions/45554639/is-there-a-way-to-check-if-a-string-can-be-a-float-in-c
 */

unsigned int test_float(const char *str)
{
    int len;
    float dummy = 0.0;
    if (sscanf(str, "%f %n", &dummy, &len) == 1 && len == (int)strlen(str))
        return 1;
    else
        return 0;
}



int check_for_one_equal_sign(char* line_in)
{
    int retval = 0;
    for(char* c = line_in; *c != '\0'; c++){
        if(*c == '='){
            retval++;
        }
    }

    return (retval > 1 ? 0 : retval);
}



char* trim(char* in)
{ 
    
    unsigned i=0;
    unsigned l=0;
    char* c = in;

    //printf("in = '%s'\n", in);

    // Count how many characters to copy
    while(*c != '\0'){
        if(!(*c == ' ' || *c == '\n')){
            i++;
        }
        c++;        
    }

    l = i;
    //printf("l = %d\n", l);

    // Allocate string
    char* out = (char*)malloc(sizeof(char) * l);

    // Copy character by character skipping spaces
    c = in;
    i=0;
    while(i<l){
        if(!(*c == ' ' || *c == '\n')){
            out[i] = *c;
            i++;
        }
        c++;        
    }

    return out;
}


void parse_line(char* line_in, char** pname, char** pvalue)
{
    char* before;
    char* after;

    if(line_in[0] != '#' && check_for_one_equal_sign(line_in) == 1){

        before = strtok_r(line_in, "=", &after);
        //printf("before = '%s'\n", before);
        //printf("after = %s\n", after);

        *pname = trim(before);
        *pvalue = trim(after);

        //printf("pname = '%s'\n",  *pname);
        //printf("pvalue = '%s'\n", *pvalue);
    }

    return;
}



unsigned int lines_with_parameters(FILE* fp)
{
    unsigned n_lines = 0;

    while(!feof(fp)){
        char line[1024];
        fgets(line, 1024, fp);

        if(check_for_one_equal_sign(line)){
            n_lines++;
        }
    }

    rewind(fp);
    return n_lines;
}



