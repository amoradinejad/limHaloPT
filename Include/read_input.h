
#ifndef _READ_INPUT_H_
#define _READ_INPUT_H_


typedef struct input_param{
    char* name;
    char* value;
    double nvalue;
} InputParam;


unsigned int test_float(const char *str);
int check_for_one_equal_sign(char* line_in);
char* trim(char* in);
void parse_line(char* line_in, char** pname, char** pvalue);
unsigned int lines_with_parameters(FILE* fp);

#endif
