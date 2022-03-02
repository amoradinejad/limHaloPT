# NOTE: Depending on your installation paths for gsl and cuba libraries, you should add the links 
# for the INCLUDES and LIBS below, replacing -I/opt/local/include and  -L/opt/local/lib which are the 
# paths on my machine. 
INCLUDES =  -I../Include/ -I/opt/local/include  -I../Class/include  
LIBS =  -L/opt/local/lib  -L../Class/lib
CFLAGS = -g -fopenmp -Wall
FLAGS =  -lgsl -lgslcblas -lcuba  -lclass -lm -o
OFILESA = main.o cubature.o utilities.o cosmology.o survey_specs.o 
OFILESB = line_ingredients.o wnw_split.o ir_res.o ps_halo_1loop.o ps_line_pt.o ps_line_hm.o 

limHaloPT: ${OFILESA} ${OFILESB}
	gcc ${OFILESA}${OFILESB} $(INCLUDES) $(LIBS) $(CFLAGS) $(FLAGS) limHaloPT       
main.o: main.c 
	gcc -c main.c -o main.o $(INCLUDES) 
cubature.o: hcubature.c ../Include/cubature.h 
	gcc -c hcubature.c -o cubature.o $(INCLUDES)
utilities.o: utilities.c ../Include/utilities.h
	gcc -c utilities.c -o utilities.o $(INCLUDES)
cosmology.o: cosmology.c ../Include/cosmology.h
	gcc -c cosmology.c -o cosmology.o $(INCLUDES)
survey_specs.o: survey_specs.c ../Include/survey_specs.h
	gcc -c survey_specs.c -o survey_specs.o $(INCLUDES)	
line_ingredients.o: line_ingredients.c ../Include/line_ingredients.h
	gcc -c line_ingredients.c -o line_ingredients.o $(INCLUDES) 
wnw_split.o: wnw_split.c ../Include/wnw_split.h
	gcc -c wnw_split.c -o wnw_split.o $(INCLUDES)
ir_res.o: ir_res.c ../Include/ir_res.h
	gcc -c ir_res.c -o ir_res.o $(INCLUDES)
ps_halo_1loop.o: ps_halo_1loop.c ../Include/ps_halo_1loop.h
	gcc -c ps_halo_1loop.c -o ps_halo_1loop.o $(INCLUDES)
ps_line_pt.o: ps_line_pt.c ../Include/ps_line_pt.h
	gcc -c ps_line_pt.c -o ps_line_pt.o $(INCLUDES)
ps_line_hm.o: ps_line_hm.c ../Include/ps_line_hm.h
	gcc -c ps_line_hm.c -o ps_line_hm.o $(INCLUDES)

clean:
	rm *.o limHaloPT
