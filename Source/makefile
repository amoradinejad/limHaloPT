INCLUDES =  -I../Include/ -I/opt/local/include  -I../Class/include  
LIBS =  -L/opt/local/lib  -L../Class/lib
FLAGS =  -lgsl -lgslcblas -lcuba  -lclass -lm -o
OFILESA = main.o cubature.o utilities.o cosmology.o survey_specs.o 
OFILESB = line_ingredients.o wnw_split.o ir_res.o ps_halo_1loop.o ps_line_pt.o ps_line_hm.o 

limHaloPT: ${OFILESA} ${OFILESB}
	mpicc ${OFILESA}${OFILESB} -fopenmp -Wall $(INCLUDES) $(LIBS)  $(FLAGS) limHaloPT       
main.o: main.c 
	mpicc -g -c -fopenmp main.c -o main.o $(INCLUDES) 
cubature.o: hcubature.c ../Include/cubature.h 
	mpicc -g -c -fopenmp hcubature.c -o cubature.o $(INCLUDES)
utilities.o: utilities.c ../Include/utilities.h
	mpicc -g -c -fopenmp utilities.c -o utilities.o $(INCLUDES)
cosmology.o: cosmology.c ../Include/cosmology.h
	mpicc -g -c -fopenmp cosmology.c -o cosmology.o $(INCLUDES)
survey_specs.o: survey_specs.c ../Include/survey_specs.h
	mpicc -g -c -fopenmp survey_specs.c -o survey_specs.o $(INCLUDES)	
line_ingredients.o: line_ingredients.o ../Include/line_ingredients.h
	mpicc -g -c  -fopenmp line_ingredients.c -o line_ingredients.o $(INCLUDES) 
wnw_split.o: wnw_split.c ../Include/wnw_split.h
	mpicc -g -c  -fopenmp wnw_split.c -o wnw_split.o $(INCLUDES)
ir_res.o: ir_res.c ../Include/ir_res.h
	mpicc -g -c  -fopenmp ir_res.c -o ir_res.o $(INCLUDES)
ps_halo_1loop.o: ps_halo_1loop.c ../Include/ps_halo_1loop.h
	mpicc -g -c  -fopenmp ps_halo_1loop.c -o ps_halo_1loop.o $(INCLUDES)
ps_line_pt.o: ps_line_pt.c ../Include/ps_line_pt.h
	mpicc -g -c  -fopenmp ps_line_pt.c -o ps_line_pt.o $(INCLUDES)
ps_line_hm.o: ps_line_hm.c ../Include/ps_line_hm.h
	mpicc -g -c  -fopenmp ps_line_hm.c -o ps_line_hm.o $(INCLUDES)

clean:
	rm *.o limHaloPT
