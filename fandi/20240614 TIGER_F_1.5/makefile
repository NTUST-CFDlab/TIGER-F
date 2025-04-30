# Makefile 
# Nom du compilateur

# ----------------- Intel compiler -------------------
FC = mpiifort

# Options de compilation: optimisation, debug etc...
OPT = -O2 -msse4.2 -axAVX,CORE-AVX2,CORE-AVX512 -heap-arrays -fopenmp -mcmodel=large -traceback


# -ipo -xCORE-AVX2 -xCORE-AVX512 -align array32byte -qopenmp -Ofast
# -mcmodel=large
# -heap-arrays 64
# ----------------------------------------------------


#   *************************************************
#   |  execute with: nohup mpirun -np <n> ./run &	|
#   *************************************************
# ----------------------------------------------------

EXE = run

LINKOPT = 

DIR_OBJ = ./obj
DIR_BIN = ./bin

# Defining the objects (OBJS) variables
OBJS =  \
    $(DIR_OBJ)/variables_module.o \
    $(DIR_OBJ)/main.o \
    $(DIR_OBJ)/MPI_division.o \
    $(DIR_OBJ)/read_data.o \
    $(DIR_OBJ)/gridder.o \
    $(DIR_OBJ)/initial_conditions.o \
    $(DIR_OBJ)/boundary_conditions.o \
    $(DIR_OBJ)/vos_function.o \
    $(DIR_OBJ)/vos_ray2d.o \
    $(DIR_OBJ)/vos_ray3d.o \
    $(DIR_OBJ)/Smagorinsky.o \
    $(DIR_OBJ)/QUICK.o \
    $(DIR_OBJ)/AdamsBashforth.o \
    $(DIR_OBJ)/prediction_correction.o \
    $(DIR_OBJ)/calcul_new_velocity.o \
    $(DIR_OBJ)/virtualForceIntegrator.o \
    $(DIR_OBJ)/filer.o \
    $(DIR_OBJ)/SOR.o \
    #$(DIR_OBJ)/check_steady.o \
    #$(DIR_OBJ)/BICGstab.o \
    #$(DIR_OBJ)/solver_1st_member.o \
    #$(DIR_OBJ)/solver_2nd_member.o \
    #$(DIR_OBJ)/solver.o \
    #$(DIR_OBJ)/solver_pressure_field.o\


# Linking object files
$(DIR_BIN)/$(EXE) :   $(OBJS)
	$(FC) $(OPT)   -o $(EXE) \
	$(OBJS) \
	$(OPT)
    
	@echo " "
	@echo "    ************** TIGER-F compiled successfully! ************** "                                                                                                                                      
	#   |  Version: 1.5                                            |
	#   |  Author : CFD Lab., Taiwan Tech.                         |
	#   |  Web    : http://smetana.me.ntust.edu.tw/                |
	#   |  Contributors:                                           |
	#   |     Prof. Ming-Jyh Chern, Zi-Hsuan Wei, Jing-Ming Chen,  |
	#   |     Desta Goytom Tewolde, Fandi Dwiputra Suprianto       |
    #    


# Rule to compile .f90 to .o
$(DIR_OBJ)/%.o:%.f90
	@mkdir -p $(DIR_OBJ)
	$(FC) $(OPT) -c -o $@ $<

main.o : variables_module.o

MPI_division.o : variables_module.o

read_data.o : variables_module.o

gridder.o : variables_module.o

initial_conditions.o : variables_module.o

boundary_conditions.o : variables_module.o

vos_function.o : variables_module.o

vos_ray2d.o : variables_module.o

vos_ray3d.o : variables_module.o

Smagorinsky.o : variables_module.o

QUICK.o : variables_module.o

AdamsBashforth.o : variables_module.o

calcul_new_velocity.o : variables_module.o

prediction_correction.o : variables_module.o

virtualForceIntegrator.o : variables_module.o

filer.o : variables_module.o

SOR.o : variables_module.o

#check_steady.o : variables_module.o

#BICGstab.o : variables_module.o

# Removing object files
clean :
	/bin/rm -rf $(DIR_OBJ)
	/bin/rm -f $(EXE)  *.mod

cleanall : 
	/bin/rm -rf $(DIR_OBJ)
	/bin/rm -f $(EXE)  *.mod
	/bin/rm -f *.dat
	/bin/rm -f *.out
	/bin/rm -rf output*/

    
config :
	if [ ! -d obj ] ; then mkdir obj ; fi


