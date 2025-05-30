# Makefile 
# Nom du compilateur


# ----------------- NVIDIA compiler ------------------
FC = mpifort

# Options de compilation: optimisation, debug etc...
# Run on GPU (OpenACC):
	OPT = -O2 -acc=gpu,multicore -gpu=pinned -cuda -mcmodel=medium -Minfo=accel -traceback
	#OPT = -O2 -acc=gpu,multicore -gpu=mem:separate:pinnedalloc -cuda -mcmodel=medium -Minfo=accel -traceback

# Run on CPU (OpenACC):
	#OPT = -O2 -acc=multicore -cuda -mcmodel=medium -Minfo=accel -traceback

# Run on CPU (OpenMP):
	#OPT = -O2 -mp=multicore -cuda -mcmodel=medium -Minfo=accel -traceback

# -fast -O3 -Munroll -Mvect (O2, Munroll & Mvect are included in fast)
# -gpu=ccxx (xx = GPU compute capability) -gpu=ccxx,managed -gpu=ccxx,pinned

#   *********************************************************************************************
#   |  execute with: 																			|
#	|		  np=1 : nohup mpirun --bind-to none -np 1 env OMP_PROC_BIND=spread,close ./run & 	|
#	|		  np>1 : nohup mpirun --bind-to none -np <np> ./run & 								|
#   *********************************************************************************************
#
#	profiling with Nsight:
#	Compile with:
#		OPT = -O2 -acc=gpu,multicore -gpu=pinned -cuda -lnvhpcwrapnvtx -mcmodel=medium -Minfo=accel -traceback
#	Nsight systems, run with:
#		nohup mpirun --bind-to none -np 1 env OMP_PROC_BIND=spread,close nsys profile -t nvtx,cuda,openacc --stats=true --force-overwrite true -o my_report ./run > nsys_log.txt 2>&1 &
#	Nsight compute, run with:
#		nohup ncu --clock-control none --set full -o profile mpirun --bind-to none -np 1 ./run > ncu_log.txt 2>&1 &
# ----------------------------------------------------



# ----------------- Intel compiler -------------------
#FC = mpiifort

# Options de compilation: optimisation, debug etc...
	#OPT = -O2 -msse4.2 -axAVX,CORE-AVX2,CORE-AVX512 -heap-arrays -fopenmp -mcmodel=large -traceback


# -ipo -xCORE-AVX2 -xCORE-AVX512 -align array32byte -qopenmp -Ofast
# -mcmodel=large
# -heap-arrays 64
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
    $(DIR_OBJ)/dynamic_coupling.o \
    $(DIR_OBJ)/vos_function.o \
    $(DIR_OBJ)/vos_ray2d.o \
    $(DIR_OBJ)/vos_ray3d.o \
    $(DIR_OBJ)/Smagorinsky.o \
    $(DIR_OBJ)/QUICK.o \
    $(DIR_OBJ)/AdamsBashforth.o \
    $(DIR_OBJ)/SOR.o \
    $(DIR_OBJ)/calcul_new_velocity.o \
    $(DIR_OBJ)/prediction_correction.o \
    $(DIR_OBJ)/virtualForceIntegrator.o \
    $(DIR_OBJ)/wall_model.o \
    $(DIR_OBJ)/filer.o \
    #$(DIR_OBJ)/SOR_fp32.o \
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
	#   |  Version: 3.2                                            |
	#   |  Author : CFD Lab., Taiwan Tech.                         |
	#   |  Web    : https://mestaiwantech.com                      |
	#   |  Contributors:                                           |
	#   |     Prof. Ming-Jyh Chern, Zi-Hsuan Wei, Jing-Ming Chen,  |
	#   |     Desta Goytom Tewolde, Fandi Dwiputra Suprianto,      |
	#   |     Adhika Satyadharma, Ming-Fang Kuan, Tai-Yi Chou	   |
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

dynamic_coupling.o : variables_module.o

vos_function.o : variables_module.o

vos_ray2d.o : variables_module.o

vos_ray3d.o : variables_module.o

Smagorinsky.o : variables_module.o

QUICK.o : variables_module.o

AdamsBashforth.o : variables_module.o

SOR.o : variables_module.o

calcul_new_velocity.o : variables_module.o

prediction_correction.o : variables_module.o

virtualForceIntegrator.o : variables_module.o

wall_model.o : variables_module.o

filer.o : variables_module.o


#SOR_fp32.o : variables_module.o

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


