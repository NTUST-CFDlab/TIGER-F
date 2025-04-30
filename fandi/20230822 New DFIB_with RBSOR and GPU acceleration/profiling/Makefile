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


# ----------------- NVIDIA compiler ------------------
#FC = mpifort

# Options de compilation: optimisation, debug etc...
#OPT = -O3 -mp=multicore -acc=gpu -gpu=pinned -fPIC -mcmodel=medium -Minfo=accel -traceback


# -fast -O2 -Munroll -Mvect (O2, Munroll & Mvect are included in fast)
# -gpu=ccxx (xx = GPU compute capability) -gpu=ccxx,managed -gpu=ccxx,pinned

#   *****************************************************************
#   |  execute with: nohup mpirun --bind-to none -np <n> ./run &	|
#   *****************************************************************
# ----------------------------------------------------

EXE = run

LINKOPT = 

DIR_OBJ = ./obj
DIR_BIN = ./bin

# Defining the objects (OBJS) variables
OBJS =  \
    variables_module.o \
    main.o \
    Mpi_division.o \
    reading_data.o \
    gridder.o \
    initial_conditions.o \
    boundary_conditions.o \
    vos_function.o \
    vos_ray2d.o \
    vos_ray3d.o \
    smagorinsky.o \
    discretisation_QUICK_centre.o \
    AdamsBashforth.o \
    calcul_new_velocity.o \
    virtualForceIntegrator.o \
    filer.o \
    GS.o \
    #prediction_force.o \
    #check_steady.o \
    #plasma.o \
    #BICGstab.o \
    #solver_1st_member.o \
    #solver_2nd_member.o \
    #solver.o \
    #solver_pressure_field.o 


# Linking object files
exe :   $(OBJS)
	$(FC) $(OPT)   -o $(EXE) \
    $(OBJS) \
    $(OPT)



	@echo "   ***************** successful *****************   "                                                                                      
	#    |  Author:  Zi-Hsuan Wei                                                 
	#    |  Version: 1.8                                                         
	#    |  Web:     http://smetana.me.ntust.edu.tw/     
	#    |  Editor : Hsuan
    #    


%.o:%.f90
	$(FC) $(OPT) -c $<

main.o : variables_module.o

Mpi_division.o : variables_module.o

reading_data.o : variables_module.o

gridder.o : variables_module.o

initial_conditions.o : variables_module.o

boundary_conditions.o : variables_module.o

vos_function.o : variables_module.o

vos_ray2d.o : variables_module.o

vos_ray3d.o : variables_module.o

smagorinsky.o : variables_module.o

discretisation_QUICK_centre.o : variables_module.o

AdamsBashforth.o : variables_module.o

calcul_new_velocity.o : variables_module.o

virtualForceIntegrator.o : variables_module.o

filer.o : variables_module.o

GS.o : variables_module.o

#prediction_force.o : variables_module.o

#check_steady.o : variables_module.o

#BICGstab.o : variables_module.o

# Removing object files
clean :
	/bin/rm -f *.dat
	/bin/rm -f *.x
	/bin/rm -f *.q

cleanall : 
	/bin/rm -f $(OBJS) $(EXE)  *.mod
	/bin/rm -f *.dat
	/bin/rm -f *.x
	/bin/rm -f *.q
	/bin/rm -f *.out
    
config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir bin ; fi


