#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ options to choose ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# -- GPU --
#arch = gpu
# extreme scale simulation results in very large array
#x_scale = enabled    

# -- CPU --
arch = cpu
# whether or not use avx512 for intel cpu, test before use 
openmp = enabled
#AVX512 = enabled

# -- MIC(KNL) --
# arch = mic
# # OpenMP and AVX512 should both be enabled for MIC
# openmp = enabled
# AVX512 = enabled

# -- compiler --
#tool = gcc
#tool = pgi
tool = intel
#tool = cray  (not implemented)

# -- optimization options --
#optimization = debug
optimization = release
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ options to choose ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

############# DIR #################################
# Directories
DIR=./
OBJ_DIR = $(DIR)1.obj/
EXEC_DIR = $(DIR)1.exec/
SRC_DIR = $(DIR)0.src/
############# DIR #################################

############# architecture ########################
ifeq ($(strip $(arch)), cpu)
	BIN_SUFFIX = cpu
	ifeq ($(strip $(AVX512)), enabled)
		ARCH_FLAGS = -xCORE-AVX512		
	endif
endif

ifeq ($(strip $(arch)), gpu)
	BIN_SUFFIX = gpu
	openmp = disabled
#compatible with pgi compiler, for cray compiler: ARCH_FLAGS = -h acc,noomp -h msgs
	ARCH_FLAGS = -acc -Mcuda -Minfo=accel -ta=tesla:pinned   
endif

ifeq ($(strip $(arch)), mic)
	BIN_SUFFIX = mic
	ifeq ($(strip $(AVX512)), enabled)
		ARCH_FLAGS = -xMIC-AVX512		
	endif
endif

############# architecture ########################

############# tool ################################
# gcc
ifeq ($(strip $(tool)), gcc)
# -J 0.obj put .mod files to 0.obj directory
	MISC_FLAGS = -ffree-line-length-300 -J $(OBJ_DIR)  
#for gcc10, add the following temporary fix: MISC_FLAGS += -fallow-argument-mismatch   

	ifeq ($(strip $(optimization)), debug)
		OPT_FLAGS = -g
	else	
		OPT_FLAGS = -O2
    endif
	
	ifeq ($(strip $(openmp)), enabled)
		OMP_FLAGS = -fopenmp
	endif
endif

# pgi
ifeq ($(strip $(tool)), pgi)
	MISC_FLAGS = -module $(OBJ_DIR) 

	ifeq ($(strip $(optimization)), debug)
		OPT_FLAGS = -g
	else 
		OPT_FLAGS = -fast
    endif

	ifeq ($(strip $(openmp)), enabled)
		OMP_FLAGS = -mp
	endif

#for very large array, add: MISC_FLAGS = -mcmodel=medium 
	ifeq ($(strip $(x_scale)), enabled)
		MISC_FLAGS += -mcmodel=medium 
	endif

endif

# intel
ifeq ($(strip $(tool)), intel)
	MISC_FLAGS = -module $(OBJ_DIR) 

	ifeq ($(strip $(optimization)), debug)
		OPT_FLAGS = -debug
	else 
		OPT_FLAGS = -O2
    endif

	ifeq ($(strip $(openmp)), enabled)
		OMP_FLAGS = -qopenmp
	endif

endif
############# tool ################################

compiler = mpif90
#compiler = ftn
#compiler = mpiifort

CFLAGS += $(ARCH_FLAGS)
CFLAGS += $(OPT_FLAGS)
CFLAGS += $(OMP_FLAGS)
CFLAGS += $(MISC_FLAGS)
CFLAGS += $(FPP_FLAGS)
LFLAGS = $(CFLAGS)

COMPILE = ${compiler} 
LINK = ${compiler}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


objects = Module.o Main_multiphase.o Init_multiphase.o Misc.o IO_multiphase.o Kernel_multiphase.o Monitor.o Mpi_misc.o Mpi.o \
Boundary_multiphase_inlet.o Boundary_multiphase_outlet.o Boundary_multiphase_other.o Phase_gradient.o Geometry_preprocessing.o

pref_objects = $(addprefix $(OBJ_DIR), $(objects))

LBM : | make_dir $(pref_objects)
	${LINK} $(LFLAGS) -o $(EXEC_DIR)MF_LBM.$(BIN_SUFFIX) $(pref_objects)


$(OBJ_DIR)%.o: $(SRC_DIR)%.F90
	${COMPILE} $(CFLAGS) -c $< -o $@
# $@ is the name of the target being generated
# $< the first prerequisite (usually a source file)


.PHONY: clean 
clean:
	rm -r $(OBJ_DIR)
	rm -r $(EXEC_DIR)

make_dir:
	@echo "Create directories."
	mkdir -p $(DIR)1.exec
	mkdir -p $(DIR)1.obj

