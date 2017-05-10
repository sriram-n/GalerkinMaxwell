#####-----------------------------------------

COMPLEX         = 1
PARALLEL        = 0

-include ../../m_options

SRC_PATH        = .
OBJ_PATH        = _obj_
MODULE_PATH     = $(OBJ_PATH)

# executable name
EXEC            = galMaxwell
EXECS           = galMaxwell

# Includes and defs
CCDEFS          = -D C_MODE=$(COMPLEX) \
                  -D PARALLEL_MODE=$(PARALLEL)

CCINCS          = -I$(HP3D_PATH)/common/hp3d \
                  -I$(HP3D_PATH)/common/graphics \
                  -I$(HP3D_PATH)/module \
                  -I$(MUMPS_MPI_INC) \
                  -I./$(MODULE_PATH)

###### User defined libraries
LIBS            = $(PROB_LIBS)
#LIBS            = $(HP3D_LIB) $(HP3D_COMMON) $(HP3D_GMP) $(HP3D_COMMON) \
#		  $(HP3D_LIB) $(VIS_LIB) \
#		  $(UHM_LIB) $(LINAL_LIB) $(FLA_LIB) \
#		  $(MUMPS_LIB) $(MUMPS_MPI_LIB) \
#		  $(PARDISO_LIB) $(WSMP_LIB) \
#		  $(METIS_LIB) \
#		  $(MKL_LIB) $(LAPACK_LIB) $(BLAS_LIB) \
#		  $(X_LIB) $(PTHREAD_LIB)

ifdef CUDA_PATH
LIBS           += $(CUDA_LIB)
endif


include make.inc

SRC         = $(addprefix $(SRC_PATH)/,$(SRC_FILE))

## Kyungjoo added.
SRC_FILE_F  = $(filter %.F, $(SRC_FILE))
SRC_FILE_F90= $(filter %.F90, $(SRC_FILE))

OBJ_F90     = $(addprefix $(OBJ_PATH)/, $(SRC_FILE_F90:.F90=.F90.o))
OBJ_F       = $(addprefix $(OBJ_PATH)/, $(SRC_FILE_F:.F=.F.o))

OBJ         = $(OBJ_F) $(OBJ_F90)

## OBJ         = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F90=.o))


FC_WORK     = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

one : $(OBJ_PATH)/.dummy $(OBJ)
	@echo " "
	@echo "- Linking - ", $(EXEC)
	@echo "- COMPLEX, PARALLEL ", $(COMPLEX) $(PARALLEL)
	@echo "------------------------------ "
	$(FC_WORK) -o $(EXEC) $(OBJ) $(LIBS) $(LDFLAGS)

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi

#$(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90
#	@mkdir -p $(dir $@)
#	@echo "Compiling $<"
#	$(FC_WORK) -o $@ -c $<

$(OBJ_PATH)/%.F.o : $(SRC_PATH)/%.F
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<


$(OBJ_PATH)/%.F90.o : $(SRC_PATH)/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<


clean :
	@rm -rf $(OBJ_PATH)
	@rm -f $(EXEC)
	@rm -f *~
	@rm *.mod
