# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= /path/to/BoxLib

NDEBUG    := t
MPI       :=
OMP       :=
PROF      :=
COMP      := ifort
MKVERBOSE := t
EXECUTABLE:=spacemarine

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak

f90sources += physics_declarations.f90 \
              solvers.f90 \
              source_module.f90 \
              advance.f90 \
              init_phi.f90 \
              main.f90 \
              regrid.f90 \
              file_io.f90

VPATH_LOCATIONS += .

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

$(EXECUTABLE): $(objects) 
	$(LINK.f90) -o $(EXECUTABLE) $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
