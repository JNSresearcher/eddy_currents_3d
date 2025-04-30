F90=gfortran
# F90=ifort

DO=./
DB=./
DS=./

DLL="C:\Program Files (x86)\Intel\oneAPI\mkl\2023.2.0\redist\intel64"

F=EC3D
EXE=$(F)

ifeq ($(F90),gfortran) 
Csw = -c  -O2  -finit-local-zero  -ffree-line-length-none  
#  -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising

Lsw =   $(DLL)\mkl_rt.2.dll  $(DLL)\mkl_core.2.dll $(DLL)\mkl_intel_thread.2.dll 
# Lsw =

OBJ =              \
$(DO)m_vxc2data.o \
$(DO)vxc2data.o   \
$(DO)solvers.o    \
$(DO)utilites.o   \
$(DO)$(F).o 

EXE: $(OBJ)
	$(F90) $(OBJ) $(Lsw)  -o $(DB)$(EXE) 

$(OBJ): $(DO)%.o: $(DS)%.f90
	$(F90) $(Csw) $< -o $@
endif

ifeq ($(F90),ifort) 
Csw = -c /QO2  /Qzero /Qmkl /heap-arrays:2000 /assume:buffered_io -nologo 
# /check:bounds 
Lsw = /Qmkl /heap-arrays:2000

OBJ =                \
$(DO)m_vxc2data.obj \
$(DO)vxc2data.obj   \
$(DO)solvers.obj    \
$(DO)utilites.obj   \
$(DO)$(F).obj 

EXE: $(OBJ)
	$(F90) $(OBJ) $(Lsw)  -o $(DB)$(EXE) 

$(OBJ): $(DO)%.obj: $(DS)%.f90
	$(F90) $(Csw) $< -o $@
endif

# Cleaning everything
clean:
ifeq ($(F90),gfortran) 
	del *.o; del *.mod; del *.exe
endif
ifeq ($(F90),ifort) 
	del *.obj; del *.mod; del *.exe
endif

cleanout:
	rmdir /s /Q out

