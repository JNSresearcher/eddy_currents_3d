F90=gfortran
# F90=ifort

DO=./
DB=./
DS=./

F=EC3D
EXE=$(F)

ifeq ($(F90),gfortran) 
Csw = -c  -O2  -finit-local-zero -ffree-line-length-none 
# -fbounds-check
Lsw = 

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
Csw = -c /QO2  /Qzero /heap-arrays:2000 /assume:buffered_io
# /check:bounds 
Lsw =  /heap-arrays:2000

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

