! 3D-Eddy Current Calculation Program
! 
! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )

PROGRAM EC3D

USE m_fparser
USE m_vxc2data
IMPLICIT NONE

! --------------internal working variables--------------!
! variables to measure calculation time
INTEGER:: k0,k1,k2
REAL(8)   T, Tcalc, Tsavedata

INTEGER nCells,     &              ! number of cells with vector field sources
        nCells0,    &              ! number of cells with scalar field sources
        nCellsGlob, &              ! total number of cells
        NcellsX, NcellsY, NcellsZ  ! number of cells with vector field sources along the axes "x", "y", "z"
 
INTEGER m1, m2,  i,j,k, L,LL, m, n, nl, nn,ios
INTEGER  kim ,kjm, kkm, kip,kjp,kkp, kkpp, kdz
INTEGER  nim,nijm,njm,njkm,nkm,nikm,nijkm, nip,njp,nkp
REAL(8) sxp,sxm,sx,dsx,   syp,sym,sy,dsy,    szp,szm,sz,dsz, a, s

! arrays for new cells of field sources as they move  new_nodes0(:),
INTEGER, ALLOCATABLE ::  new_nodesX(:), new_nodesY(:), new_nodesZ(:) 
INTEGER  nij, Lnew, jnew, inew, flag_shift, flag_move,  movestop(3)

! arrays of parameters of difference schemes
REAL(8), ALLOCATABLE :: QxschP(:), QyschP(:), QzschP(:), QxschM(:), QyschM(:),QzschM(:),QschP(:)
REAL(8) WschP

! an array of column references for each row and a list of columns for BCG method
! the compressed sparse row  (CSR) format storage
INTEGER, ALLOCATABLE::  irow(:), jcol(:)   
REAL(8), ALLOCATABLE :: valA(:)  !,matrix(:,:)  
INTEGER  num_nz                ! number of non-zero elements

! arrays for storing numbers of boundary cells
! for vector potential 
integer num_bndX, num_bndY, num_bndZ
integer, allocatable:: cel_bndX(:),  cel_bndY(:),  cel_bndZ(:)
! for scalar potential
integer num_bndUx, num_bndUy, num_bndUz
integer, allocatable::  cel_bndUx(:),  cel_bndUy(:),  cel_bndUz(:)

integer num_nzX, num_nzY, num_nzZ, num_nzU

! variables to control the output of calculation results
INTEGER n_new, Ntime, Npoint, Nout, Nprint, Nprint_display, Nout_display
CHARACTER (LEN = 4)  ch_sd

! Uaf - solution vector, Jaf(:) - right hand side vector
REAL(8),ALLOCATABLE ::  Uaf(:),   Jaf(:),  Jafbuf(:)   ! arrays for vector fields (dimension nCells)
                       ! Jafbuf(:) - buffer arrays for storing inertial sources of vector field
!---------------------------------------!
! input data from SUBROUTINE vxc2data   !
!---------------------------------------!
REAL(8) :: delta(3)         ! array grid spacing along X, Y and Z   
REAL(8) :: dt,   &          ! time step
           Time, &          ! stop time
           dtt              ! jump duration
INTEGER :: sdx,sdy,sdz      ! number of cells along X,  Y and  Z   
INTEGER :: nsub,     &      ! number of physical domains
           nsub_air, &      ! number of environment domains
           nsub_glob        ! total number of domains
INTEGER :: numfun,   &      ! number of functions for calculating external sources
           numMech          ! number of functions for calculating movements of external sources

INTEGER :: size_PHYS_C      ! number of domains with diffusion properties

! solv - character string corresponding to the name of methods: 'BCG' or 'SOR'
CHARACTER(LEN=3) :: solv 
REAL(8)      :: BND(3,2)  ! values  boundary conditions

! character string for boundary conditions on 6 faces: N - Neumann D - Dirichlet A - Absorption
CHARACTER(6) :: bound     

CHARACTER(16):: files      ! name for output files
REAL(8) ::  tolerance      ! convergence criterion
INTEGER ::  iter, itmax    ! actual number of iterations,  maximum number of iterations  
!-------------------------------!
! call to data entry subroutine !
!-------------------------------!
CALL vxc2data ( delta,    dt,        Time,     dtt,                   &
                sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, size_PHYS_C,&
                numfun, numMech,       &
                solv,   files,    tolerance, itmax,    bound,  BND ) 
!----------------------!
! START of calculation !
!----------------------!
nCells = sdx*sdy*sdz; kdz=sdx*sdy
sz = 1.d0/(delta(3)*delta(3))
sy = 1.d0/(delta(2)*delta(2))
sx = 1.d0/(delta(1)*delta(1)) 
dsx = 0.5d0/delta(1) 
dsy = 0.5d0/delta(2)
dsz = 0.5d0/delta(3)

! calculated of number of cells with scalar field
Ncells0 = 0
do n=1,size_PHYS_C
    Ncells0 = Ncells0 + PHYS_C(n)%siznod
enddo
nCellsGlob = 3*nCells + Ncells0

!---------------------------------------------------------------------------!
! formation of a sparse matrix for the BCG solver                       !
! (BiConjugate Gradient stabilized method with restert)                     !
!---------------------------------------------------------------------------!
allocate (irow(nCellsGlob + 1),  source=1)
print*, 'nCells=',nCells, ' nCellsGlob=',nCellsGlob, ' Ncells0=',Ncells0, ' numfun=' ,numfun

CALL gen_sparse_matrix         ! generation of sparse matrix

    ! allocate(matrix(nCellsGlob,nCellsGlob), source=0.0d0 )
    ! do i=1,nCellsGlob
        ! do m = irow(i), irow(i+1)-1
            ! j = jcol(m)
            ! matrix(i,j) = valA(m)
        ! enddo
    ! enddo

    ! print*,'matrix:'
    ! write(*, '(  *(i12,a)  )' )  ( i,',',i=1,nCellsGlob)
    ! do i=1,nCellsGlob
        ! write(*,'(i3, a,$ )')  i, ' [,'
        ! ! write(*,'( *(g12.3) ,a)') ( matrix(i,j), j=1,Ncells0),'],'
        ! write(*,'( *(g12.3, a) ,a)') ( matrix(i,j),',' , j=1,nCellsGlob),'],'
    ! enddo
    ! stop

!-------------------------------!
! preparation of time intervals !
!-------------------------------!
T=0.
Ntime=0                              ! counter of each point with step dt
Nout_display=nint(Time/(100.0*DT))   ! jump size for display of symbol  ">"
Nprint_display=Ntime + Nout_display  ! counter of points for display symbol ">" with jump Nout_display
Npoint=0                             ! counter of calculated points with step dtt
Nout=nint(DTT/DT)                    ! jump size for output to files
Nprint = Ntime + Nout                ! point counter with step Nout for output to files with step Nout

!-----------------------------------!
! preparing arrays of vector fields !
!-----------------------------------!
allocate (Jaf(nCellsGlob), Uaf(nCellsGlob), source=0.0_8)
   
IF ( size_PHYS_C /= 0 ) THEN
    allocate (Jafbuf(nCellsGlob), source=0.0_8)
ENDIF

!---------------------------------------------------------------!
! preparation of arrays for calculating the movement of sources !
!---------------------------------------------------------------!
flag_move = 0 ! global flag of motion in input data

IF (numfun /=0 ) THEN
    DO i=1,numfun
    ! motion search processing
        fun_nod(i)%Distance = 0.d0 
    
        ! unlock
        IF (fun_nod(i)% num_Vmech(1) == 0  .and. fun_nod(i)%move(1) /=0) THEN ! it is not function
        !  distance traveled in time dt with speed Vx in fractions of cell size along dX
            fun_nod(i)%shift(1) = fun_nod(i)%vel_Vmech(1)*dt/delta(1)   !Vx*dt/dX
            flag_move = 1
        ELSEIF (  fun_nod(i)%move(1) /=0) THEN 
            flag_move = 1
        ENDIF
        
        IF (fun_nod(i)% num_Vmech(2) == 0  .and. fun_nod(i)%move(2) /=0) THEN ! it is not function
            fun_nod(i)%shift(2) = fun_nod(i)%vel_Vmech(2)*dt/delta(2)
            flag_move = 1; 
        ELSEIF (  fun_nod(i)%move(2) /=0) THEN 
            flag_move = 1
        ENDIF
    
        IF (fun_nod(i)% num_Vmech(3) == 0 .and. fun_nod(i)%move(3) /=0 ) THEN ! it is not function, but not a shift either
            fun_nod(i)%shift(3) = fun_nod(i)%vel_Vmech(3)*dt/delta(3)
            flag_move = 1; 
        ELSEIF (  fun_nod(i)%move(3) /=0) THEN 
            flag_move = 1
        ENDIF
    ENDDO

! initial allocation of array memory for source nodes, regardless of whether there is movement or not
    NcellsX = 0; NcellsY = 0; NcellsZ = 0; 
    DO n=1,numfun
        IF     (Fun(n)%ex == 'X') THEN
            NcellsX = NcellsX + fun_nod(n)%numnod_Fx
        ELSEIF (Fun(n)%ex == 'Y') THEN
            NcellsY = NcellsY + fun_nod(n)%numnod_Fy
        ELSEIF (Fun(n)%ex == 'Z') THEN
            NcellsZ = NcellsZ + fun_nod(n)%numnod_Fz
        ENDIF 
    ENDDO
    IF (NcellsX /= 0) ALLOCATE (new_nodesX(NcellsX), source=0)
    IF (NcellsY /= 0) ALLOCATE (new_nodesY(NcellsY), source=0)
    IF (NcellsZ /= 0) ALLOCATE (new_nodesZ(NcellsZ), source=0)
ENDIF !numfun

! filling source node arrays if there is no movement
IF (flag_move == 0 .and. numfun /=0 ) THEN
    NcellsX = 0; NcellsY = 0; NcellsZ = 0; ! counters
    DO n=1,numfun
        IF (Fun(n)%ex == 'X') THEN
            DO k = 1,fun_nod(n)%numnod_Fx
                m = fun_nod(n)%nods_Fx(k)
                NcellsX = NcellsX + 1
                new_nodesX(NcellsX) = m
            ENDDO
        ELSEIF (Fun(n)%ex == 'Y') THEN
            DO k = 1,fun_nod(n)%numnod_Fy
                m = fun_nod(n)%nods_Fy(k) - nCells
                NcellsY = NcellsY + 1
                new_nodesY(NcellsY) = m !+ 
            ENDDO
        ELSEIF (Fun(n)%ex == 'Z') THEN
            DO k = 1,fun_nod(n)%numnod_Fz
                m = fun_nod(n)%nods_Fz(k) - 2*nCells
                NcellsZ = NcellsZ + 1
                new_nodesZ(NcellsZ) = m !+ 
            ENDDO
        ELSE
            STOP 'err calc filling  array src, if not move'
        ENDIF  
    ENDDO
ENDIF

call execute_command_line ('mkdir '//trim(files), exitstat=i)
PRINT *, 'created dir '//trim(files)//': status=', i

!!!!!!  START time counter
CALL system_clock(k1,k0)
Tsavedata = 0.0
movestop = 1  ! flag of local start-stop movement in all coordinates
write(*,'( "start solver: ", a)') solv

2000 CONTINUE   
!-------------------------------------!
! Calculation of independent sources  !
!-------------------------------------!
IF (numfun /= 0) THEN
    CALL initf (numfun) 
    DO i=1,numfun
        DO m=1,Fun(i)%args
            IF ( trim(Fun(i)%namex(m)) == 'T') THEN
                Fun(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(fun(i)%eqn), fun(i)%namex(1:fun(i)%args) )
        fun(i)%vely = evalf (i, fun(i)%velx(1:fun(i)%args)) * 0.12566370964050292d-5  ! * mu0
    ENDDO   
ENDIF
!-------------------------------------------------!
! Calculation of the speed of movement of sources !
!-------------------------------------------------!
IF (numMech /= 0) THEN
    CALL initf (numMech) 
    DO i=1,numMech
        DO m=1,Vmech(i)%args
            IF ( trim(Vmech(i)%namex(m)) == 'T') THEN
                Vmech(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(Vmech(i)%eqn), Vmech(i)%namex(1:Vmech(i)%args) )
        Vmech(i)%vely = evalf (i, Vmech(i)%velx(1:Vmech(i)%args))
    ENDDO
ENDIF
!--------------------------------------------------!
! Distribution of independent sources into domains !
!--------------------------------------------------!
IF (flag_move == 1) THEN 
    ! there are moving sources ! 
    IF (size_PHYS_C /= 0  ) THEN      ! only inertial sources 
        Jafbuf = 0.0d0
        ! recording inertial sources to the buffer
        DO m=1,size_PHYS_C 
            ! n = PHYS_C(m)%siznod
            DO n=1, PHYS_C(m)%siznod
                L  = PHYS_C(m)%nod(n)  ! nodesL(m)  
                k  = L + nCells
                nl = L + 2*nCells
                Jafbuf(L) =  Jaf(L) 
                Jafbuf(k) =  Jaf(k)
                Jafbuf(nl) =  Jaf(nl)
            ENDDO
        enddo
        Jaf = 0.0d0  ! reset all sources
        Jaf = Jafbuf  ! restore of inertial sources
    ELSE
        Jaf = 0.0d0
    ENDIF
    
    ! calculation of new coordinates of nodes of independent sources 
    IF (numfun /= 0) THEN
        NcellsX=0; NcellsY=0; NcellsZ=0
        
        DO n=1,numfun
            CALL motion_calc    ! motion calculation
        
            a = Fun(n)%vely
            IF (Fun(n)%ex == 'X') THEN
                DO k = 1,fun_nod(n)%numnod_Fx
                    m = fun_nod(n)%nods_Fx(k)
                    CALL new_m !(m)
                    Jaf(m) = a
                    
                    NcellsX = NcellsX + 1
                    new_nodesX(NcellsX) = m
                ENDDO
            ELSEIF (Fun(n)%ex == 'Y') THEN
                DO k = 1,fun_nod(n)%numnod_Fy
                    m = fun_nod(n)%nods_Fy(k) - nCells
                    CALL new_m 
                    
                    NcellsY = NcellsY + 1
                    new_nodesY(NcellsY) = m 
                    
                    m = m + nCells
                    Jaf(m) = a
                ENDDO
            ELSEIF (Fun(n)%ex == 'Z') THEN
                DO k = 1,fun_nod(n)%numnod_Fz
                    m = fun_nod(n)%nods_Fz(k) - 2*nCells
                    CALL new_m 
                    
                    NcellsZ = NcellsZ + 1
                    new_nodesZ(NcellsZ) = m 
                    
                    m = m + 2*nCells
                    Jaf(m) = a
                ENDDO
            ELSE
                STOP 'err calc Fx or Fy or Fz '
            ENDIF                    
        ENDDO
    ENDIF
ELSE
    !--------------------!
    ! no moving sources  !
    !--------------------!
    DO n=1,numfun
        a = Fun(n)%vely
        IF (Fun(n)%ex == 'X') THEN
            DO k = 1,fun_nod(n)%numnod_Fx
                m = fun_nod(n)%nods_Fx(k)
                Jaf(m) = a
            ENDDO
        ELSEIF (Fun(n)%ex == 'Y') THEN
            DO k = 1,fun_nod(n)%numnod_Fy
                m = fun_nod(n)%nods_Fy(k)
                Jaf(m) = a
            ENDDO
        ELSEIF (Fun(n)%ex == 'Z') THEN
            DO k = 1,fun_nod(n)%numnod_Fz
                m = fun_nod(n)%nods_Fz(k)
                Jaf(m) = a
            ENDDO
        ELSE
            STOP 'err calc Fx or Fy or Fz '
        ENDIF  
    ENDDO

ENDIF


IF (size_PHYS_C /=0) THEN
    !---------------------------------!
    ! calculation of inertial sources !
    !---------------------------------!
    DO m=1,size_PHYS_C !kol_nL !
        a = PHYS_C(m)%valdom   ! 2.d0*valPHYS(kp,2)/dt   2*mu0*sigma/dt
        DO n=1, PHYS_C(m)%siznod
            L  = PHYS_C(m)%nod(n)  ! nodesL(m)  
            k  = L + nCells
            nl = L + 2*nCells
            
            Jaf(L)  = a * Uaf(L)  + Jaf(L) 
            Jaf(k)  = a * Uaf(k)  + Jaf(k)   ! Jb = 2С/dt * U_ + Ib  
            Jaf(nl) = a * Uaf(nl) + Jaf(nl)
            
            s=0.d0
            do i = irow(3*nCells + n), irow(3*nCells + n +1)-1   
                k = jcol(i)  !
                if ( k < 3*nCells+1) then !
                    s = s + valA(i) * Uaf(k) 
                endif
            enddo
            Jaf( 3*nCells + n ) =  s
        ENDDO
    ENDDO
    
    Jaf(cel_bndUx) = 0.0d0    
    Jaf(cel_bndUy) = 0.0d0 
    Jaf(cel_bndUz) = 0.0d0    
    
    Jaf(cel_bndX) = 0.0d0    
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0  
    
ENDIF
        !---------------------------!
        ! SOLVERS for vector fields !
        !---------------------------!
    CALL sprsBCGstabwr (valA, irow, jcol, nCellsGlob, Jaf, Uaf,  tolerance, itmax, iter)
        !---------------------------------------------!
        ! calculation of sources of inertial branches !
        !---------------------------------------------!
IF (size_PHYS_C /=0) THEN 
    DO m=1,size_PHYS_C !kol_nL !  PHYS_C( )%valdom  = 2.d0*valPHYS(kp,2)/dt
        a = PHYS_C(m)%valdom   ! 2.d0*valPHYS(kp,2)/dt   2*mu0*sigma/dt
        DO n=1, PHYS_C(m)%siznod
            L  = PHYS_C(m)%nod(n)  ! nodesL(m)  
            k  = L + nCells
            nl = L + 2*nCells
            
            Jaf(L)  = a * Uaf(L)  - Jaf(L) 
            Jaf(k)  = a * Uaf(k)  - Jaf(k)   ! Jb = 2С/dt * U_ + Ib  
            Jaf(nl) = a * Uaf(nl) - Jaf(nl)
        enddo
    enddo
        
    Jaf(cel_bndX) = 0.0d0    
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0   

    Uaf(cel_bndX) = 0.0d0    
    Uaf(cel_bndY) = 0.0d0 
    Uaf(cel_bndZ) = 0.0d0 
endif

!===========================================================================================
IF (Ntime >= Nprint .and. Ntime /=0 ) THEN  ! output with step dtt skip 1st point
    Nprint = Ntime  + Nout
    Npoint= Npoint + 1                      ! point counter with step dtt
    ios=0

    call writeVtk_field (Npoint,  sdx, sdy, sdz, nCells, delta, Uaf, Jaf, geoPHYS,geoPHYS_C, typPHYS, size_PHYS_C, files)
    CALL writeVtk_src ( Npoint, numfun,  NcellsX, NcellsY, NcellsZ,  new_nodesX, new_nodesY,  new_nodesZ,  &
                                         sdx, sdy, sdz, delta, files)
ENDIF

IF (Ntime >= Nprint_display) THEN
    Nprint_display = Ntime + Nout_display

    write(*,'( a,$ )') '>'
ENDIF

Ntime = Ntime + 1   ! this is every point
T = T + DT

IF (T < Time) GOTO 2000
PRINT*,'|'

CALL system_clock(k2);
Tcalc =  REAL(k2-k1)/REAL(k0)
   
write(*,'( "solve complet. Tcalc= ", g10.3)') Tcalc

CONTAINS

SUBROUTINE gen_sparse_matrix

integer colX(10), colY(10), colZ(10), colU(13)
real(8) valX(10), valY(10), valZ(10), valU(13)

real(8) s, a, b, sz, sy, sx, dsx, dsy, dsz

integer countU, nAx, nAy, nAz, Lx, Ly, Lz, Lfi, nFix, nFiy, nFiz,  kFi,  nc

type espm
    integer im
    DOUBLE PRECISION em
    type (espm),pointer ::prec
end type espm
type(espm),pointer :: sp_colX,next_X, sp_colY,next_Y,  sp_colZ,next_Z,   sp_colU,next_U  !, sp_nL, next_nL

nullify(sp_colX, next_X,sp_colY, next_Y,sp_colZ, next_Z,sp_colU,next_U ) !, sp_F, next_F )

if (size_PHYS_C /= 0) then
    num_bndX=0; num_bndY=0; num_bndZ=0
    cel_bndX = [integer:: ];  cel_bndY = [integer:: ];  cel_bndZ = [integer:: ] 

    num_bndUx=0; num_bndUy=0; num_bndUz=0
    cel_bndUx = [integer:: ];  cel_bndUy = [integer:: ];  cel_bndUz = [integer:: ]
endif

num_nzX=0;  !irow(        1) = 1;  
num_nzY=0;  !irow(  nCells+1) = 1;  
num_nzZ=0;  !irow(2*nCells+1) = 1;
num_nzU=0;  !irow(3*nCells+1) = 1;  

sz = 1.d0/(delta(3)*delta(3))
sy = 1.d0/(delta(2)*delta(2))
sx = 1.d0/(delta(1)*delta(1)) 
dsx = 0.5d0/delta(1)  ! 1/2dx
dsy = 0.5d0/delta(2) ! 1/2dy
dsz = 0.5d0/delta(3) ! 1/2dy

!  X-1:nCells  Y-nCells+1: 2*nCells  Fi-2*nCells+1 : nCellsGlob (2*nCells + nCellsFi)
nn = 0; 
countU = 0 
do k=1,sdz
  do j=1,sdy
    do i=1,sdx
        n = geoPHYS(i,j,k)   
        nn =  nn + 1    
        valX = 0.0d0;   colX = 0;  Lx=0; 
        valY = 0.0d0;   colY = 0;  Ly=0; 
        valZ = 0.0d0;   colZ = 0;  Lz=0;
        nAx=0; nAy=0; nAz=0  
        valU = 0.0d0;  colU = 0; Lfi=0; 
        nFix=0; nFiy=0; nFiz=0  
        kFi=0;

        IF (geoPHYS_C(i,j,k) /=0) THEN
            kFi=1;                           
            countU = countU + 1 
        endif

        kim = nn-1;  kjm = nn-sdx;  kkm = nn-kdz;
        kip = nn+1;  kjp = nn+sdx;  kkp = nn+kdz;
        ! 8 corners
        
    IF (i==1 .or.j==1 .or.k==1 .or.i==sdx .or.j==sdy .or.k==sdz ) THEN
!boundary X
        if     (i==1 .and. j==1 .and. k==1  ) then
            Lx = 4; 
            colX(1:Lx) = [kip, kjp, kkp,     nn]
            valX(1:Lx) = [BND(1,2)*sx, BND(2,2)*sy ,BND(3,2)*sz, (sx + sy + sz)  ]
        elseif ( i==sdx .and. j==1 .and. k==1) then
            Lx = 4; 
            colX(1:Lx) = [kim, kjp, kkp,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, BND(2,2)*sy ,BND(3,2)*sz, (sx + sy + sz)  ]
        elseif (i==1 .and. j==sdy .and. k==1) then
            Lx = 4; 
            colX(1:Lx) = [kip, kjm, kkp,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, BND(2,1)*sy ,BND(3,2)*sz, (sx + sy + sz)  ]
        elseif (i==sdx .and. j==sdy .and. k==1 ) then
            Lx = 4; 
            colX(1:Lx) = [kim, kjm, kkp,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, BND(2,1)*sy ,BND(3,2)*sz, (sx + sy + sz)  ]
        elseif ( i==1 .and. j==1   .and. k==sdz) then
            Lx = 4; 
            colX(1:Lx) = [kip, kjp, kkm,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, BND(2,2)*sy ,BND(3,1)*sz, (sx + sy + sz)  ]
        elseif (i==sdx .and. j==1   .and. k==sdz ) then
            Lx = 4; 
            colX(1:Lx) = [kim, kjp, kkm,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, BND(2,2)*sy ,BND(3,1)*sz, (sx + sy + sz)  ]
        elseif (i==1   .and. j==sdy .and. k==sdz) then
            Lx = 4; 
            colX(1:Lx) = [kip, kjm, kkm,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, BND(2,1)*sy ,BND(3,1)*sz, (sx + sy + sz)  ]
        elseif (i==sdx .and. j==sdy .and. k==sdz) then
            Lx = 4; 
            colX(1:Lx) = [kim, kjm, kkm,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, BND(2,1)*sy ,BND(3,1)*sz, (sx + sy + sz)  ]
            
            !12 edges
            !edge along  x
        elseif ((i>1 .and. i<sdx ) .and. j==1 .and. k==1 ) then
            Lx = 5; 
            colX(1:Lx) = [kim, kip, kjp, kkp,     nn]
            valX(1:Lx) =   [-sx, -sx, BND(2,2)*sy ,BND(3,2)*sz, (2.d0*sx + sy + sz)  ]
        elseif ((i>1 .and. i<sdx ) .and. j==sdy .and.  k==1 ) then
            Lx = 5; 
            colX(1:Lx) = [kim, kip, kjm, kkp,     nn]
            valX(1:Lx) =   [-sx, -sx, BND(2,1)*sy ,BND(3,2)*sz, (2.d0*sx + sy + sz)  ]
        elseif ((i>1 .and. i<sdx ) .and. j==1 .and. k==sdz  ) then
            Lx = 5; 
            colX(1:Lx) = [kim, kip, kjp, kkm,     nn]
            valX(1:Lx) =   [-sx, -sx, BND(2,2)*sy ,BND(3,1)*sz, (2.d0*sx + sy + sz)  ]
        elseif ((i>1 .and. i<sdx ) .and. j==sdy .and.  k==sdz) then
            Lx = 5; 
            colX(1:Lx) = [kim, kip, kjm, kkm,     nn]
            valX(1:Lx) =   [-sx, -sx, BND(2,1)*sy ,BND(3,1)*sz, (2.d0*sx + sy + sz)  ]
            
            !edge along y
        elseif (  i==1  .and. (j>1 .and. j<sdy ) .and. k==1) then
            Lx = 5; 
            colX(1:Lx) = [kip, kjm, kjp, kkp,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, -sy, -sy ,BND(3,2)*sz, (sx + 2.d0*sy + sz)  ]
        elseif ( i==sdx .and. (j>1 .and. j<sdy ) .and. k==1 ) then
            Lx = 5; 
            colX(1:Lx) = [kim, kjm, kjp, kkp,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, -sy, -sy ,BND(3,2)*sz, (sx + 2.d0*sy + sz)  ]
        elseif ( i==1   .and. (j>1 .and. j<sdy ) .and. k==sdz ) then
            Lx = 5; 
            colX(1:Lx) = [kip, kjm, kjp, kkm,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, -sy, -sy ,BND(3,1)*sz, (sx + 2.d0*sy + sz)  ]
        elseif ( i==sdx .and. (j>1 .and. j<sdy ) .and. k==sdz ) then
            Lx = 5; 
            colX(1:Lx) = [kim, kjm, kjp, kkm,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, -sy, -sy ,BND(3,1)*sz, (sx + 2.d0*sy + sz)  ]
            
            !edge along Z
        elseif ( i==1   .and. j==1   .and. (k>1 .and. k<sdz ) ) then
            Lx = 5; 
            colX(1:Lx) = [kip, kjp, kkm, kkp,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, BND(2,2)*sy, -sz, -sz, (sx + sy + 2.d0*sz)  ]
        elseif ( i==sdx .and. j==1   .and. (k>1 .and. k<sdz ) ) then
            Lx = 5; 
            colX(1:Lx) = [kim, kjp, kkm, kkp,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, BND(2,2)*sy, -sz, -sz, (sx + sy + 2.d0*sz)  ]
        elseif ( i==1   .and. j==sdy .and. (k>1 .and. k<sdz ) ) then
            Lx = 5; 
            colX(1:Lx) = [kip, kjm, kkm, kkp,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, BND(2,1)*sy, -sz, -sz, (sx + sy + 2.d0*sz)  ]
        elseif ( i==sdx .and. j==sdy .and. (k>1 .and. k<sdz ) ) then
            Lx = 5; 
            colX(1:Lx) = [kim, kjm, kkm, kkp,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, BND(2,1)*sy, -sz, -sz, (sx + sy + 2.d0*sz)  ]

            !6 faces
        elseif ( (i>1 .and. i<sdx ) .and. (j>1 .and. j<sdy ) .and. k==1 ) then ! face XY z
            Lx = 6; 
            colX(1:Lx) = [kim, kip, kjm, kjp, kkp,     nn]
            valX(1:Lx) =   [-sx, -sx, -sy, -sy ,BND(3,2)*sz, (2.d0*sx + 2.d0*sy + sz)  ]
        elseif ( (i>1 .and. i<sdx ) .and. (j>1 .and. j<sdy ) .and. k==sdz ) then
            Lx = 6; 
            colX(1:Lx) = [kim, kip, kjm, kjp, kkm,     nn]
            valX(1:Lx) =   [-sx, -sx, -sy, -sy ,BND(3,1)*sz, (2.d0*sx + 2.d0*sy + sz)  ]
        elseif (  i==1  .and. (j>1 .and. j<sdy ) .and. (k>1 .and. k<sdz ) ) then
            Lx = 6; 
            colX(1:Lx) = [kip, kjm, kjp, kkm, kkp,     nn]
            valX(1:Lx) =   [BND(1,2)*sx, -sy, -sy, -sz, -sz, (sx + 2.d0*sy + 2.d0*sz)  ]
        elseif ( i==sdx .and. (j>1 .and. j<sdy ) .and. (k>1 .and. k<sdz ) ) then
            Lx = 6; 
            colX(1:Lx) = [kim, kjm, kjp, kkm, kkp,     nn]
            valX(1:Lx) =   [BND(1,1)*sx, -sy, -sy, -sz, -sz, (sx + 2.d0*sy + 2.d0*sz)  ]
        elseif ( (i>1 .and. i<sdx ) .and. j==1   .and.(k>1 .and. k<sdz ) ) then
            Lx = 6; 
            colX(1:Lx) = [kim, kip, kjp, kkm, kkp,     nn]
            valX(1:Lx) =   [-sx, -sx, BND(2,2)*sy, -sz, -sz, (2.d0*sx + sy + 2.d0*sz)  ]
        elseif ( (i>1 .and. i<sdx ) .and. j==sdy .and. (k>1 .and. k<sdz ) ) then
            Lx = 6; 
            colX(1:Lx) = [kim, kip, kjm, kkm, kkp,     nn]
            valX(1:Lx) =   [-sx, -sx, BND(2,1)*sy, -sz, -sz, (2.d0*sx + sy + 2.d0*sz)  ]
        endif
        
        Ly = Lx; Lz = Lx
        colY = nCells+colX; colZ = 2*nCells+colX; valY = valX; valZ = valX
    else
!-----------------------------------  
        Lx=7 ! for Ax
        colX(1:Lx) =  [kim, kip, kjm, kjp, kkm, kkp,     nn]
        valX(1:Lx) =  [-sx,-sx, -sy, -sy, -sz, -sz, 2.d0*(sx + sy + sz)  ]
        
        Ly = Lx; Lz = Lx
        colY = nCells+colX; colZ = 2*nCells+colX; valY = valX; valZ = valX
        
        if (kFi /= 0) then ! add diffus component
            valX(1) = valX(1) - valPHYS(n,3)/(2.d0*delta(1)) !  Vex *valPHYS(n,2) valPHYS(n,3) = Vex='mu0*37e6*15'
            valX(2) = valX(2) + valPHYS(n,3)/(2.d0*delta(1)) 
            valX(3) = valX(3) - valPHYS(n,4)/(2.d0*delta(2)) !  Vey PHYS_C()%valdom * valPHYS(n,4
            valX(4) = valX(4) + valPHYS(n,4)/(2.d0*delta(2)) 
            valX(5) = valX(5) - valPHYS(n,5)/(2.d0*delta(3)) !  Vez
            valX(6) = valX(6) + valPHYS(n,5)/(2.d0*delta(3)) 
            valX(7) = valX(7) + 2.d0*valPHYS(n,2)/dt

            valY = valX; valZ = valX
                                ! PHYS_C( )%valdom  = 2.d0*valPHYS(kp,2)/dt
            if     (geoPHYS_C(i+1,j,k) == 0) then
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i,j,k);   valX(Lx) = -3.d0*valPHYS(n,2)*dsx
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i-1,j,k); valX(Lx) = +4.d0*valPHYS(n,2)*dsx
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i-2,j,k); valX(Lx) = -1.d0*valPHYS(n,2)*dsx
                nAx  = 1
            elseif (geoPHYS_C(i-1,j,k) == 0) then
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i,j,k);   valX(Lx) = +3.d0*valPHYS(n,2)*dsx
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i+1,j,k); valX(Lx) = -4.d0*valPHYS(n,2)*dsx
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i+2,j,k); valX(Lx) = +1.d0*valPHYS(n,2)*dsx
                nAx  = 1
            else
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i+1,j,k); valX(Lx) = -valPHYS(n,2)*dsx
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i-1,j,k); valX(Lx) = +valPHYS(n,2)*dsx
            endif

            if     (geoPHYS_C(i,j+1,k) == 0) then
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j,k);   valY(Ly) = -3.d0*valPHYS(n,2)*dsy
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j-1,k); valY(Ly) = +4.d0*valPHYS(n,2)*dsy
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j-2,k); valY(Ly) = -1.d0*valPHYS(n,2)*dsy
                nAy = 1
            elseif (geoPHYS_C(i,j-1,k) == 0) then
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j,k);   valY(Ly) = +3.d0*valPHYS(n,2)*dsy
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j+1,k); valY(Ly) = -4.d0*valPHYS(n,2)*dsy
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j+2,k); valY(Ly) = +1.d0*valPHYS(n,2)*dsy
                nAy = 1
            else
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j+1,k); valY(Ly) = -valPHYS(n,2)*dsy
                Ly = Ly + 1; colY(Ly) = geoPHYS_C(i,j-1,k); valY(Ly) = +valPHYS(n,2)*dsy
            endif
            
            if     (geoPHYS_C(i,j,k+1) == 0) then
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k);   valZ(Lz) = -3.d0*valPHYS(n,2)*dsz
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k-1); valZ(Lz) = +4.d0*valPHYS(n,2)*dsz
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k-2); valZ(Lz) = -1.d0*valPHYS(n,2)*dsz
                nAz = 1
            elseif (geoPHYS_C(i,j,k-1) == 0) then
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k);  valZ(Lz) = +3.d0*valPHYS(n,2)*dsz
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k+1); valZ(Lz) = -4.d0*valPHYS(n,2)*dsz
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k+2); valZ(Lz) = +1.d0*valPHYS(n,2)*dsz
                nAz = 1
            else
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k+1); valZ(Lz) = -valPHYS(n,2)*dsz
                Lz = Lz + 1; colZ(Lz) = geoPHYS_C(i,j,k-1); valZ(Lz) = +valPHYS(n,2)*dsz
            endif
        endif
        ! end form row for A
    endif
        !--------------------------X------------------------------------
        call full_sort(colX, valX, Lx, 1,1)
        do m = 1, Lx
            if (colX(m) <= 0) then
                print*, 'cell=',nn, 'm',m,'colX(m)',colX(m)
                stop
            endif
            allocate(sp_colX)
            sp_colX=espm(colX(m),valX(m),next_X) 
            num_nzX = num_nzX + 1
            next_X => sp_colX
        enddo
        irow(nn+1) = irow(nn) + Lx  !   
        
        !---------------------------------Y----------------------
        call full_sort(colY, valY, Ly, 1,1)
        do m = 1, Ly
            if (colY(m) <= 0) then
                print*, 'cell=',nn, 'm',m,'colY(m)',colY(m)
                stop
            endif
            allocate(sp_colY)
            sp_colY=espm(colY(m),valY(m),next_Y) 
            num_nzY = num_nzY + 1
            next_Y => sp_colY
        enddo

        irow( nCells + nn+1) = irow( nCells + nn) + Ly  
        
        !-------------------------------Z---------------------------------------
        call full_sort(colZ, valZ, Lz, 1,1)
        do m = 1, Lz
            if (colZ(m) <= 0) then
                print*, 'cell=',nn, 'm',m,'colZ(m)',colZ(m)
                stop
            endif
            allocate(sp_colZ)
            sp_colZ=espm(colZ(m),valZ(m),next_Z) 
            num_nzZ = num_nzZ + 1
            next_Z => sp_colZ
        enddo

        irow( 2*nCells + nn+1) = irow( 2*nCells + nn) + Lz  !  
        
        if ( nAx==1 )  cel_bndX = [cel_bndX, nn]
        if ( nAy==1 )  cel_bndY = [cel_bndY, (nn + nCells) ]
        if ( nAz==1 )  cel_bndZ = [cel_bndZ, (nn + 2*nCells) ]

            
!================================================
!===========END vect cells A
!=================================================
    if (kFi /=0) then  
        nc = geoPHYS_C(i,j,k) ! 3*nCells
        if (nc /= 0 ) then
            nim = geoPHYS_C(i-1,j,k);   nip = geoPHYS_C(i+1,j,k)
            njm = geoPHYS_C(i,j-1,k);   njp = geoPHYS_C(i,j+1,k)
            nkm = geoPHYS_C(i,j,k-1);   nkp = geoPHYS_C(i,j,k+1) 
            ! 8 corners
            if     ( nim == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i-1 j-1 k-1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nip,       njp,      nkp,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,-2.d0*sz,  s,          a,      b,       -2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i+1 j-1 k-1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nim,       njp,      nkp,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s,          a,      b,       -2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
            else if ( nim == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i-1 j+1 k-1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(1));    b = +2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nip,       njm,      nkp,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s,          a,      b,       -2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i+1 j+1 k-1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = +2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nim,       njm,      nkp,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s,          a,      b,       -2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nim == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i-1 j-1 k+1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nip,       njp,      nkm,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s,          a,      b,       +2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i+1 j-1 k+1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nim,       njp,      nkm,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s,          a,      b,       +2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nim == 0 .and. njp == 0 .and. nkp == 0 ) then ! not i-1 j+1 k+1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nip,       njm,      nkm,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s,          a,      b,       +2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njp == 0 .and. nkp == 0 ) then  ! not  i+1 j+1 k+1     
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = +2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nim,       njm,      nkm,     nc,        nn,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s,          a,      b,       +2.0d0/(dt*delta(3)) ]
                nFix = 1; nFiy = 1; nFiz = 1; 
                !                        12 edges
                ! edges  along X
            elseif ( njp == 0  .and. nkm == 0 ) then  ! not  j+1  k-1  
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(2));    b = -2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nip,       nim,      njm,       nkp,        nc,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-1.d0*sx, -sx,  -2.d0*sy,  -2.d0*sz,    s,      a,      b ]
                nFiy = 1; nFiz = 1; 
            elseif ( njm == 0  .and. nkm == 0 ) then   ! not  j-1  k-1   
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(2));    b = -2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nip,       nim,      njp,       nkp,        nc,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-sx, -sx,  -2.d0*sy,  -2.d0*sz,    s,      a,      b ]
                nFiy = 1; nFiz = 1; 
            elseif ( njp == 0  .and. nkp == 0 ) then   ! not  j+1  k+1
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(2));    b = +2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nip,       nim,      njm,       nkm,        nc,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-sx, -sx,  -2.d0*sy,  -2.d0*sz,    s,      a,      b ]
                nFiy = 1; nFiz = 1; 
            elseif ( njm == 0  .and. nkp == 0 ) then   ! not  j-1  k+1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(2));    b = +2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nip,       nim,      njp,       nkm,        nc,   nCells+nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-sx, -sx,  -2.d0*sy,  -2.d0*sz,    s,      a,      b ]
                nFiy = 1; nFiz = 1; 
                
                ! edges  along  Y
            elseif ( nip == 0  .and. nkm == 0 ) then  ! not  i+1  k-1  
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nim,       njm,      njp,       nkp,        nc,     nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -sy,  -sy,  -2.d0*sz,    s,      a,      b ]
                nFix = 1; nFiz = 1; 
            elseif ( nim == 0  .and. nkm == 0 ) then   ! not  i-1  k-1  
                Lfi = 7; s = 2.d0*(sx + sy + sz); a = -2.0d0/(dt*delta(1)); b = -2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nip,       njm,      njp,     nkp,      nc,    nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -sy,-sy, -2.d0*sz,  s,      a,      b ]
                nFix = 1; nFiz = 1; 
            elseif ( nip == 0  .and. nkp == 0 ) then  ! not  i+1  k+1  
                Lfi = 7; s = 2.d0*(sx + sy + sz); a = +2.0d0/(dt*delta(1)); b = +2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nim,       njm,      njp,     nkm,      nc,    nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -sy,-sy, -2.d0*sz,  s,      a,      b ]
                nFix = 1;  nFiz = 1; 
            elseif ( nim == 0  .and. nkp == 0 ) then   ! not  i-1  k+1   
                Lfi = 7; s = 2.d0*(sx + sy + sz); a = -2.0d0/(dt*delta(1)); b = +2.0d0/(dt*delta(3))
                colU(1:Lfi) = [nip,       njm,      njp,     nkm,      nc, nn, (2*nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -sy,-sy, -2.d0*sz,  s,   a,      b ]
                nFix = 1;  nFiz = 1; 

                ! edges  along Z 
            elseif ( nim == 0  .and. njm == 0 ) then   ! not  i-1 j-1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nip,       njp,      nkp,       nkm,        nc,   nn, (nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -sz,  -sz,    s,      a,      b ]
                nFix = 1; nFiy = 1; 
            elseif ( nip == 0 .and. njm == 0 ) then    ! not  i+1 j-1 ! 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = -2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nim,       njp,      nkp,       nkm,        nc,   nn, (nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -sz,  -sz,    s,      a,      b ]
                nFix = 1; nFiy = 1;  
            elseif ( nim == 0  .and. njp == 0 ) then   ! not  i-1 j+1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(1));    b = +2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nip,       njm,      nkp,       nkm,        nc,   nn, (nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -sz,  -sz,    s,      a,      b ]
                nFix = 1; nFiy = 1; 
            elseif ( nip == 0  .and. njp == 0 ) then   ! not   i+1 j+1  
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));    b = +2.0d0/(dt*delta(2))
                colU(1:Lfi) = [nim,       njm,      nkm,       nkp,        nc,   nn, (nCells + nn) ]
                valU(1:Lfi) =   [-2.d0*sx, -2.d0*sy,  -sz,  -sz,    s,    a,      b ]
                nFix = 1; nFiy = 1;  
                
                ! 6 faces  
            elseif (  nim == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                ! not   i-1 x-
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(1));   
                colU(1:Lfi) = [nip,       njm,      njp,       nkm,        nkp,        nc, nn ]
                valU(1:Lfi) =   [-2.d0*sx, -sy,  -sy,  -sz,  -sz,    s,   a ]
                nFix = 1; 
            elseif (  nip == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                 ! not   i+1
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(1));   
                colU(1:Lfi) = [nim,       njm,      njp,       nkm,        nkp,       nc, nn ]
                valU(1:Lfi) =   [-2.d0*sx, -sy,  -sy,  -sz,  -sz,    s,  a ]
                nFix = 1;  
            elseif (  njp == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                 ! not   j+1
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(2));   
                colU(1:Lfi) = [nim,       nip,        njm,       nkm,        nkp,    nc, nCells + nn ]
                valU(1:Lfi) =   [-sx, -sx,  -2.d0*sy,  -sz,  -sz,  s,   a ]
                nFiy = 1; 
            elseif (  njm == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                ! not   j-1
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(2));   
                colU(1:Lfi) = [nim,       nip,        njp,       nkm,        nkp,    nc, nCells + nn ]
                valU(1:Lfi) =   [-sx, -sx,  -2.d0*sy,  -sz,  -sz,  s,   a ]
                nFiy = 1;
            elseif (  nkp == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
                ! not   k+1 
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = +2.0d0/(dt*delta(3));   
                colU(1:Lfi) = [nim,       nip,        njm,       njp,        nkm,    nc, 2*nCells + nn ]
                valU(1:Lfi) =   [-sx, -sx,  -sy,  -sy,  -2.d0*sz,  s,   a ]
                nFiz = 1; 
            elseif (  nkm == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
                ! not   k-1
                Lfi = 7;       s = 2.d0*(sx + sy + sz);    a = -2.0d0/(dt*delta(3));   
                colU(1:Lfi) = [nim,       nip,        njm,       njp,        nkp,    nc, 2*nCells + nn ]
                valU(1:Lfi) =   [-sx, -sx,  -sy,  -sy,  -2.d0*sz,  s,   a ]
                nFiz = 1; 
            else
                Lfi = 13; s = 2.d0*(sx + sy + sz);
                colU(1:13) = [nim,  nip, njm, njp, nkm, nkp, nc, kip, kim, nCells+kjp, nCells+kjm, 2*nCells+kkp, 2*nCells+kkm ]
                valU(1:7)   = [-sx, -sx, -sy, -sy, -sz, -sz,  s ]
                valU(8:13) = 0.5d0/dt* [-1.d0/delta(1), 1.d0/delta(1), -1.d0/delta(2), 1.d0/delta(2), -1.d0/delta(3), 1.d0/delta(3)] 
            endif
            
            k0=0 
            mm: do k1=1,Lfi-1  
                do k2=k1+1, Lfi
                    if  (colU(k1) == colU (k2) ) then 
                        k0 = colU(k2)
                        exit mm
                    endif
                enddo
            enddo mm
            if (k0 /=0 ) then
                print*,'node Fi double', k0, 'i=',i, 'j=',j, 'k=',k
                stop
            endif
            
            if ( nFix == 1 ) cel_bndUx = [cel_bndUx, nc]
            if ( nFiy == 1 ) cel_bndUy = [cel_bndUy, nc]
            if ( nFiz == 1 ) cel_bndUz = [cel_bndUz, nc]
            !
            call full_sort(colU, valU, Lfi, 1,1)
            
            do m = 1, Lfi
                if (colU(m) <= 0) then
                    print*, 'cell=',nn, 'm',m,'colU(m)',colU(m)
                    stop
                endif
                allocate(sp_colU)
                sp_colU=espm(colU(m),valU(m),next_U) 
                num_nzU = num_nzU + 1 
                next_U => sp_colU
            enddo
            
            irow( 3*nCells + countU + 1) = irow( 3*nCells + countU) + Lfi
            
        endif !geoPHYS_C /=0
        
    endif ! kFi/=0  
            
    enddo !i
  enddo ! j
enddo ! k

print*,'end form Y'

num_nz = num_nzX + num_nzY + num_nzZ + num_nzU
print*, 'num_nzX=', num_nzX, 'num_nzY=', num_nzY, 'num_nzZ=', num_nzZ, 'num_nzU=', num_nzU

num_bndX = size(cel_bndX );  num_bndY = size(cel_bndY );  num_bndZ = size(cel_bndZ );  
print*, 'num_bndX=', num_bndX, 'num_bndY=', num_bndY, 'num_bndZ=', num_bndZ

irow(nCells+1) = num_nzX + 1
do i=nCells+2, 2*nCells
    irow(i) = irow(i) + num_nzX  !
enddo
m = num_nzX + num_nzY 
irow(2*nCells+1) = m + 1
do i=2*nCells+2, 3*nCells
    irow(i) = irow(i) + m  ! 
enddo
m =  num_nzX + num_nzY + num_nzZ
irow(3*nCells+1) = m + 1
do i=3*nCells+2, nCellsGlob+1                ! 
    irow(i) = irow(i) + m
enddo
! write(*, '(  a, *(i11) )')     '          row=', (m, m=1,nCellsGlob+1)      !ia
! write(*, '(  a, *(i11) )')     'irow row=', (irow(m),m=1,nCellsGlob+1)      !ia
! stop

allocate(jcol(num_nz), source=0)
allocate(valA(num_nz),source=0.0d0)
write(*,'( a,i9)')'start enter param, kol nonzero num_nz=',num_nz

i=num_nzX 
DO WHILE (associated(sp_colX))
    jcol(i) = sp_colX%im  !
    valA(i) = sp_colX%em   !
    i=i-1
    sp_colX => sp_colX%prec
END DO
nullify(sp_colX, next_X)

i=  num_nzX + num_nzY 
DO WHILE (associated(sp_colY))
    jcol(i) = sp_colY%im  
    valA(i) = sp_colY%em   !
    i=i-1
    sp_colY => sp_colY%prec
END DO
nullify(sp_colY, next_Y)

i=  num_nzX + num_nzY + num_nzZ
DO WHILE (associated(sp_colZ))
    jcol(i) = sp_colZ%im  
    valA(i) = sp_colZ%em   !
    i=i-1
    sp_colZ => sp_colZ%prec
END DO
nullify(sp_colZ, next_Z)
! tau = 18*0.005 = 0.09   2*tau*50 = 9 m/c
i=num_nz
DO WHILE (associated(sp_colU))
    jcol(i) = sp_colU%im  
    valA(i) = sp_colU%em   !
    i=i-1
    sp_colU => sp_colU%prec
END DO
nullify(sp_colU, next_U)

! write(*, '(  a, *(i11) )')     '          row=', (m, m=1,nCellsGlob+1)      !ia
! write(*, '(  a, *(i11) )')     'irow row=', (irow(m),m=1,nCellsGlob+1)      !ia

! write(*, '(  a )')     '          '
! write(*, '(  a, *(i11) )')     '          col=', (m, m=1,num_nz)   
! write(*, '(  a, *(i11) )')     '        jcol=', (jcol(m),m=1,num_nz)                 !ja
! write(*, '(  a, *(a,e10.3) )') '          val=', (' ',valA(m),m=1,num_nz)                 !aa
! stop

! do i=1,nCellsGlob
    ! write(*, '(  a, (i11) )')      '  i=', i    
    ! write(*, '(  a, *(i11) )')     '  j=', (jcol(j),    j=irow(i), irow(i+1)-1)
    ! write(*, '(  a, *(a,e10.3) )') 'val=', (' ',valA(j),j=irow(i), irow(i+1)-1)
! enddo

PRINT '(a,i3,a,i3,a,i3,     a,i9, a,g12.5,a/)', 'sparse matrix generation completed on grid (', sdx,' x ',sdy,' x ',sdz, &
                 ' ), Non zero elem= ',num_nz, ' Density of matrix:', 100.0* REAL(num_nz)/REAL(nCells)/REAL(nCells),'%'

END SUBROUTINE gen_sparse_matrix


SUBROUTINE motion_calc 
    DO i=1,3
        IF (fun_nod(n)% num_Vmech(i) == 0 ) THEN     !
            fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) + movestop(1)*fun_nod(n)%shift(i)
            fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i))
        ELSE                                         ! 
            fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) +  Vmech(fun_nod(n)%num_Vmech(i))%vely*dt/delta(i)
            fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i))
        ENDIF
    ENDDO
END SUBROUTINE motion_calc 

SUBROUTINE new_m 
    L = ceiling(REAL(m)/( REAL(sdx*sdy) ) ) 
    Lnew = L + fun_nod(n)%length(3)
    
    IF  ( Lnew > sdz-2  ) THEN 
        movestop(3) =0; Lnew = sdz-2
    ELSEIF(Lnew < 2    ) THEN
        movestop(3) =0;  Lnew = 2
    ELSEIF ( movestop(3) == 0 .and. (Lnew < sdz-2 .or. Lnew > 2)  ) THEN
             movestop(3) = 1 
    ENDIF
    IF (L == 1) THEN
        nij = m
    ELSE
        nij = m - (L-1)*sdx*sdy
    ENDIF
    j = ceiling( REAL(nij) / REAL(sdx) )  
    jnew = j + fun_nod(n)%length(2)
!-------------------------------------------------------------------------------------------
! comment/uncomment to check for out of bounds along the y-axis
!  1 variant
    IF  ( jnew > sdy-2  ) THEN 
        movestop(2) =0; jnew = sdy-2
    ELSEIF(jnew < 2    ) THEN
        movestop(2) =0;  jnew = 2
    ELSEIF ( movestop(2) == 0 .and. (jnew < sdy-2 .or. jnew > 2)  ) THEN
             movestop(2) = 1 
    ENDIF
!======================================================
!  2 variant
    ! IF  ( jnew > sdy  ) THEN 
        ! movestop(2) =0; jnew = sdy
    ! ELSEIF(jnew < 0    ) THEN
        ! movestop(2) =0;  jnew = 1
    ! ELSEIF ( movestop(2) == 0 .and. (jnew < sdy .or. jnew > 0 )  ) THEN
             ! movestop(2) = 1 
    ! ENDIF
!-----------------------------------------------------------------------------------------
    i = nij - (j - 1) * sdx  
    inew = i + fun_nod(n)%length(1)

    IF  ( inew > sdx-2  ) THEN  !
        movestop(1) =0;  inew = sdx-2
    ELSEIF(inew < 2    ) THEN !
        movestop(1) =0; inew = 2
    ELSEIF ( movestop(1) == 0 .and. (inew < sdx-2 .or. inew > 2)  ) THEN
             movestop(1) = 1 
    ENDIF

    m  = inew + sdx*(jnew-1) + sdx*sdy*(Lnew-1)
END SUBROUTINE new_m

END program EC3D


