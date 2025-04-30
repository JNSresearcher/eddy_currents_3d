! 3D-Eddy Current Calculation Program
! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! 
PROGRAM EC3D

USE m_fparser
USE m_vxc2data
IMPLICIT NONE

! --------------internal working variables--------------!
! variables to measure calculation time
INTEGER:: k0,k1,k2
REAL(8)   T, Tcalc, Tsavedata

INTEGER nCells,     &              ! number of cells with vector field sources
        nCellsGlob, &              ! total number of cells
        NcellsX, NcellsY, NcellsZ  ! number of cells with vector field sources along the axes "x", "y", "z"
 
INTEGER m1, m2,  i,j,k, L,LL, m, n, nl, nn,ios
INTEGER  kim ,kjm, kkm, kip,kjp,kkp,  kdz 
INTEGER  kim2 ,kjm2, kkm2, kip2,kjp2,kkp2
INTEGER  nim,njm,nkm,  nip,njp,nkp,  nijm
INTEGER  nim2,njm2,nkm2,  nip2,njp2,nkp2
INTEGER  kipjp, kipjm, kimjp, kimjm, kipkp, kipkm, kimkp, kimkm, kjpkp, kjpkm,  kjmkp, kjmkm

REAL(8) Vex,Ve_dx,sx,dsx,   Vey,Ve_dy,sy,dsy,    Vez,Ve_dz,sz,dsz, a, s
REAL(8) sxy, sxz, syz

! arrays for new cells of field sources as they move  new_nodes0(:),
INTEGER, ALLOCATABLE ::  new_nodesX(:), new_nodesY(:), new_nodesZ(:) 
INTEGER  nij, nX, nY, nZ, Lnew, jnew, inew, flag_move,  movestop(3)

!                               base                  for div                 for rhs
INTEGER, ALLOCATABLE::  irow(:), jcol(:),    irowU(:), jcolU(:),   irowUrhs(:), jcolUrhs(:)
REAL(8), ALLOCATABLE :: valA(:),             valAu(:),             valAUrhs(:) 
INTEGER                 num_nz,num_nzU,      num_nzU2,             num_nzU3   ! number of non-zero elements

! arrays for storing numbers of boundary cells
! for vector potential 
integer num_bndX, num_bndY, num_bndZ
integer, allocatable:: cel_bndX(:),  cel_bndY(:),  cel_bndZ(:)
! for scalar potential
integer num_bndUx, num_bndUy, num_bndUz
integer, allocatable::  cel_bndUx(:),  cel_bndUy(:),  cel_bndUz(:)

integer num_nzX, num_nzY, num_nzZ

! variables to control the output of calculation results
INTEGER  Ntime, Npoint, Nout, Nprint, Nprint_display, Nout_display

! Uaf - solution vector, Jaf(:) - right hand side vector
REAL(8),ALLOCATABLE ::  Uaf(:),   Jaf(:)  ! arrays for vector fields (dimension nCells)

! Jxsum - the sum of the speed, transformer and scalar components of the eddy current
REAL(8),ALLOCATABLE ::  Uf(:),   divJ(:), Jxsum(:), Jysum(:), Jzsum(:)
                       
!---------------------------------------!
! input data from SUBROUTINE vxc2data   !
!---------------------------------------!
REAL(8) :: delta(3)         ! array grid spacing along X, Y and Z   
REAL(8) :: dt,   &          ! time step
           Time, &          ! stop time
           dtt              ! jump duration

! for update rhs
REAL(8) :: valrhs(13), valU(40)
INTEGER :: colrhs(13), colU(40), Lfi,  Lx

INTEGER :: sdx,sdy,sdz      ! number of cells along X,  Y and  Z   
INTEGER :: nsub,     &      ! number of physical domains
           nsub_air, &      ! number of environment domains
           nsub_glob        ! total number of domains
INTEGER :: numfun,   &      ! number of functions for calculating external sources
           numMech,  &      ! number of functions for calculating movements of external sources
           numVenv

! solv - character string corresponding to the name of methods: 'BCG' or 'PRD'
CHARACTER(LEN=3) :: solv 
REAL(8)      :: BND(3,2)  ! values  boundary conditions

! character string for boundary conditions on 6 faces: N - Neumann D - Dirichlet A - Absorption
! CHARACTER(6) :: bound     ! not use

CHARACTER(16):: files      ! name for output files
REAL(8) ::  tolerance      ! convergence criterion
INTEGER ::  iter, itmax    ! actual number of iterations,  maximum number of iterations  

!// =====================the variables of pardiso====================
    integer(kind=8)           :: pt(64)   !// kind=8:x64; kind=4:win32
    integer(kind=8)           :: pt2(64)
    integer                   :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
    integer, allocatable      :: perm(:), perm2(:)
    integer                   :: iparm(64), iparm2(64)
!// =====================the variables of pardiso====================

!-------------------------------!
! call to data entry subroutine !
!-------------------------------!
CALL vxc2data ( delta,    dt,        Time,     dtt,                   &
                sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, &
                numfun, numMech,  numVenv,     &
                solv,   files,    tolerance, itmax,   BND ) 
!----------------------!
! START of calculation !
!----------------------!
nCells = sdx*sdy*sdz; kdz=sdx*sdy
sx = 1.d0/(delta(1)*delta(1)) 
sy = 1.d0/(delta(2)*delta(2))
sz = 1.d0/(delta(3)*delta(3))

dsx = 0.5d0/delta(1) ! 1/2dx
dsy = 0.5d0/delta(2) ! 1/2dy
dsz = 0.5d0/delta(3) ! 1/2dz

sxy = 1.d0/(delta(1)*delta(2)) 
sxz = 1.d0/(delta(1)*delta(3)) 
syz = 1.d0/(delta(2)*delta(3)) 

! calculated of number of cells with scalar field
nCellsGlob = 3*nCells + siznod_C

!----------------------------------------------------!
! formation of a sparse matrix                       !
!----------------------------------------------------!
allocate (irow(nCellsGlob + 1), irowU(siznod_C + 1), irowUrhs(siznod_C + 1), source=1)
! print'(*(a,i8,1x))', 'nCells=',nCells, ' nCellsGlob=',nCellsGlob, ' siznod_C=',siznod_C, ' numfun=' ,numfun
CALL gen_sparse_matrix         ! generation of sparse matrix

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
allocate (Jaf(nCellsGlob), Uaf(nCellsGlob),  source=0.0_8)
   
IF ( siznod_C /= 0 ) THEN
    allocate ( divJ(siznod_C), Uf(siznod_C), Jxsum(siznod_C),Jysum(siznod_C),Jzsum(siznod_C),  source=0.0_8)
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

!detect speed
nX=0;nY=0;nZ=0
IF (numMech /= 0) THEN
    Ve_dx=0.d0; Ve_dy=0.d0; Ve_dz=0.d0
    Vex=0.d0;  Vey=0.d0;  Vez=0.d0
    do n=1,numMech
        if     (Vmech(n)%ex == 'X' .and. nX==0 ) then 
            nX = n
        elseif (Vmech(n)%ex == 'Y' .and. nY==0 ) then 
            nY = n
        elseif (Vmech(n)%ex == 'Z' .and. nZ==0 ) then 
            nZ = n
        endif
    enddo
elseIF (numVenv /= 0) THEN
    DO n=1,numVenv
        if     (Venv(n)%ex == 'X' .and. nX==0 ) then 
            nX = n
        elseif (Venv(n)%ex == 'Y' .and. nY==0 ) then 
            nY = n
        elseif (Venv(n)%ex == 'Z' .and. nZ==0 ) then 
            nZ = n
        endif
    enddo
endif

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

if (solv == 'PRD') then
!// ===========================ini solve equations=============================
    allocate( perm( nCellsGlob ), perm2( siznod_C ),  source=0  )   !, permFi( NodesGlob )
    ! !allocate( perm( nCellsGlob ),  source=0  )   !, permFi( NodesGlob )
    pt = 0  !// pointer initialization
    maxfct = 1; mnum = 1; 
    ! !// mtype = -2: symmetric nonpositive definite matrix, mtype = 2: symmetric positive definite matrix
    mtype = 11  !real and unsymmetric matrix
    nrhs = 1
    iparm(1) = 0  !// iparm use default values
    iparm(5) = 1  
    iparm2(1) = 0  !// iparm use default values
    iparm2(5) = 1  
    error = 0  
    msglvl = 0
    phase = 12  !// LU decompose
    call pardiso( pt, maxfct, mnum, mtype, phase, nCellsGlob, valA, irow, jcol, perm, nrhs, iparm, msglvl, Jaf, Uaf, error )
    write(*, '(a, i9, a, i3)' ) 'pardiso for base: iparm15+16+63', iparm(15) + iparm(16) + iparm(63), ' error=',error
    call pardiso( pt2, maxfct, mnum, mtype, phase, siznod_C, valAu, irowU, jcolU, perm2, nrhs, iparm2, msglvl, divJ, Uf, error )
    write(*, '(a, i9, a, i3)' ) 'pardiso for div : iparm15+16+63', iparm2(15) + iparm2(16) + iparm2(63), ' error=',error
    phase = 33  !// solve equations
endif
!// ============================================================

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
!-------------------------------------------------!
! Calculation of the speed of movement of environment !
!-------------------------------------------------
IF (numVenv /= 0) THEN
    CALL initf (numVenv) 
    DO i=1,numVenv
        DO m=1,Venv(i)%args
            IF ( trim(Venv(i)%namex(m)) == 'T') THEN
                Venv(i)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (i, trim(Venv(i)%eqn), Venv(i)%namex(1:Venv(i)%args) )
        Venv(i)%vely = evalf (i, Venv(i)%velx(1:Venv(i)%args))
        
        if (Venv(i)%ex == 'X') then ! Vex 
            valPHYS( Venv(i)%nomsch,3)  = Venv(i)%vely * valPHYS( Venv(i)%nomsch,2)
        elseif (Venv(i)%ex == 'Y') then ! Vey 
            valPHYS( Venv(i)%nomsch,4)  = Venv(i)%vely * valPHYS( Venv(i)%nomsch,2)
        elseif (Venv(i)%ex == 'Z') then ! Vez
            valPHYS( Venv(i)%nomsch,5)  = Venv(i)%vely * valPHYS( Venv(i)%nomsch,2)
        endif
    ENDDO
ENDIF

!--------------------------------------------------!
! Distribution of independent sources into domains !
!--------------------------------------------------!
IF (flag_move == 1) THEN 
    ! there are moving sources ! 
    IF (siznod_C /= 0  ) THEN      ! only inertial sources 
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

IF (siznod_C /=0) THEN
    DO n=1, siznod_C
        L   = nods_C(n)  ! nodesL(m)  
        nn  = L + nCells
        nl  = L + 2*nCells
        ! calculation of inertial sources, valdom_C= 2.d0*valPHYS(m,2)/dt = 2*mu0*sigma/dt
        Jaf(L)  = valdom_C * Uaf(L)  + Jxsum(n) 
        Jaf(nn) = valdom_C * Uaf(nn) + Jysum(n)  ! Jb = 2C/dt * U_ + Ib  
        Jaf(nl) = valdom_C * Uaf(nl) + Jzsum(n) 
        m1=irowUrhs(n); m2=irowUrhs(n+1)-1 ! influence of a vector field on a scalar field
        Jaf( 3*nCells + n ) = dot_product(valAUrhs(m1:m2), Uaf(jcolUrhs(m1:m2) ) )
            ! alternativ variant without valAUrhs and jcolUrhs :
            ! s=0.d0
            ! do i = irow(3*nCells + n), irow(3*nCells + n +1)-1   
                ! k = jcol(i)
                ! if ( k < 3*nCells+1) then
                    ! s = s + valA(i) * Uaf(k) 
                ! endif
            ! enddo
            ! Jaf( 3*nCells + n ) =  s 
    enddo
    Jaf(cel_bndX) = 0.0d0 
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0  
 
!======================= update right parts of vector field=== and scalar field =========================
IF (numMech /= 0) THEN
    ! Ve_dx=0.d0; Ve_dy=0.d0; Ve_dz=0.d0
    ! Vex=0.d0;  Vey=0.d0;  Vez=0.d0
    ! Automatically sets Venv environment speed ( valPHYS(numdom_C,3) ) if field source speed ( Vmech(nX) ) is set
    if (nX /= 0 )  valPHYS(numdom_C,3) = -Vmech(nX)%vely*valPHYS(numdom_C,2) !  valPHYS(numdom_C,2) = sigma*mu0
    if (nY /= 0 )  valPHYS(numdom_C,4) = -Vmech(nY)%vely*valPHYS(numdom_C,2)
    if (nZ /= 0 )  valPHYS(numdom_C,5) = -Vmech(nZ)%vely*valPHYS(numdom_C,2)
endif

IF (numVenv /= 0  .or. numMech /= 0) THEN
!  constants for the velocity shift of the vector valrhs(1:7) and scalar valU field       
    if ( nx /=0) then
        Ve_dx=valPHYS(numdom_C,3)/(2.d0*delta(1)) ! valPHYS(m,3) = Vex * muo*sigma
        Vex = valPHYS(numdom_C,3)/valPHYS(numdom_C,2)    ! Vex = (Vex * muo*sigma)/muo*sigma
    endif
    if ( ny /=0) then
        Ve_dy=valPHYS(numdom_C,4)/(2.d0*delta(2))
        Vey = valPHYS(numdom_C,4)/valPHYS(numdom_C,2)
    endif
    if ( nz /=0) then
        Ve_dz=valPHYS(numdom_C,5)/(2.d0*delta(3))
        Vez = valPHYS(numdom_C,5)/valPHYS(numdom_C,2)
    endif

    DO n=1, siznod_C
        L  = nods_C(n)
        nn = L + nCells
        nl = L + 2*nCells
        
        k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
        IF (k == 1) THEN
            nij = L
        ELSE
            nij = L - (k-1)*sdx*sdy
        ENDIF
        j = ceiling( REAL(nij) / REAL(sdx) )  
        i = nij - (j - 1) * sdx  
        
        nim  = geoPHYS_C(i-1,j,k);   nip  = geoPHYS_C(i+1,j,k)
        njm  = geoPHYS_C(i,j-1,k);   njp  = geoPHYS_C(i,j+1,k)
        nkm  = geoPHYS_C(i,j,k-1);   nkp  = geoPHYS_C(i,j,k+1)
           
        kim = L-1;  kjm = L-sdx;  kkm = L-kdz;
        kip = L+1;  kjp = L+sdx;  kkp = L+kdz;

                 ! vector field stretching due to velocity
        if (  nim /= 0 .and. nip /= 0 .and. njm /= 0 .and. njp /= 0 .and. nkm /= 0 .and. nkp /= 0 ) then 
            Lx=0;  colrhs = 0;  valrhs = 0.d0 ! 
            if ( nx /=0 ) then
                colrhs(Lx+1:Lx+2) =       [kip,  kim ]
                valrhs(Lx+1:Lx+2) = Ve_dx*[-1d0, 1d0]
                Lx = Lx + 2  
            endif
            if ( ny /=0 ) then
                colrhs(Lx+1:Lx+2) =       [kjp, kjm  ]
                valrhs(Lx+1:Lx+2) = Ve_dy*[-1d0, 1d0]
                Lx = Lx + 2
            endif
            if ( nz /=0 ) then
                colrhs(Lx+1:Lx+2) =      [kkp, kkm  ]
                valrhs(Lx+1:Lx+2) = Ve_dz*[-1d0, 1d0]
                Lx = Lx + 2
            endif

            Jaf( colrhs(1:Lx)) = Jaf( colrhs(1:Lx) ) + valrhs(1:Lx) * Uaf( colrhs(1:Lx) )

            colrhs(1:Lx) = colrhs(1:Lx) + nCells
            Jaf( colrhs(1:Lx)) = Jaf( colrhs(1:Lx) ) + valrhs(1:Lx) * Uaf( colrhs(1:Lx) )

            colrhs(1:Lx) = colrhs(1:Lx) + nCells
            Jaf( colrhs(1:Lx)) = Jaf( colrhs(1:Lx) ) + valrhs(1:Lx) * Uaf( colrhs(1:Lx) )
        endif
        
        kim2= L-2;  kjm2 = L-2*sdx; kkm2 = L-2*kdz;
        kip2= L+2;  kjp2 = L+2*sdx; kkp2 = L+2*kdz; 
        kimjm = L - 1 - sdx; kimjp = L - 1 + sdx; kimkm = L - 1 - kdz;  kimkp = L - 1 + kdz;
        kipjm = L + 1 - sdx; kipjp = L + 1 + sdx; kipkm = L + 1 - kdz;  kipkp = L + 1 + kdz; 
        kjmkm = L - sdx - kdz;  kjmkp = L - sdx + kdz;
        kjpkm = L + sdx - kdz;  kjpkp = L + sdx + kdz;

        Lfi = 0; colU = 0; valU = 0.d0
        call update_rhs
        
        !adding velocity components of the vector field to the central node Vx*ddA/ddx of the scalar field
        Jaf( 3*nCells + n ) = Jaf( 3*nCells + n ) + dot_product(valU(1:Lfi), Uaf(colU(1:Lfi)) ) 
    ENDDO

ENDIF  ! numVenv /=0
ENDIF  ! siznod_C /=0

        !---------------------------!
        ! SOLVERS for vector fields !
        !---------------------------!
if (solv == 'BCG') then
    CALL sprsBCGstabwr (valA, irow, jcol, nCellsGlob, Jaf, Uaf,  tolerance, itmax, iter)
else
    call pardiso( pt, maxfct, mnum, mtype, phase, nCellsGlob, valA, irow, jcol, perm, nrhs, iparm, msglvl, Jaf, Uaf, error )
endif
        !---------------------------------------------!
        ! calculation of sources of inertial branches !
        !---------------------------------------------!
IF (siznod_C /=0) THEN 
    do n=1,siznod_C
        L  = nods_C(n)  ! nodesL(m)  
        nn = L + nCells
        nl = L + 2*nCells
        
        ! valdom_C  = 2.d0*valPHYS(kp,2)/dt  = 2*mu0*sigma/dt
        Jaf(L)  = valdom_C * Uaf(L)  - Jaf(L) 
        Jaf(nn) = valdom_C * Uaf(nn) - Jaf(nn)   ! Jb = 2C/dt * U_ + Ib  
        Jaf(nl) = valdom_C * Uaf(nl) - Jaf(nl)

        k= ceiling(REAL(L)/( REAL(sdx*sdy) ) )
        IF (k == 1) THEN
            nij = L
        ELSE
            nij = L - (k-1)*sdx*sdy
        ENDIF
        j = ceiling( REAL(nij) / REAL(sdx) )  
        i = nij - (j - 1) * sdx  
        ! a  = valPHYS(nijm,2)  !  valPHYS(nijm,2) = sigma*mu0
        ! a = PHYS_C(LL)%valdom  ! param C  domain a = 2*mu0*sigma/dt
        
        LL  = geoPHYS_C(i,j,k) !L + 3*nCells
        kim  = geoPHYS_C(i-1,j,k);   kip  = geoPHYS_C(i+1,j,k)
        kjm  = geoPHYS_C(i,j-1,k);   kjp  = geoPHYS_C(i,j+1,k)
        kkm  = geoPHYS_C(i,j,k-1);   kkp  = geoPHYS_C(i,j,k+1)
        kim2 = geoPHYS_C(i-2,j,k);   kip2 = geoPHYS_C(i+2,j,k)
        kjm2 = geoPHYS_C(i,j-2,k);   kjp2 = geoPHYS_C(i,j+2,k)
        kkm2 = geoPHYS_C(i,j,k-2);   kkp2 = geoPHYS_C(i,j,k+2)
        if     (kim == 0) then
            sx = 0.5d0*( -Uaf(kip2) + 4.d0*Uaf(kip) - 3.d0*Uaf(LL) )/delta(1)
        elseif (kip == 0) then
            sx = 0.5d0*(  Uaf(kim2) - 4.d0*Uaf(kim) + 3.d0*Uaf(LL) )/delta(1)
        else
            sx = 0.5d0*( Uaf(kip) - Uaf(kim) )/delta(1)
        endif
        
        if (kjm == 0) then
            sy = 0.5d0*( -Uaf(kjp2) + 4.d0*Uaf(kjp) - 3.d0*Uaf(LL) )/delta(2)
        elseif (kjp == 0) then
            sy = 0.5d0*(  Uaf(kjm2) - 4.d0*Uaf(kjm) + 3.d0*Uaf(LL) )/delta(2)
        else
            sy = 0.5d0*( Uaf(kjp) - Uaf(kjm) )/delta(2)
        endif
        
        if (kkm == 0) then
            sz = 0.5d0*( -Uaf(kkp2) + 4.d0*Uaf(kkp) - 3.d0*Uaf(LL) )/delta(3)
        elseif (kkp == 0) then
            sz = 0.5d0*(  Uaf(kkm2) - 4.d0*Uaf(kkm) + 3.d0*Uaf(LL) )/delta(3)
        else
            sz = 0.5d0*( Uaf(kkp) - Uaf(kkm) )/delta(3)
        endif
        ! current from scalar potential
        Jxsum(n) =  valPHYS(numdom_C,2)*sx
        Jysum(n) =  valPHYS(numdom_C,2)*sy
        Jzsum(n) =  valPHYS(numdom_C,2)*sz
    
    IF (numVenv /= 0  .or. numMech /= 0) THEN
        sx=0d0; sy=0d0; sz=0d0
        ! if (valPHYS(numdom_C,3) /= 0.d0) then
        if ( nx /=0) then
            nim = L - 1; ! njm  = L - sdx;   nkm  = L - kdz;
            nip = L + 1; ! njp  = L + sdx;   nkp  = L + kdz; 
            nim2= L - 2; ! njm2 = L-2*sdx;   nkm2 = L-2*kdz;
            nip2= L + 2; ! njp2 = L+2*sdx;   nkp2 = L+2*kdz; 
            if     (kim == 0) then
                sx = 0.5d0*( -Uaf(nip2) + 4.d0*Uaf(nip) - 3.d0*Uaf(L) )/delta(1)
                sy = 0.5d0*( -Uaf(nip2+nCells) + 4.d0*Uaf(nip+nCells) - 3.d0*Uaf(L+nCells) )/delta(1)
                sz = 0.5d0*( -Uaf(nip2+2*nCells) + 4.d0*Uaf(nip+2*nCells) - 3.d0*Uaf(L+2*nCells) )/delta(1)
            elseif (kip == 0) then
                sx = 0.5d0*(  Uaf(nim2) - 4.d0*Uaf(nim) + 3.d0*Uaf(L) )/delta(1)
                sy = 0.5d0*(  Uaf(nim2+nCells) - 4.d0*Uaf(nim+nCells) + 3.d0*Uaf(L+nCells) )/delta(1)
                sz = 0.5d0*(  Uaf(nim2+2*nCells) - 4.d0*Uaf(nim+2*nCells) + 3.d0*Uaf(L+2*nCells) )/delta(1)
            else
                sx = 0.5d0*( Uaf(nip)          - Uaf(nim) )        /delta(1)
                sy = 0.5d0*( Uaf(nCells+nip)   - Uaf(nCells+nim))  /delta(1)
                sz = 0.5d0*( Uaf(2*nCells+nip) - Uaf(2*nCells+nim))/delta(1)
            endif
            
            ! velocity component of current for Vx
            Jxsum(n) = Jxsum(n) +  sx*valPHYS(numdom_C,3)
            Jysum(n) = Jysum(n) +  sy*valPHYS(numdom_C,3)
            Jzsum(n) = Jzsum(n) +  sz*valPHYS(numdom_C,3) 
        endif
        
        sx=0d0; sy=0d0; sz=0d0
        ! if (valPHYS(numdom_C,4) /= 0.d0) then
        if ( ny /=0) then
            njm  = L - sdx;   !nkm  = L - kdz; nim = L - 1; 
            njp  = L + sdx;   !nkp  = L + kdz; nip = L + 1; 
            njm2 = L-2*sdx;   !nkm2 = L-2*kdz; nim2= L - 2; 
            njp2 = L+2*sdx;   !nkp2 = L+2*kdz; nip2= L + 2; 
            if     (kjm == 0) then
                sx = 0.5d0*( -Uaf(njp2) + 4.d0*Uaf(njp) - 3.d0*Uaf(L) )/delta(2)
                sy = 0.5d0*( -Uaf(njp2+nCells) + 4.d0*Uaf(njp+nCells) - 3.d0*Uaf(L+nCells) )/delta(2)
                sz = 0.5d0*( -Uaf(njp2+2*nCells) + 4.d0*Uaf(njp+2*nCells) - 3.d0*Uaf(L+2*nCells) )/delta(2)
            elseif (kjp == 0) then
                sx = 0.5d0*(  Uaf(njm2) - 4.d0*Uaf(njm) + 3.d0*Uaf(L) )/delta(2)
                sy = 0.5d0*(  Uaf(njm2+nCells) - 4.d0*Uaf(njm+nCells) + 3.d0*Uaf(L+nCells) )/delta(2)
                sz = 0.5d0*(  Uaf(njm2+2*nCells) - 4.d0*Uaf(njm+2*nCells) + 3.d0*Uaf(L+2*nCells) )/delta(2)
            else
                sx = 0.5d0*( Uaf(njp)          - Uaf(njm) )        /delta(2)
                sy = 0.5d0*( Uaf(nCells+njp)   - Uaf(nCells+njm))  /delta(2)
                sz = 0.5d0*( Uaf(2*nCells+njp) - Uaf(2*nCells+njm))/delta(2)
            endif
            
            ! velocity component of current for Vy
            Jxsum(n) = Jxsum(n) +  sx*valPHYS(numdom_C,4)
            Jysum(n) = Jysum(n) +  sy*valPHYS(numdom_C,4)
            Jzsum(n) = Jzsum(n) +  sz*valPHYS(numdom_C,4) 
        endif
        
        sx=0d0; sy=0d0; sz=0d0
        ! if (valPHYS(numdom_C,5) /= 0.d0) then
        if ( nz /=0) then
            nkm  = L - kdz; !nim = L - 1; njm  = L - sdx;   
            nkp  = L + kdz; !nip = L + 1; njp  = L + sdx;   
            nkm2 = L-2*kdz; !nim2= L - 2; njm2 = L-2*sdx;   
            nkp2 = L+2*kdz; !nip2= L + 2; njp2 = L+2*sdx;   
            if     (kkm == 0) then
                sx = 0.5d0*( -Uaf(nkp2) + 4.d0*Uaf(nkp) - 3.d0*Uaf(L) )/delta(3)
                sy = 0.5d0*( -Uaf(nkp2+nCells) + 4.d0*Uaf(nkp+nCells) - 3.d0*Uaf(L+nCells) )/delta(3)
                sz = 0.5d0*( -Uaf(nkp2+2*nCells) + 4.d0*Uaf(nkp+2*nCells) - 3.d0*Uaf(L+2*nCells) )/delta(3)
            elseif (kkp == 0) then
                sx = 0.5d0*(  Uaf(nkm2) - 4.d0*Uaf(nkm) + 3.d0*Uaf(L) )/delta(3)
                sy = 0.5d0*(  Uaf(nkm2+nCells) - 4.d0*Uaf(nkm+nCells) + 3.d0*Uaf(L+nCells) )/delta(3)
                sz = 0.5d0*(  Uaf(nkm2+2*nCells) - 4.d0*Uaf(nkm+2*nCells) + 3.d0*Uaf(L+2*nCells) )/delta(3)
            else
                sx = 0.5d0*( Uaf(nkp)          - Uaf(nkm) )        /delta(3)
                sy = 0.5d0*( Uaf(nCells+nkp)   - Uaf(nCells+nkm))  /delta(3)
                sz = 0.5d0*( Uaf(2*nCells+nkp) - Uaf(2*nCells+nkm))/delta(3)
            endif
            
            ! velocity component of current for Vz
            Jxsum(n) = Jxsum(n) +  sx*valPHYS(numdom_C,5)
            Jysum(n) = Jysum(n) +  sy*valPHYS(numdom_C,5)
            Jzsum(n) = Jzsum(n) +  sz*valPHYS(numdom_C,5) 
        endif
    endif

    enddo

    Jaf(cel_bndX) = 0.0d0 
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0  

    !=====================================================================
    !== cycle for Jaf!  Calculating Divergence Jaf
    !=====================================================================
    Uf=0.d0
    do nijm = 1,10
        DO n=1, siznod_C
        ! the global number of cell numbers for each np from the list PHYS_C(i)%numdom
            m  = nods_C(n)
            nn = m + nCells
            nl = m + 2*nCells
        
            k= ceiling(REAL(m)/( REAL(sdx*sdy) ) )
            IF (k == 1) THEN
                nij = m
            ELSE
                nij = m - (k-1)*sdx*sdy
            ENDIF
            j = ceiling( REAL(nij) / REAL(sdx) )  
            i = nij - (j - 1) * sdx  
            
            nim = m - 1;  njm  = m - sdx;   nkm  = m - kdz;
            nip = m + 1;  njp  = m + sdx;   nkp  = m + kdz; 

            kim = geoPHYS_C(i-1,j,k);   kip = geoPHYS_C(i+1,j,k)
            kjm = geoPHYS_C(i,j-1,k);   kjp = geoPHYS_C(i,j+1,k)
            kkm = geoPHYS_C(i,j,k-1);   kkp = geoPHYS_C(i,j,k+1) 
                       
            sx=0.5d0;sy=0.5d0;sz=0.5d0;
            if (kim == 0) then
                nim = m; sx=1.d0
            endif
            if (kip == 0)  then
                nip = m; sx=1.d0
            endif
            if (kjm == 0) then 
                njm = m; sy=1.d0
            endif
            if (kjp == 0) then 
                njp = m; sy=1.d0
            endif
            if (kkm == 0) then 
                nkm = m; sz=1.d0
            endif
            if (kkp == 0) then 
                nkp = m;  sz=1.d0
            endif
                        ! 0.5 - optim
            divJ(n) =  -0.5d0*( sx*( Jaf(nip) - Jaf(nim) )/delta(1) + &  ! dJx/dx + dJy/dy + dJz/dz
                          sy*( Jaf(njp+  nCells) - Jaf(njm+   nCells) )/delta(2) + &
                          sz*( Jaf(nkp + 2*nCells) - Jaf(nkm+ 2*nCells) )/delta(3) ) 
            
        enddo
        
        Uf(cel_bndUx- 3*nCells)=0.d0
        Uf(cel_bndUy- 3*nCells)=0.d0
        Uf(cel_bndUz- 3*nCells)=0.d0
        
        if (solv == 'BCG') then
            CALL sprsBCGstabwr (valAu, irowU, jcolU, siznod_C, divJ, Uf,  tolerance, itmax, iter)
        else
            call pardiso( pt2, maxfct, mnum, mtype, phase, siznod_C, valAu, irowU, jcolU, perm2, nrhs, iparm2, msglvl, divJ, Uf, error )
        endif
        
        DO n=1, siznod_C
            m  = nods_C(n)
            nn = m + nCells
            nl = m + 2*nCells
            
            k= ceiling(REAL(m)/( REAL(sdx*sdy) ) )
            IF (k == 1) THEN
                nij = m
            ELSE
                nij = m - (k-1)*sdx*sdy
            ENDIF
            j = ceiling( REAL(nij) / REAL(sdx) )  
            i = nij - (j - 1) * sdx 
                
            LL = geoPHYS_C(i,j,k)- 3*nCells ! 3*nCells
            nim  = geoPHYS_C(i-1,j,k)- 3*nCells;   nip  = geoPHYS_C(i+1,j,k)- 3*nCells
            njm  = geoPHYS_C(i,j-1,k)- 3*nCells;   njp  = geoPHYS_C(i,j+1,k)- 3*nCells
            nkm  = geoPHYS_C(i,j,k-1)- 3*nCells;   nkp  = geoPHYS_C(i,j,k+1)- 3*nCells
            
            sx=0.5d0;sy=0.5d0;sz=0.5d0;
            if (nim <= 0) then
                nim = LL;sx=1.d0
            endif
            if (nip <= 0) then
                nip = LL;sx=1.d0
            endif
            if (njm <= 0) then
                njm = LL;sy=1.d0
            endif
            if (njp <= 0) then
                njp = LL;sy=1.d0
            endif
            if (nkm <= 0) then
                nkm = LL;sz=1.d0
            endif
            if (nkp <= 0) then
                nkp = LL;sz=1.d0
            endif
            ! Subtracting the gradient of a scalar field
            Jaf(m)  =  Jaf(m)  - sx*(Uf(nip) - Uf(nim) )/delta(1)
            Jaf(nn) =  Jaf(nn) - sy*(Uf(njp) - Uf(njm) )/delta(2)
            Jaf(nl) =  Jaf(nl) - sz*(Uf(nkp) - Uf(nkm) )/delta(3)
            
        enddo ! end cycle n
    enddo ! end cycle nijm
    !=====================================================================
    !== end cycle for calc divJaf
    !=====================================================================
    
    ! adding current  eddy from scalar potential + eddy from speed component
    DO n=1, siznod_C
        m  = nods_C(n) 
        nn = m + nCells
        nl = m + 2*nCells
        Jxsum(n) = Jxsum(n) + Jaf(m)   !Jxsum - is the sum of transformer, scalar and speed currents
        Jysum(n) = Jysum(n) + Jaf(nn)
        Jzsum(n) = Jzsum(n) + Jaf(nl) 
    enddo
    Jxsum(cel_bndUx - 3*nCells)=0.d0
    Jysum(cel_bndUy - 3*nCells)=0.d0
    Jzsum(cel_bndUz - 3*nCells)=0.d0
    
    !=====================================================================
    !== cycle for calc divJsum
    !=====================================================================
    s=0.d0; !Uf=0.d0
    do nijm = 1,10 
        DO n=1, siznod_C
        ! the global number of cell numbers for each np from the list PHYS_C(i)%numdom
            m  = nods_C(n)
            
            k= ceiling(REAL(m)/( REAL(sdx*sdy) ) )
            IF (k == 1) THEN
                nij = m
            ELSE
                nij = m - (k-1)*sdx*sdy
            ENDIF
            j = ceiling( REAL(nij) / REAL(sdx) )  
            i = nij - (j - 1) * sdx  

            nim = geoPHYS_C(i-1,j,k)- 3*nCells;   nip = geoPHYS_C(i+1,j,k)- 3*nCells
            njm = geoPHYS_C(i,j-1,k)- 3*nCells;   njp = geoPHYS_C(i,j+1,k)- 3*nCells
            nkm = geoPHYS_C(i,j,k-1)- 3*nCells;   nkp = geoPHYS_C(i,j,k+1)- 3*nCells
            
            sx=0.5d0;sy=0.5d0;sz=0.5d0;
            if (nim <= 0) then
                nim = n;sx=1.d0
            endif
            if (nip <= 0) then
                nip = n;sx=1.d0
            endif
            if (njm <= 0) then
                njm = n;sy=1.d0
            endif
            if (njp <= 0) then
                njp = n;sy=1.d0
            endif
            if (nkm <= 0) then
                nkm = n;sz=1.d0
            endif
            if (nkp <= 0) then
                nkp = n;sz=1.d0
            endif
                    ! 0.5 - optim;  0.25 
            divJ(n) =  -0.5d0*( sx*( Jxsum(nip) - Jxsum(nim) )/delta(1) + &  ! dJx/dx + dJy/dy + dJz/dz
                            sy*( Jysum(njp) - Jysum(njm) )/delta(2) + &
                            sz*( Jzsum(nkp) - Jzsum(nkm) )/delta(3) ) 
        enddo
   
        Uf(cel_bndUx- 3*nCells)=0.d0
        Uf(cel_bndUy- 3*nCells)=0.d0
        Uf(cel_bndUz- 3*nCells)=0.d0

        if (solv == 'BCG') then
            CALL sprsBCGstabwr (valAu, irowU, jcolU, siznod_C, divJ, Uf,  tolerance, itmax, iter)
        else
            call pardiso( pt2, maxfct, mnum, mtype, phase, siznod_C, valAu, irowU, jcolU, perm2, nrhs, iparm2, msglvl, divJ, Uf, error )
        endif
        
        DO n=1, siznod_C
            m  = nods_C(n)  ! nodesL(m)  
            nn = m + nCells
            nl = m + 2*nCells
        
            k= ceiling(REAL(m)/( REAL(sdx*sdy) ) )
            IF (k == 1) THEN
                nij = m
            ELSE
                nij = m - (k-1)*sdx*sdy
            ENDIF
            j = ceiling( REAL(nij) / REAL(sdx) )  
            i = nij - (j - 1) * sdx 
        
            nim = geoPHYS_C(i-1,j,k)- 3*nCells;   nip = geoPHYS_C(i+1,j,k)- 3*nCells
            njm = geoPHYS_C(i,j-1,k)- 3*nCells;   njp = geoPHYS_C(i,j+1,k)- 3*nCells
            nkm = geoPHYS_C(i,j,k-1)- 3*nCells;   nkp = geoPHYS_C(i,j,k+1)- 3*nCells

            sx=0.5d0;sy=0.5d0;sz=0.5d0;
            if (nim <= 0) then
                nim = n;sx=1.d0
            endif
            if (nip <= 0) then
                nip = n;sx=1.d0
            endif
            if (njm <= 0) then
                njm = n;sy=1.d0
            endif
            if (njp <= 0) then
                njp = n;sy=1.d0
            endif
            if (nkm <= 0) then
                nkm = n;sz=1.d0
            endif
            if (nkp <= 0) then
                nkp = n;sz=1.d0
            endif
            ! Subtracting the gradient of a scalar field
            Jxsum(n) =  Jxsum(n) - sx*(Uf(nip) - Uf(nim) )/delta(1)
            Jysum(n) =  Jysum(n) - sy*(Uf(njp) - Uf(njm) )/delta(2)
            Jzsum(n) =  Jzsum(n) - sz*(Uf(nkp) - Uf(nkm) )/delta(3)
        enddo
        Jxsum(cel_bndUx- 3*nCells)=0.d0
        Jysum(cel_bndUy- 3*nCells)=0.d0
        Jzsum(cel_bndUz- 3*nCells)=0.d0
        
    enddo 
    !=====================================================================
    !== end cycle for calc divJsum
    !=====================================================================
    Jaf(cel_bndX) = 0.0d0 
    Jaf(cel_bndY) = 0.0d0 
    Jaf(cel_bndZ) = 0.0d0  

    Uaf(cel_bndx) = 0.0d0    
    Uaf(cel_bndy) = 0.0d0 
    Uaf(cel_bndz) = 0.0d0  
endif 

!===========================================================================================
IF (Ntime >= Nprint .and. Ntime /=0 ) THEN  ! output with step dtt skip 1st point
    Nprint = Ntime  + Nout
    Npoint= Npoint + 1                      ! point counter with step dtt
    ios=0

    call writeVtk_field (Npoint,  sdx, sdy, sdz, nCells, delta, Uaf, Jaf,Jxsum,Jysum,Jzsum, geoPHYS_C,  siznod_C, files)
    IF (numMech /= 0) THEN
        CALL writeVtk_src ( Npoint, numfun,  NcellsX, NcellsY, NcellsZ,  new_nodesX, new_nodesY,  new_nodesZ,  &
                            sdx, sdy, delta, files)  
    endif
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

integer colX(10), colY(10), colZ(10), colU(40), colU2(10), colU3(6)
real(8) valX(10), valY(10), valZ(10), valU(40), valU2(10), valU3(6)

real(8) s, a, b, c,  sz, sy, sx, sxy, sxz, syz, dsx, dsy, dsz

integer countU, nAx, nAy, nAz, Lx, Ly, Lz, Lfi,Lfi2, Lfi3, nFix, nFiy, nFiz,  kFi,  nc
INTEGER  kim, kip,  kjm, kkm,  kjp, kkp
INTEGER  kipjp, kipjm, kimjp, kimjm, kipkp, kipkm, kimkp, kimkm, kjpkp, kjpkm,  kjmkp, kjmkm
! 
type espm
    integer im
    DOUBLE PRECISION em
    type (espm),pointer ::prec
end type espm
type(espm),pointer :: sp_colX,next_X, sp_colY,next_Y, sp_colZ,next_Z, sp_colU, next_U, sp_colU2,next_U2, sp_colU3,next_U3

nullify(sp_colX,next_X,  sp_colY,next_Y,  sp_colZ,next_Z,  sp_colU,next_U, sp_colU2,next_U2, sp_colU3,next_U3) 

if (siznod_C /= 0) then
    num_bndX=0; num_bndY=0; num_bndZ=0
    cel_bndX = [integer:: ];  cel_bndY = [integer:: ];  cel_bndZ = [integer:: ] 

    num_bndUx=0; num_bndUy=0; num_bndUz=0
    cel_bndUx = [integer:: ];  cel_bndUy = [integer:: ];  cel_bndUz = [integer:: ]
endif

num_nzX=0;  !irow(        1) = 1;  
num_nzY=0;  !irow(  nCells+1) = 1;  
num_nzZ=0;  !irow(2*nCells+1) = 1;
num_nz =0; ! base
num_nzU=0;  ! base for U
num_nzU2=0; ! irow(siznod_C+1) = 1;     div
num_nzU3=0;  ! rhs

sz = 1.d0/(delta(3)*delta(3))
sy = 1.d0/(delta(2)*delta(2))
sx = 1.d0/(delta(1)*delta(1)) 
dsx = 0.5d0/delta(1) ! 1/2dx
dsy = 0.5d0/delta(2) ! 1/2dy
dsz = 0.5d0/delta(3) ! 1/2dz
s = 2.d0*(sx + sy + sz);

a = 2.0d0/(dt*delta(1))  
b = 2.0d0/(dt*delta(2))
c = 2.0d0/(dt*delta(3))
sxy = 1.d0/(delta(1)*delta(2)) 
sxz = 1.d0/(delta(1)*delta(3)) 
syz = 1.d0/(delta(2)*delta(3)) 

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
        
        valU = 0.0d0;   colU = 0;  Lfi=0; 
        valU2 = 0.0d0;  colU2 = 0; Lfi2=0; 
        valU3 = 0.0d0;  colU3 = 0; Lfi3=0;
        
        nFix=0; nFiy=0; nFiz=0  
        kFi=0;

        IF (geoPHYS_C(i,j,k) /=0) THEN
            kFi=1;                           
            countU = countU + 1 
        endif
        ! nn = i + sdx*(j-1) + kdz*(k-1) 
        ! (i+1,j+1,k) nn + 1 + sdx
        ! (i-1,j-1,k  nn - 1 - sdx 
        ! (i+1,j,k+1) nn + 1 + kdz
        ! (i-1,j,k-1  nn - 1 - kdz
        kim = nn-1;  kjm = nn-sdx;  kkm = nn-kdz;
        kip = nn+1;  kjp = nn+sdx;  kkp = nn+kdz;
        kimjm = nn - 1 - sdx; kimjp = nn - 1 + sdx; kimkm = nn - 1 - kdz;  kimkp = nn - 1 + kdz;
        kipjm = nn + 1 - sdx; kipjp = nn + 1 + sdx; kipkm = nn + 1 - kdz;  kipkp = nn + 1 + kdz; 
        kjmkm = nn - sdx - kdz;  kjmkp = nn - sdx + kdz;
        kjpkm = nn + sdx - kdz;  kjpkp = nn + sdx + kdz;

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
        !  valPHYS(n,3) =  Vex * mu0 * sigma;  valPHYS(n,2)=sigms*mu0;   Vex=alPHYS(n,3)/valPHYS(n,2)
        if (kFi /= 0) then ! add diffus component
            valX(7) = valX(7) + 2.d0*valPHYS(n,2)/dt    ! 2*sigms*mu0/dt
            valY = valX; valZ = valX
                                ! PHYS_C( )%valdom  = 2.d0*valPHYS(kp,2)/dt = 2*sigma*mu0/dt
            if     (geoPHYS_C(i+1,j,k) == 0) then
                Lx = Lx + 1; colX(Lx) = geoPHYS_C(i,j,k);   valX(Lx) = -3.d0*valPHYS(n,2)*dsx ! -3*sigms*mu0/2dx
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
        nc = geoPHYS_C(i,j,k) ! 3*nCells + local
        if (nc /= 0 ) then
            nim = geoPHYS_C(i-1,j,k);   nip = geoPHYS_C(i+1,j,k)
            njm = geoPHYS_C(i,j-1,k);   njp = geoPHYS_C(i,j+1,k)
            nkm = geoPHYS_C(i,j,k-1);   nkp = geoPHYS_C(i,j,k+1) 
            
            ! 8 corners
            if     ( nim == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i-1 j-1 k-1 
                Lfi2 = 4;
                colU2(1:Lfi2) =   [nip, njp,  nkp, nc ] - 3*nCells
                valU2(1:Lfi2) =   [-sx, -sy,  -sz,  s ]

                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ -a,      -b,             -c ]
                
                colU(1:4) = [nip,        njp,    nkp,      nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,  -2.d0*sz,  s]       
                   colU(5:7) = colU3(1:3)
                   valU(5:7) = valU3(1:3)
                   Lfi = 7

                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i+1 j-1 k-1 
                Lfi2 = 4;
                colU2(1:Lfi2) = [nim,   njp,   nkp,   nc ] - 3*nCells
                valU2(1:Lfi2) =  [-sx,  -sy,   -sz,    s ]
                
                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ a,      -b,             -c ]
                
                colU(1:4) = [nim,           njp,    nkp,    nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,   -2.d0*sz,  s]       
                   colU(5:7) = colU3(1:3) 
                   valU(5:7) = valU3(1:3) 
                   Lfi = 7

                nFix = 1; nFiy = 1; nFiz = 1; 
                
            elseif ( nim == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i-1 j+1 k-1 
                Lfi2 = 4;
                colU2(1:Lfi2) = [nip, njm, nkp,  nc ] - 3*nCells  
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  s ]
                
                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ -a,      b,             -c ]
                
                colU(1:4) = [nip,           njm,     nkp,         nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,   -2.d0*sz,  s]       
                    colU(5:7) = colU3(1:3) 
                    valU(5:7) = valU3(1:3) 
                    Lfi = 7

                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i+1 j+1 k-1 
                
                Lfi2 = 4;
                colU2(1:Lfi2) = [nim,  njm,  nkp,  nc ] - 3*nCells
                valU2(1:Lfi2) =   [-sx, -sy,  -sz,  s ]
                
                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ a,      b,             -c ]
                
                colU(1:4) = [nim,       njm,          nkp,   nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,   -2.d0*sz,  s]       
                    colU(5:7) = colU3(1:3) 
                    valU(5:7) = valU3(1:3) 
                    Lfi = 7

                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nim == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i-1 j-1 k+1 
                Lfi2 = 4;
                colU2(1:Lfi2) = [nip,   njp,  nkm,  nc ] - 3*nCells
                valU2(1:Lfi2) =   [-sx, -sy,  -sz,  s]

                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ -a,      -b,             c ]
                
                colU(1:4) = [nip,           njp,     nkm,    nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,   -2.d0*sz,  s]       
                    colU(5:7) = colU3(1:3) 
                    valU(5:7) = valU3(1:3) 
                    Lfi = 7

                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i+1 j-1 k+1 

                Lfi2 = 4;
                colU2(1:Lfi2) = [nim, njp,  nkm,  nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  s ]

                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ a,      -b,             c ]
                
                colU(1:4) = [nim,       njp,          nkm,   nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,   -2.d0*sz,  s]       
                    colU(5:7) = colU3(1:3) 
                    valU(5:7) = valU3(1:3) 
                    Lfi = 7

                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nim == 0 .and. njp == 0 .and. nkp == 0 ) then ! not i-1 j+1 k+1 
                Lfi2 = 4;
                colU2(1:Lfi2) = [nip, njm,   nkm, nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  s ]

                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ -a,      b,             c ]
                
                colU(1:4) = [nip,        njm,        nkm,    nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,   -2.d0*sz,  s]       
                    colU(5:7) = colU3(1:3) 
                    valU(5:7) = valU3(1:3) 
                    Lfi = 7

                nFix = 1; nFiy = 1; nFiz = 1; 
            elseif ( nip == 0 .and. njp == 0 .and. nkp == 0 ) then  ! not  i+1 j+1 k+1     

                Lfi2 = 4;
                colU2(1:Lfi2) = [nim, njm,  nkm,  nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  s ]

                Lfi3 = 3  ! for rhs
                colU3(1:3) = [ nn,   nCells+nn, 2*nCells + nn]
                valU3(1:3) = [ a,      b,             c ]
                
                colU(1:4) = [nim,       njm,          nkm,   nc]       
                valU(1:4) = [-2.d0*sx, -2.d0*sy,   -2.d0*sz,  s]       
                    colU(5:7) = colU3(1:3) 
                    valU(5:7) = valU3(1:3) 
                   Lfi = 7
                   
                nFix = 1; nFiy = 1; nFiz = 1;  
                !                        12 edges
                ! edges  along X
            elseif ( njp == 0  .and. nkm == 0 ) then  ! not  j+1  k-1  
                Lfi2 = 5;
                colU2(1:Lfi2) = [nip, nim,  njm,  nkp,  nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sx,  -sy,  -sz,    s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [ nCells+nn, 2*nCells + nn]
                valU3(1:2) = [    b,               -c ]
                
                colU(1:5) = [nip,      nim,     njm,          nkp,         nc]       
                valU(1:5) = [-sx,     -sx,       -2.d0*sy,   -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7

                nFiy = 1; nFiz = 1; 
            elseif ( njm == 0  .and. nkm == 0 ) then   ! not  j-1  k-1   
                Lfi2 = 5;
                colU2(1:Lfi2) = [nip, nim, njp,   nkp,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sx,  -sy,  -sz,    s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [ nCells+nn, 2*nCells + nn]
                valU3(1:2) = [    -b,               -c ]
                
                colU(1:5) = [nip,      nim,     njp,          nkp,         nc]       
                valU(1:5) = [-sx,     -sx,       -2.d0*sy,   -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7

                nFiy = 1; nFiz = 1; 
            elseif ( njp == 0  .and. nkp == 0 ) then   ! not  j+1  k+1
                Lfi2 = 5
                colU2(1:Lfi2) = [nip,    nim,   njm,    nkm,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx,    -sx,  -sy,     -sz,    s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [ nCells+nn, 2*nCells + nn]
                valU3(1:2) = [    b,               c ]
                
                colU(1:5) = [nip,      nim,     njm,          nkm,         nc]       
                valU(1:5) = [-sx,     -sx,       -2.d0*sy,   -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7

                nFiy = 1; nFiz = 1; 
            elseif ( njm == 0  .and. nkp == 0 ) then   ! not  j-1  k+1 
                Lfi2 = 5;
                colU2(1:Lfi2) = [nip,  nim,  njp,  nkm,  nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx,  -sx,  -sy,  -sz,   s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [ nCells+nn, 2*nCells + nn]
                valU3(1:2) = [    -b,               c ]
                
                colU(1:5) = [nip,      nim,     njp,       nkm,      nc]       
                valU(1:5) = [-sx,     -sx,     -2.d0*sy,  -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7
                nFiy = 1; nFiz = 1; 
                
                ! edges  along  Y
            elseif ( nip == 0  .and. nkm == 0 ) then  ! not  i+1  k-1  
                Lfi2 = 5;
                colU2(1:Lfi2) = [nim, njm,  njp,  nkp,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sy,  -sz,    s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     2*nCells + nn]
                valU3(1:2) = [ a,         -c       ]
                
                colU(1:5) = [nim,         njp,     njm,     nkp,      nc]       
                valU(1:5) = [-2.d0*sx,   -sy,     -sy,     -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                   Lfi = 7

                nFix = 1; nFiz = 1; 
            elseif ( nim == 0  .and. nkm == 0 ) then   ! not  i-1  k-1  
                Lfi2 = 5;
                colU2(1:Lfi2) =   [ nip, njm,  njp,  nkp, nc ] - 3*nCells
                valU2(1:Lfi2) =   [-sx,  -sy,  -sy,  -sz,  s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     2*nCells + nn]
                valU3(1:2) = [ -a,      -c ]
                
                colU(1:5) = [nip,       njp,    njm,     nkp,      nc]       
                valU(1:5) = [-2.d0*sx,  -sy,    -sy,     -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                   Lfi = 7
                   
                nFix = 1; nFiz = 1; 
            elseif ( nip == 0  .and. nkp == 0 ) then  ! not  i+1  k+1  
                Lfi2 = 5;
                colU2(1:Lfi2) = [nim, njm, njp, nkm, nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,-sy, -sz,  s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     2*nCells + nn]
                valU3(1:2) = [ a,               c          ]
                
                colU(1:5) = [nim,        njp,   njm,     nkm,      nc]       
                valU(1:5) = [-2.d0*sx,  -sy,    -sy,     -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7

                nFix = 1;  nFiz = 1; 
            elseif ( nim == 0  .and. nkp == 0 ) then   ! not  i-1  k+1   
                Lfi2 = 5;
                colU2(1:Lfi2) = [nip, njm,  njp,  nkm,  nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy, -sy,   -sz,  s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     2*nCells + nn]
                valU3(1:2) = [ -a,          c      ]
                
                colU(1:5) = [nip,        njp,   njm,     nkm,      nc]       
                valU(1:5) = [-2.d0*sx,  -sy,    -sy,     -2.d0*sz,   s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                   Lfi = 7

                nFix = 1;  nFiz = 1; 
                ! edges  along Z 
            elseif ( nim == 0  .and. njm == 0 ) then   ! not  i-1 j-1 
                Lfi2 = 5;
                colU2(1:Lfi2) = [nip,  njp, nkp,  nkm,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  -sz,    s  ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     nCells + nn]
                valU3(1:2) = [ -a,               -b          ]
                
                colU(1:5) = [nip,         njp,          nkp,     nkm,   nc]       
                valU(1:5) = [-2.d0*sx,  -2.d0*sy,   -sz,     -sz,     s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7

                nFix = 1; nFiy = 1; 
            elseif ( nip == 0 .and. njm == 0 ) then    ! not  i+1 j-1 ! 
                Lfi2 = 5;
                colU2(1:Lfi2) = [nim,  njp, nkp,  nkm,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  -sz,    s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     nCells + nn]
                valU3(1:2) = [ a,         -b          ]
                
                colU(1:5) = [nim,         njp,     nkp,     nkm,   nc]       
                valU(1:5) = [-2.d0*sx,  -2.d0*sy,   -sz,     -sz,     s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7
                nFix = 1; nFiy = 1;  
            elseif ( nim == 0  .and. njp == 0 ) then   ! not  i-1 j+1 
                Lfi2 = 5;
                colU2(1:Lfi2) = [nip, njm,  nkp,  nkm,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  -sz,    s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     nCells + nn]
                valU3(1:2) = [ -a,         b     ]
                
                colU(1:5) = [nip,         njm,     nkp,     nkm,   nc]       
                valU(1:5) = [-2.d0*sx,  -2.d0*sy,   -sz,     -sz,     s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7

                nFix = 1; nFiy = 1; 
            elseif ( nip == 0  .and. njp == 0 ) then   ! not   i+1 j+1  
                Lfi2 = 5;
                colU2(1:Lfi2) = [nim, njm,   nkm, nkp,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sz,  -sz,    s ]

                Lfi3 = 2  ! for rhs
                colU3(1:2) = [nn,     nCells + nn]
                valU3(1:2) = [ a,               b ]
                
                colU(1:5) = [nim,         njm,     nkp,     nkm,   nc]       
                valU(1:5) = [-2.d0*sx,  -2.d0*sy,   -sz,     -sz,     s  ]       
                    colU(6:7) = colU3(1:2) 
                    valU(6:7) = valU3(1:2) 
                    Lfi = 7

                nFix = 1; nFiy = 1;  
                
                ! 6 faces  
            elseif (  nip == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                 ! not   i+1
                Lfi2 = 6;
                colU2(1:Lfi2) = [nim, njm,  njp,  nkm,  nkp,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sy,  -sz,  -sz,    s ]

                Lfi3 = 1  ! for rhs
                colU3(1:1) = [  nn  ]
                valU3(1:1) = [   a   ]

                colU(1:6) = [nim,        njp,    njm,    nkp,     nkm,    nc]       
                valU(1:6) = [-2.d0*sx,   -sy,   -sy,      -sz,      -sz,   s  ]       
                    colU(7) = colU3(1) 
                    valU(7) = valU3(1) 
                    Lfi = 7

                nFix = 1;
            elseif (  nim == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                ! not   i-1 x-
                Lfi2 = 6;
                colU2(1:Lfi2) = [nip, njm,  njp,  nkm,  nkp,   nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sy,  -sy,  -sz,  -sz,    s ]

                Lfi3 = 1  ! for rhs
                colU3(1:1) = [  nn  ]
                valU3(1:1) = [   -a   ]

                colU(1:6) = [nip,        njp,    njm,    nkp,     nkm,    nc]       
                valU(1:6) = [-2.d0*sx,   -sy,   -sy,      -sz,      -sz,        s  ]       
                    colU(7) = colU3(1) 
                    valU(7) = valU3(1) 
                    Lfi = 7

                nFix = 1; 
            elseif (  njp == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                 ! not   j+1 
                Lfi2 = 6;
                colU2(1:Lfi2) = [nim, nip,   njm,  nkm,  nkp,  nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx, -sx,  -sy,  -sz,  -sz,  s ]

                Lfi3 = 1  ! for rhs
                colU3(1:1) = [  nCells + nn]
                valU3(1:1) = [          b  ]

                colU(1:6) = [nip,   nim,    njm,        nkp,     nkm,    nc]       
                valU(1:6) = [-sx,    -sx,   -2.d0*sy,   -sz,      -sz,        s  ]       
                    colU(7) = colU3(1) 
                    valU(7) = valU3(1) 
                   Lfi = 7

                nFiy = 1; 
            elseif (  njm == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
                ! not   j-1
                Lfi2 = 6;
                colU2(1:Lfi2) = [nim,   nip,    njp,     nkm,   nkp,  nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx,  -sx,    -sy,     -sz,  -sz,   s ]

                Lfi3 = 1  ! for rhs
                colU3(1:1) = [  nCells + nn]
                valU3(1:1) = [          -b ]

                colU(1:6) = [nip,   nim,    njp,        nkp,     nkm,    nc]       
                valU(1:6) = [-sx,    -sx,   -2.d0*sy,   -sz,      -sz,        s  ]       
                    colU(7) = colU3(1) 
                    valU(7) = valU3(1) 
                    Lfi = 7

                nFiy = 1;
            elseif (  nkp == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
                ! not   k+1 
                Lfi2 = 6;
                colU2(1:Lfi2) = [nim,    nip,  njm,  njp,    nkm,    nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx,    -sx,  -sy,  -sy,    -sz,     s ]

                Lfi3 = 1  ! for rhs
                colU3(1:1) = [  2*nCells + nn]
                valU3(1:1) = [          c    ]

                colU(1:6) = [nip,   nim,   njp,  njm,     nkm,        nc]       
                valU(1:6) = [-sx,   -sx,   -sy,  -sy,    -2.d0*sz,     s  ]       
                    colU(7) = colU3(1) 
                    valU(7) = valU3(1) 
                    Lfi = 7

                nFiz = 1; 
            elseif (  nkm == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
                ! not   k-1
                Lfi2 = 6;
                colU2(1:Lfi2) = [nim,   nip,   njm,  njp,   nkp,    nc ] - 3*nCells
                valU2(1:Lfi2) = [-sx,   -sx,  -sy,   -sy,  -sz,      s]

                Lfi3 = 1  ! for rhs
                colU3(1:1) = [  2*nCells + nn]
                valU3(1:1) = [          -c   ]

                colU(1:6) = [nip,   nim,   njp,  njm,     nkp,        nc]       
                valU(1:6) = [-sx,    -sx,   -sy,   -sy,    -2.d0*sz,     s  ]       
                    colU(7) = colU3(1) 
                    valU(7) = valU3(1) 
                    Lfi = 7
                nFiz = 1; 
            else
                Lfi2 = 7;  ! for calc div
                colU2(1:Lfi2) = [nim,   nip,   njm,  njp, nkm,  nkp,    nc ] - 3*nCells
                valU2(1:Lfi2) =   [-sx, -sx,  -sy,  -sy,  -sz,  -sz,     s]
                
                Lfi3 = 6   ! for rhs 
                colU3(1:Lfi3) =        [ kip,  kim,  nCells+kjp,  nCells+kjm,  2*nCells+kkp,  2*nCells+kkm  ]
                valU3(1:Lfi3) = 0.25d0*[-a,    a,     -b,                   b,      -c,              c]
                
                colU(1:7) = [nim,  nip, njm, njp, nkm, nkp, nc]
                valU(1:7) = [-sx, -sx, -sy, -sy, -sz, -sz,  s ]
                colU(8:13) = colU3(1:Lfi3)
                valU(8:13) = valU3(1:Lfi3)
                Lfi = 13
            endif
            
            k0=0 ! for base
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
            
            k0=0 ! for div
            mm2: do k1=1,Lfi2-1  
                do k2=k1+1, Lfi2
                    if  (colU2(k1) == colU2(k2) ) then 
                        k0 = colU2(k2)
                        exit mm2
                    endif
                enddo
            enddo mm2
            if (k0 /=0 ) then
                print*,'node Fi2 double', k0, 'i=',i, 'j=',j, 'k=',k
                stop
            endif
            
            k0=0  ! for rhs
            mm3: do k1=1,Lfi3-1  
                do k2=k1+1, Lfi3
                    if  (colU3(k1) == colU3(k2) ) then 
                        k0 = colU3(k2)
                        exit mm3
                    endif
                enddo
            enddo mm3
            if (k0 /=0 ) then
                print*,'node Fi3 double', k0, 'i=',i, 'j=',j, 'k=',k
                stop
            endif
            
            if ( nFix == 1 ) cel_bndUx = [cel_bndUx, nc]
            if ( nFiy == 1 ) cel_bndUy = [cel_bndUy, nc]
            if ( nFiz == 1 ) cel_bndUz = [cel_bndUz, nc]
            !
            call full_sort(colU, valU, Lfi, 1,1)
            call full_sort(colU2, valU2, Lfi2, 1,1)
            call full_sort(colU3, valU3, Lfi3, 1,1)
            
            do m = 1, Lfi      ! for base U
                if (colU(m) <= 0) then
                    print*, 'cell=',nn, 'm',m,'colU(m)',colU(m)
                    stop
                endif
                allocate(sp_colU)
                sp_colU=espm(colU(m),valU(m),next_U) 
                num_nzU = num_nzU + 1 
                next_U => sp_colU
            enddo
            
            do m = 1, Lfi2 ! for div
                if (colU2(m) <= 0) then
                    print*, 'cell=',nn, 'm',m,'colU2(m)',colU2(m)
                    stop
                endif
                allocate(sp_colU2)
                sp_colU2=espm(colU2(m),valU2(m),next_U2) 
                num_nzU2 = num_nzU2 + 1 
                next_U2 => sp_colU2
            enddo
            
            do m = 1, Lfi3 ! for rhs
                if (colU3(m) <= 0) then
                    print*, 'cell=',nn, 'm',m,'colU3(m)',colU2(m)
                    stop
                endif
                allocate(sp_colU3)
                sp_colU3=espm(colU3(m),valU3(m),next_U3) 
                num_nzU3 = num_nzU3 + 1 
                next_U3 => sp_colU3
            enddo

            irow( 3*nCells + countU + 1) = irow( 3*nCells + countU) + Lfi ! base U
            irowU(  countU + 1) = irowU( countU) + Lfi2  ! for div
            irowUrhs(  countU + 1) = irowUrhs( countU) + Lfi3 ! for rhs
        endif !geoPHYS_C /=0
        
    endif ! kFi/=0  
            
    enddo !i
  enddo ! j
enddo ! k

! print*,'end form Y'

num_nz = num_nzX + num_nzY + num_nzZ + num_nzU
num_bndX = size(cel_bndX );  num_bndY = size(cel_bndY );  num_bndZ = size(cel_bndZ );  

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

allocate(jcol(num_nz), jcolU(num_nzU2),  jcolUrhs(num_nzU3),source=0)
allocate(valA(num_nz), valAu(num_nzU2),  valAUrhs(num_nzU3),source=0.0d0)

if (size(jcol) < irow( 3*nCells + countU + 1)-1 ) then
    print*,'base ', i,j,k, countU, ' irow=', irow( 3*nCells + countU + 1) , ' size jcol=', size(jcol)
    stop
endif

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

! PRINT '(a,i3,a,i3,a,i3,     a,i9, a,g12.5,a/)', 'sparse matrix generation completed on grid (', sdx,' x ',sdy,' x ',sdz, &
                 ! ' ), Non zero elem= ',num_nz, ' Density of matrix:', 100.0* REAL(num_nz)/REAL(nCells)/REAL(nCells),'%'
i= num_nzU2
DO WHILE (associated(sp_colU2))
    jcolU(i) = sp_colU2%im  
    valAu(i) = sp_colU2%em   !  valAu
    i=i-1
    sp_colU2 => sp_colU2%prec  
END DO
nullify(sp_colU2, next_U2)

! PRINT '(a,    a,i9, a,g12.5,a/)', 'sparse matrix U2 generation completed ',  &
                 ! ' Non zero elem= ',num_nzU2, ' Density of matrix:', 100.0* REAL(num_nzU2)/REAL(siznod_C)/REAL(siznod_C),'%'
i= num_nzU3
DO WHILE (associated(sp_colU3))
    jcolUrhs(i) = sp_colU3%im  
    valAUrhs(i) = sp_colU3%em   
    i=i-1
    sp_colU3 => sp_colU3%prec  
END DO
nullify(sp_colU3, next_U3)

! PRINT '(a,    a,i9, a,g12.5,a/)', 'sparse matrix U3 rhs generation completed ',  &
                 ! ' Non zero elem= ',num_nzU3, ' Density of matrix:', 100.0* REAL(num_nzU3)/REAL(siznod_C)/REAL(siznod_C),'%'
PRINT '(a)', 'sparse matrixs generation completed!'

END SUBROUTINE gen_sparse_matrix
!==============================================
!================================================
SUBROUTINE update_rhs
    ! 8 corners
    if     ( nim == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i-1 j-1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy+Vez*sxz,  Vex*sxy+Vey*sy+Vez*syz,  Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjp2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njm == 0 .and. nkm == 0  ) then ! not i+1 j-1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy-Vez*sxz,  -Vex*sxy+Vey*sy+Vez*syz,  -Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjp2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nim == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i-1 j+1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy+Vez*sxz,  -Vex*sxy+Vey*sy-Vez*syz,  Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njp == 0 .and. nkm == 0 ) then ! not i+1 j+1 k-1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy-Vez*sxz,  Vex*sxy+Vey*sy-Vez*syz,  -Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkp,       kkp2, nCells+kkp, nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nim == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i-1 j-1 k+1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy-Vez*sxz,  Vex*sxy+Vey*sy-Vez*syz,  -Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjp2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njm == 0 .and. nkp == 0 ) then ! not i+1 j-1 k+1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy+Vez*sxz,  -Vex*sxy+Vey*sy-Vez*syz,  Vex*sxz-Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =     [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =     [kjp,       kjm2, nCells+kjp, nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        4d0*syz,      -syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =     [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   4d0*syz,    -syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nim == 0 .and. njp == 0 .and. nkp == 0 ) then ! not i-1 j+1 k+1 
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx-Vey*sxy-Vez*sxz,  -Vex*sxy+Vey*sy+Vez*syz,  -Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =      [kip,       kip2, nCells+kip, nCells+kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   4d0*sxy,    -sxy,        4d0*sxz,      -sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =      [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ 4d0*sxy,   -sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =      [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ 4d0*sxz,   -sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
    elseif ( nip == 0 .and. njp == 0 .and. nkp == 0 ) then  ! not  i+1 j+1 k+1     
        if (Nx /=0 .or. Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+3) =      [    L,                      L+nCells,                     L+2*nCells ]
            valU(Lfi+1:Lfi+3) = -3d0*[ Vex*sx+Vey*sxy+Vez*sxz,  Vex*sxy+Vey*sy+Vez*syz,  Vex*sxz+Vey*syz+Vez*sz] 
            Lfi = Lfi + 3
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+6) =      [kim,       kim2, nCells+kim, nCells+kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+6) = -Vex*[ -4d0*sx,   sx,   -4d0*sxy,    sxy,        -4d0*sxz,      sxz    ]
            Lfi = Lfi + 6
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+6) =      [kjm,       kjm2, nCells+kjm, nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+6) = -Vey*[ -4d0*sxy,   sxy,   -4d0*sy,    sy,        -4d0*syz,      syz    ]
            Lfi = Lfi + 6
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+6) =      [kkm,       kkm2, nCells+kkm, nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2 ] 
            valU(Lfi+1:Lfi+6) = -Vez*[ -4d0*sxz,   sxz,   -4d0*syz,    syz,         -4d0*sz,      sz    ]
            Lfi = Lfi + 6
        endif
        !                        12 edges
        ! edges  along X
    elseif ( njp == 0  .and. nkm == 0 ) then  ! not  j+1  k-1  
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy - Vez*sxz,  -Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ sxy,        -sxy,        -sxz,          sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [nCells+kjm,  nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         4d0*syz,     -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kkp,   nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*syz,      -syz,        -4d0*sz,       sz ]
            Lfi = Lfi + 4
        endif
    elseif ( njm == 0  .and. nkm == 0 ) then   ! not  j-1  k-1   
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy + Vez*sxz,  Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -sxy,        sxy,        -sxz,          sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [nCells+kjp,  nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         -4d0*syz,     syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kkp,   nCells+kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*syz,      syz,        -4d0*sz,       sz ]
            Lfi = Lfi + 4
        endif
    elseif ( njp == 0  .and. nkp == 0 ) then   ! not  j+1  k+1
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy + Vez*sxz,  Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ sxy,        -sxy,        sxz,          -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [nCells+kjm,  nCells+kjm2, 2*nCells+kjm, 2*nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         -4d0*syz,      syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [nCells+kkm,   nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*syz,      syz,        -4d0*sz,        sz ]
            Lfi = Lfi + 4
        endif
    elseif ( njm == 0  .and. nkp == 0 ) then   ! not  j-1  k+1 
        if (Ny /= 0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L+nCells,    L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vey*sy - Vez*sxz,  -Vey*syz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =      [nCells+kip, nCells+kim, 2*nCells+kip, 2*nCells+kim ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -sxy,        sxy,        sxz,          -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =      [nCells+kjp,  nCells+kjp2, 2*nCells+kjp, 2*nCells+kjp2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sy,       sy,         4d0*syz,     -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =      [nCells+kkm,   nCells+kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*syz,      -syz,        -4d0*sz,       sz ]
            Lfi = Lfi + 4
        endif
        ! edges  along  Y
    elseif ( nip == 0  .and. nkm == 0 ) then  ! not  i+1  k-1  
        if (Nx /=0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L,             L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx - Vez*sxz,  -Vex*sxz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kim,     kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxz,      -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ sxy,  -sxy,   -syz,         syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,       kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*sxz,  -sxz,    -4d0*sz,            sz ]
            Lfi = Lfi + 4
        endif
    elseif ( nim == 0  .and. nkm == 0 ) then   ! not  i-1  k-1  
        if (Nx /=0 .or. Nz /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L,             L+2*nCells ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx + Vez*sxz,  +Vex*sxz + Vez*sz] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxz,      sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -sxy,  sxy,   -syz,         syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,       kkp2, 2*nCells+kkp, 2*nCells+kkp2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*sxz,  sxz,    -4d0*sz,        sz ]
            Lfi = Lfi + 4
        endif
    elseif ( nip == 0  .and. nkp == 0 ) then  ! not  i+1  k+1  
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx + Vey*sxy,  +Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kim,     kim2, 2*nCells+kim, 2*nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxy,    sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ sxy,  -sxy,   syz,         -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkm,     kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -4d0*sxz,   sxz,    -4d0*sz,            sz ]
            Lfi = Lfi + 4
        endif
    elseif ( nim == 0  .and. nkp == 0 ) then   ! not  i-1  k+1   
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx - Vey*sxy,  -Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, 2*nCells+kip, 2*nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxz,      -sxz ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjp,  kjm, 2*nCells+kjp, 2*nCells+kjm ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -sxy,  sxy,   syz,         -syz ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkm,       kkm2, 2*nCells+kkm, 2*nCells+kkm2  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ 4d0*sxz,  -sxz,    -4d0*sz,        sz ]
            Lfi = Lfi + 4
        endif
        ! edges  along Z 
    elseif ( nim == 0  .and. njm == 0 ) then   ! not  i-1 j-1 
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx + Vey*sxy,  +Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, nCells+kip, nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxy,    sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjm,     kjm2, nCells+kjm, nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sxy,  sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,     kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -sxz,   sxz,    -syz,       syz ]
            Lfi = Lfi + 4
        endif
    elseif ( nip == 0 .and. njm == 0 ) then    ! not  i+1 j-1 ! 
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[ Vex*sx - Vey*sxy,  -Vex*sxy - Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =      [ kim,     kim2, nCells+kim, nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxy,    -sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =      [ kjp,     kjp2, nCells+kjp, nCells+kjp2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ 4d0*sxy, -sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =      [ kkp,   kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ sxz, -sxz,    -syz,       syz ]
            Lfi = Lfi + 4
        endif
    elseif ( nim == 0  .and. njp == 0 ) then   ! not  i-1 j+1 
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =     [        L,             L+nCells    ]
            valU(Lfi+1:Lfi+2) = -3d0*[Vex*sx - Vey*sxy,  -Vex*sxy + Vey*sy] 
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =     [ kip,     kip2, nCells+kip, nCells+kip2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   4d0*sxy,    -sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =     [ kjm,     kjm2, nCells+kjm, nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[  4d0*sxy, -sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =     [ kkp,   kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ -sxz, sxz,    syz,       -syz ]
            Lfi = Lfi + 4
        endif
    elseif ( nip == 0  .and. njp == 0 ) then   ! not   i+1 j+1  
        if (Nx /=0 .or. Ny /= 0 ) then
            colU(Lfi+1:Lfi+2) =      [        L,              L+nCells   ]
            valU(Lfi+1:Lfi+2) = -3d0*[ (Vex*sx+Vey*sxy),  (Vex*sxy+Vey*sy) ]
            Lfi = Lfi + 2
        endif
        if (Nx /=0) then
            colU(Lfi+1:Lfi+4) =      [ kim,     kim2, nCells+kim, nCells+kim2 ] 
            valU(Lfi+1:Lfi+4) = -Vex*[ -4d0*sx,  sx,   -4d0*sxy,    sxy ]
            Lfi = Lfi + 4
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =      [ kjm,     kjm2, nCells+kjm, nCells+kjm2 ] 
            valU(Lfi+1:Lfi+4) = -Vey*[ -4d0*sxy,  sxy,   -4d0*sy,    sy ]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =      [ kkp,   kkm, nCells+kkp, nCells+kkm  ]
            valU(Lfi+1:Lfi+4) = -Vez*[ sxz, -sxz,    syz,       -syz ]
            Lfi = Lfi + 4
        endif
        ! 6 faces  
    elseif (  nip == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
         ! not   i+1
        if (Nx /=0) then
            colU(Lfi+1:Lfi+3) =          [ L,  kim,  kim2 ] 
            valU(Lfi+1:Lfi+3) = Vex*sx * [ -3d0, 4d0, -1d0 ]
            Lfi = Lfi + 3
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =           [  kjp,  kjm]
            valU(Lfi+1:Lfi+2) = Vey*sxy * [ -1d0,  1d0]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =           [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) = Vez*sxz * [ -1d0, 1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  nim == 0 .and. njp /= 0 .and. njm /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
        ! not   i-1 x-
        if (Nx /=0) then
            colU(Lfi+1:Lfi+3) =          [ L,  kip,  kip2 ] 
            valU(Lfi+1:Lfi+3) = Vex*sx * [ -3d0, 4d0, -1d0 ]
            Lfi = Lfi + 3
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =           [  kjp,  kjm]
            valU(Lfi+1:Lfi+2) = Vey*sxy * [ 1d0,  -1d0]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =          [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) =Vez*sxz * [ 1d0,  -1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  njp == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
         ! not   j+1 
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) = Vex*sxy * [ -1d0, 1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+3) = nCells + [ L,  kjm,  kjm2 ]
            valU(Lfi+1:Lfi+3) = Vey*sy * [ -3d0, 4d0,  -1d0  ]
            Lfi = Lfi + 3
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) = Vez*syz * [ -1d0,  1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  njm == 0 .and. nip /= 0 .and. nim /= 0 .and. nkp /= 0 .and. nkm /= 0  ) then
        ! not   j-1
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) = Vex*sxy * [ 1d0, -1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+3) = nCells + [ L,  kjp,  kjp2 ]
            valU(Lfi+1:Lfi+3) = Vey*sy * [ -3d0, 4d0,  -1d0  ]
            Lfi = Lfi + 3
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+2) =  nCells + [ kkp,  kkm ]
            valU(Lfi+1:Lfi+2) = Vez*syz * [ 1d0,  -1d0 ]
            Lfi = Lfi + 2
        endif
    elseif (  nkp == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
        ! not   k+1 
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) = 2*nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) =  Vex*sxz * [ -1d0, 1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =2*nCells + [  kjp,  kjm ]
            valU(Lfi+1:Lfi+2) = Vey*syz * [ -1d0,  1d0 ]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+3) =2*nCells + [L,  kkm,  kkm2 ]
            valU(Lfi+1:Lfi+3) =  Vez*sz * [-3d0, 4d0,  -1d0 ]
            Lfi = Lfi + 3
        endif
            
     elseif (  nkm == 0 .and. nip/= 0 .and. nim /= 0 .and. njp /= 0 .and. njm /= 0  ) then
        ! not   k-1
        if (Nx /=0) then
            colU(Lfi+1:Lfi+2) = 2*nCells + [ kip,  kim ] 
            valU(Lfi+1:Lfi+2) =  Vex*sxz * [ 1d0, -1d0 ]
            Lfi = Lfi + 2
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+2) =2*nCells + [  kjp,  kjm ]
            valU(Lfi+1:Lfi+2) = Vey*syz * [ 1d0,  -1d0 ]
            Lfi = Lfi + 2
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+3) =2*nCells + [ L,  kkp,  kkp2 ]
            valU(Lfi+1:Lfi+3) =  Vez*sz * [ -3d0, 4d0,  -1d0 ]
            Lfi = Lfi + 3
        endif
    else 
        Lfi = 0
        if (Nx /=0) then
            colU(Lfi+1:Lfi+3)=         [  L,  kip, kim ]
            valU(Lfi+1:Lfi+3)= Vex*sx* [ -2d0, 1d0,  1d0] 
            Lfi= Lfi + 3
            
            colU(Lfi+1:Lfi+4) =         nCells + [ kipjp, kipjm, kimjp, kimjm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vex*sxy * [ 1d0,   -1d0,   -1d0,   1d0 ]
            Lfi = Lfi + 4 
            
            colU(Lfi+1:Lfi+4) =       2*nCells + [kipkp, kipkm, kimkp, kimkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vex*sxz * [ 1d0,   -1d0,   -1d0,   1d0]
            Lfi = Lfi + 4     
        endif
        if (Ny /= 0) then  
            colU(Lfi+1:Lfi+4) =                  [ kipjp, kipjm, kimjp, kimjm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vey*sxy * [ 1d0,   -1d0,   -1d0,   +1d0 ]
            Lfi = Lfi + 4
            
            colU(Lfi+1:Lfi+3)= nCells + [ L,    kjp, kjm ]
            valU(Lfi+1:Lfi+3)=  Vey*sy* [ -2.d0, 1d0,  1d0] 
            Lfi= Lfi + 3

            colU(Lfi+1:Lfi+4) =       2*nCells + [kjpkp, kjpkm, kjmkp, kjmkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vey*syz * [ 1d0,   -1d0,   -1d0,   1d0]
            Lfi = Lfi + 4
        endif
        if (Nz /= 0) then
            colU(Lfi+1:Lfi+4) =                  [ kipkp, kimkp, kipkm, kimkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vez*sxz * [ 1d0,   -1d0,   -1d0,   1d0 ]
            Lfi = Lfi + 4
            
            colU(Lfi+1:Lfi+4) =         nCells + [kjpkp, kjpkm, kjmkp, kjmkm ]
            valU(Lfi+1:Lfi+4) = 0.25d0*Vez*syz * [ 1d0,   -1d0,   -1d0,   1d0]
            Lfi = Lfi + 4
            
            colU(Lfi+1:Lfi+3)= 2*nCells + [ L,    kjp, kjm ]
            valU(Lfi+1:Lfi+3)=   Vez*sz * [ -2.d0, 1d0,  1d0] 
            Lfi= Lfi + 3
        endif
    endif

END SUBROUTINE update_rhs

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


