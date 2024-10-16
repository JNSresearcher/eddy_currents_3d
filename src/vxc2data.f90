INCLUDE 'm_fparser.f90'

! This program converts "VoxCad" format for calculating 3D fields.

! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )

SUBROUTINE vxc2data ( delta,    dt,        Time,     dtt,                   &
                      sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, size_PHYS_C, &
                      numfun, numMech,     &
                      solv,   files,    tolerance, itmax,    bound,  BND ) !
USE m_vxc2data
USE m_fparser
IMPLICIT NONE

REAL(8), INTENT(OUT) :: delta(3),  &         ! Grid spacing along X, Y and Z
                        dt,        &         ! time step
                        Time,      &         ! stop time
                        dtt                  ! jump duration
INTEGER, INTENT(OUT) :: sdx,  sdy,  sdz,  &  ! Number of cells along X,  Y and  Z
                        nsub,             &  ! number of physical domains
                        nsub_air,         &  ! number of environment domains
                        nsub_glob            ! total number of domains
INTEGER, INTENT(OUT) :: numfun,       &      ! number of functions for calculating external sources
                        numMech,      &      ! number of functions for calculating movements of external sources
                        size_PHYS_C

CHARACTER(LEN=3), INTENT(OUT) :: solv        ! character string corresponding to the name of methods: 'BCG' or 'SOR'
CHARACTER(16),    INTENT(OUT) :: files       ! name for output files
REAL(8),          INTENT(OUT) :: tolerance   ! convergence criterion
INTEGER,          INTENT(OUT) :: itmax       ! maximum number of iterations

! character string for boundary conditions on 6 faces: 
CHARACTER(6),     INTENT(OUT) :: bound
REAL(8),          INTENT(OUT) :: BND(3,2)   ! values of boundary conditions on 6 faces
!-------------------!
! working variables !
!-------------------!
REAL(8) delta0
CHARACTER (LEN=6), ALLOCATABLE ::  nameFun(:), nameVmech(:)
INTEGER(1),        ALLOCATABLE::   v(:), v_e(:)

TYPE tmp
    INTEGER jm                   ! function number
    INTEGER im                   ! cell number
    CHARACTER (len = 1)  ch      ! function type: X Y Z or 0
    TYPE (tmp),pointer ::prec
END TYPE tmp
TYPE(tmp),pointer :: sp_F,next_F 

CHARACTER (len = 10) ch_e
CHARACTER (len = 50) words(20) !
CHARACTER (len = 74) letter
CHARACTER (len=:),ALLOCATABLE :: st, ch

INTEGER m1, m2, kdz, Cells, nn, i,j,k,  L, m, n, k0, k1, k2, ios
INTEGER num_nodx, num_nody, num_nodz, num_nodV! sign of the presence of a vector field
INTEGER lst, neww, np, kp, ier, idev
LOGICAL uncompress
CHARACTER (LEN=4)::  varconst(12)
REAL(8)          ::  valconst(12)

! namelist /WRITEdata/  &
    ! delta,    dt,        Time,     dtt,                   &
    ! sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, size_PHYS_C, &
    ! numfun, numMech, num_nodx, num_nody,  num_nodz, num_nodV,     &
    ! solv,   files,    tolerance, itmax,    bound 

REAL(8)  numeric
EXTERNAL numeric

letter='123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz'

! init data DEFAULT for solver
solv='BCG';    tolerance=1.d-3; itmax=10000;    bound='DDDDDD'; BND = -0.95d0;  files='out';  

OPEN(1,file='in.vxc')
! input structure
k0=0;k1=0;k2=0; kp=0;  

numfun = 0;  numMech=0; size_PHYS_C=0;

m=0 !  m - it is counter CDATA - not USE in DO WHILE
convert_strings0: &
DO WHILE ( 1==1)
    idev = 1
    CALL readline(st,ch,ier,idev)             !============================= 1 idev= 1
    
    ! if(ier.ne.0) exit convert_strings0
    IF (ier == -1) EXIT convert_strings0
    
    IF (st(1:6)=='</VXC>' ) EXIT
    lst = len_trim(st)
    
    m1=index(st,'<Lattice_Dim>')
    IF (m1>=1) THEN
        m2=index(st,'</Lattice_Dim>')
        ch_e=trim(st(m1+13:m2-1))
        delta0=numeric(ch_e)
        PRINT '( a,g10.3,$ )', 'delta0=',delta0
    ENDIF
    m1=index(st,'<X_Dim_Adj>')
    IF (m1>=1) THEN
        m2=index(st,'</X_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(1)=numeric(ch_e) * delta0
        WRITE(*,'( a,g10.3,$ )') ' deltaX=',REAL(delta(1))
    ENDIF
    m1=index(st,'<Y_Dim_Adj>')
    IF (m1>=1) THEN
        m2=index(st,'</Y_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(2)=numeric(ch_e) * delta0
        WRITE(*,'( a,g10.3,$ )') ' deltaY=',REAL(delta(2))
    ENDIF
    m1=index(st,'<Z_Dim_Adj>')
    IF (m1>=1) THEN
        m2=index(st,'</Z_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(3)=numeric(ch_e) * delta0
        WRITE(*,'( a,g10.3 )') ' deltaZ=',REAL(delta(3))
    ENDIF
    m1=index(st,'Material ID')
    IF (m1>=1) THEN
        kp=kp+1
    ENDIF

    k0=index(st,'<Name>')
    IF (k0>0 ) THEN ! 
        n=index(st,'</')
        st=st(k0+6:n-1)
        i=1
        DO
            k = index( st(i:),'=')
            IF (k/=0) THEN
                i=i+k-1; st(i:i) = ' '
            ELSE
                EXIT
            ENDIF
        ENDDO
        
        CALL Upp (st, st)  
        CALL string2words(trim(st), words, neww)
        
        IF (neww > 20 ) THEN
            PRINT*,'error: numbers words >20  newws=: ', neww,' correct dim words(20) string 64 in vxc2data.f90 '
            STOP
        ENDIF
        !  calc numbers L-domain
        DO i=2,neww
            IF (index( words(i)(1:1),'D') /=0 ) THEN
                IF (i+2 >= neww) CYCLE ! it is last word
                
                DO j=i+2, neww
                    IF ( index( words(j)(1:1),'C')/=0) THEN ! with environ
                        size_PHYS_C = size_PHYS_C + 1
                        ! PRINT*,'size_PHYS_C=',size_PHYS_C
                    ENDIF
                ENDDO

                IF ( index( words(i+2),'SRC') /=0) THEN 
                    DO j=i+2, neww
                        IF     ( index( words(j),'SRCY') /=0 ) THEN ! what is SRC?
                            numfun = numfun + 1
                        ELSEIF ( index( words(j),'SRCX') /=0 ) THEN 
                            numfun = numfun + 1
                        ELSEIF ( index( words(j),'SRCZ') /=0) THEN
                            numfun = numfun + 1
                        ELSEIF ( index( words(j),'VSY') /=0 ) THEN 
                            numMech = numMech + 1
                        ELSEIF ( index( words(j),'VSX') /=0 ) THEN 
                            numMech = numMech + 1
                        ELSEIF ( index( words(j),'VSZ') /=0) THEN
                            numMech = numMech + 1
                        ELSE
                            !STOP ' err no define exitation X Y or Z '
                        ENDIF
                    ENDDO
                ENDIF    
            ENDIF
            
            IF (trim(words(i)) == 'TRAN') THEN
                DO j=i+1, neww-1,2
                    IF ( index( words(j),'STOP') /=0) THEN
                        ch_e=trim(words(j+1) ) !(k+1:))
                        Time=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' Time=',REAL(Time)
                    ELSEIF ( index(words(j),'STEP') /=0) THEN
                        ch_e=trim(words(j+1))
                        DT=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' DT=',REAL(DT)
                    ELSEIF ( index(words(j),'JUMP') /=0) THEN
                        ch_e=trim(words(j+1))
                        DTT=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' DTT=',REAL(DTT)
                    ENDIF
                ENDDO
            ENDIF
            
            IF (trim(words(i)) == 'SOLVER') THEN 
                DO j=i+1,neww-1
                    IF ( index(words(j),'TOL') /=0) THEN
                        ch_e=trim(words(j+1))
                        tolerance=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' tolerance=',REAL(tolerance)
                    ELSEIF ( index(words(j),'ITMAX') /=0) THEN
                        ch_e=trim(words(j+1))
                        itmax= nint(numeric(ch_e))
                        WRITE(*,'( a,i8,$ )') ' itmax=',itmax
                    ELSEIF ( index(words(j),'SOLV') /=0) THEN
                        solv = trim(words(j+1))
                        WRITE(*,'( a,a,$ )') ' solver=',solv
                    ELSEIF ( index(words(j),'DIR') /=0) THEN
                        files = trim(words(j+1))
                        WRITE(*,'( a,a,$ )') ' directory=',trim(files)
                    ELSEIF ( index(words(j),'BOUND') /=0) THEN
                        bound = trim(words(j+1))
                        WRITE(*,'( a,6a,$ )') ' boundary=',bound
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDIF !Name
    
    IF (index(st,'Structure Compression="ZLIB"') >= 1) THEN
        uncompress = .true.
    ELSEIF (index(st,'Structure Compression="ASCII_READABLE"') >= 1) THEN
        uncompress = .false.
    ENDIF
    m1=index(st,'<X_Voxels>')
    IF (m1>=1) THEN
        m2=index(st,'</X_Voxels>')
        ch_e=st(m1+10:m2-1)
        READ(ch_e,'(i3)')sdx
        WRITE(*,'( a,i3,$ )') ' sdx=',sdx
    ENDIF
    m1=index(st,'<Y_Voxels>')
    IF (m1>=1) THEN
        m2=index(st,'</Y_Voxels>')
        ch_e=st(m1+10:m2-1)
        READ(ch_e,'(i3)')sdy
        WRITE(*,'( a,i3,$ )') ' sdy=',sdy
    ENDIF
    m1=index(st,'<Z_Voxels>')
    IF (m1>=1 .and. m==0) THEN
        m2=index(st,'</Z_Voxels>')
        ch_e=st(m1+10:m2-1)
        READ(ch_e,'(i3)')sdz
        WRITE(*,'( a,i3 )') ' sdz=',sdz

        ALLOCATE(v(sdx*sdy*sdz), v_e(sdx*sdy*sdz), source=0_1)

        READ(1,'(a)',end=91,iostat=ios) ch_e ! ïðîïóñê ñòðîêè   <Data>
        IF ( uncompress ) THEN
        ! START UNCOMPRESS
        PRINT*,'start uncompress'
            ios=0
            open(2,file='compress.txt', iostat=ios) 
            
            idev=1
            DO j=1,sdz
                CALL readline(st,ch, ier, idev)             !================== 2 dev=1
                k=index(st,'CDATA')
                k2=index(st,']]></')
                IF (k2 == 0) THEN
                    STOP ' Error: long string CDATA'
                ENDIF
                k1=k+6
                ch=st(k1:k2-1)                        
                lst = len_trim(ch)
                WRITE(2, '( a )') ch(1:lst)
                m = m + lst
            ENDDO
            CLOSE(2)
            ! CALL system('python uncompress_zlib.py', i)
            call execute_command_line ('python uncompress_zlib.py', exitstat=i)
            
            
            PRINT *, "uncompress.txt created, status: ", i 
            open(2,file='uncompress.txt', iostat=ios)
            st=''; m=0
            idev = 2
            DO j=1,sdz
                CALL readline(st,ch, ier, idev)               !=======     3 idev= 2
                lst = len_trim(st)
                lst=lst-1
                DO i=1,lst
                    v(m+i) = index(letter, st(i:i))
                ENDDO
                m = m + lst
            ENDDO
            CLOSE(2)
            ! PRINT*, ' formed exitation array  v() '
            
            ! CALL system('del compress.txt; del uncompress.txt ',i)
            call execute_command_line ('del compress.txt; del uncompress.txt ', exitstat=i)
            PRINT *, "delete files uncompress.txt and compress.txt, status: ", i 
        ELSE
            idev=1
            DO j=1,sdz
                CALL readline(st,ch,ier, idev)              !--------------4  idev=1
                k=index(st,'CDATA')
                k2=index(st,']]></')
                IF (k2 == 0)   STOP ' Error: long string CDATA'
                k1=k+6
                ch=st(k1:k2-1)
                lst = len_trim(ch)
                DO i=1,lst
                    v(m+i) = index(letter, ch(i:i))
                ENDDO
                m = m + lst
            ENDDO
        ENDIF
    ENDIF

91 IF (ios/=0) EXIT

END DO convert_strings0
!=============================

nsub=maxval(v)        ! only the number of domains
Cells = sdx*sdy*sdz
kdz = sdx*sdy

j=0; k=1
DO i=1,Cells
    IF (v(i) == 0) THEN
        j = j + 1                ! new cell
        IF (j == 500000) THEN    ! 3500 2200  300000  
            j=0                  ! reset cells counter
            k = k + 1            ! new domain
        ENDIF
        v(i) = nsub + k          ! another one added for environment
    ENDIF
ENDDO
IF (j==0) k = k - 1              ! there was no enlargement of the cells
nsub_air = k

nsub_glob = nsub + nsub_air
 
PRINT '( "Cells= ",i9,  " numfun= ",i3, " numMechanical= ",i3, " nsub= ", i3, " nsub_air= ", i3 ," nsub+nsub_air=", i3 )', &
            Cells ,       numfun,        numMech,            nsub,          nsub_air,          nsub_glob
IF (numfun /=0 ) THEN
    ALLOCATE ( nameFun(numfun) )
    nameFun='';
ENDIF

IF ( numMech /=0 ) THEN
    ALLOCATE (  nameVmech(numMech) )
    nameVmech='';
    ALLOCATE (Vmech(numMech) ) 
    DO i=1,numMech
        Vmech(i)%vely = 0.d0
        Vmech(i)%args = 0
        Vmech(i)%nomsch = 0
        Vmech(i)%ex = ' '
        Vmech(i)%eqn = ' '
    ENDDO
ELSE
    ! PRINT*,  'var mech not use'
ENDIF

ALLOCATE ( typPHYS(nsub_glob), namePHYS(nsub_glob ) )
 typPHYS=''; namePHYS='';
ALLOCATE ( valPHYS(nsub_glob,5), source=0.d0 )

! DEFAULT param for ENVIRONMENT
IF (nsub_air /=0) THEN
    DO j=nsub+1, nsub_glob
        typPHYS(j)   = 'R     '
        namePHYS(j)  = 'AIR   '
        valPHYS(j,1) = 1.d0 ! D
    ENDDO
ENDIF

IF (  size_PHYS_C /=0  )  THEN
    PRINT*,'numbers domains dinamic with C = ',size_PHYS_C
    ALLOCATE ( PHYS_C(size_PHYS_C) )
ENDIF

if(numfun /=0) THEN
    ALLOCATE (Fun(numfun), fun_nod(numfun) )
    DO i=1,numfun
        Fun(i)%vely = 0.d0
        Fun(i)%args = 0
        Fun(i)%nomsch = 0
        Fun(i)%ex = ' '
        Fun(i)%eqn = ' '
        fun_nod(i)%vel_Vmech(1:3) = 0.d0; 
        fun_nod(i)%num_Vmech(1:3) = 0;  
        fun_nod(i)%move(1:3) = 0
    ENDDO
    fun_nod(:)%length(1)   = 0;    fun_nod(:)%length(2)   = 0;    fun_nod(:)%length(3)   = 0;   
    fun_nod(:)%shift(1)    = 0.d0; fun_nod(:)%shift(2)    = 0.d0; fun_nod(:)%shift(3)    = 0.d0;
    fun_nod(:)%Distance(1) = 0.d0; fun_nod(:)%Distance(2) = 0.d0; fun_nod(:)%Distance(3) = 0.d0; 
ENDIF

!           1       2        3      4       5         6        7       8          9       10     11       12
varconst=['PI  ','E   ','MU0 ','E0  ','DT  ', 'DX  ', 'DY  ', 'DZ  ', 'TIME', 'NX  ','NY  ','NZ  ']

valconst(1) = 3.1415926535897932384626433832795_8 ! 'PI  '
valconst(2) = 0.27182818284590451e+001_8          ! 'E   '
valconst(3) = 0.12566370964050292e-005_8          ! 'MU0 ' 
valconst(4) = 0.88541878176203908e-011_8          ! 'E0  '
valconst(5) = dt        ! 'DT  '
valconst(6) = delta(1)  ! 'DX  '
valconst(7) = delta(2)  ! 'DY  '
valconst(8) = delta(3)  ! 'DZ  '
valconst(9) = Time      ! 'TIME'
valconst(10)= sdx       ! 'NX  '
valconst(11)= sdy       ! 'NY  '
valconst(12)= sdz       ! 'NZ  '

REWIND(1)

!filling arrays
kp=0; ios=0; numfun=0; numMech=0
size_PHYS_C = 0;! size_PHYS_U0=0; size_PHYS_Ubnd = 0;

idev = 1
convert_strings1: &
DO WHILE ( 1==1)
    CALL readline(st,ch, ier, idev)              !-------------- 5
    IF (st(1:6)=='</VXC>' ) EXIT
    k0=index(st,'<Name>')
    IF ( k0 > 0 ) THEN 
        kp = kp + 1      ! count all, including lines with tran param etc
        m2=index(st,'</')
        st=st(k0+6:m2-1)
        i=1
        DO
            k = index( st(i:),'=')
            IF (k/=0) THEN
                i=i+k-1; st(i:i) = ' '
            ELSE
                EXIT
            ENDIF
        ENDDO
    
        CALL Upp (st, st)
        CALL string2words(trim(st),words,neww)
        
        DO i=2,neww
            IF ( words(i)(1:1)=='D' .and. kp <= nsub ) THEN  ! Key word D out environ
                valPHYS(kp,1) = evaluate( words(i+1) )
                typPHYS(kp) = 'R'
                namePHYS(kp) = trim(words(1))
                                  
                IF (i+1 == neww) THEN
                    PRINT '( "nsub= ",i3,  " valPHYS1= ",e10.3,  " typPHYS= ", a,   " namPHYS= ", a  )', &
                                  kp ,  valPHYS(kp,1),       trim(typPHYS(kp)),   trim(namePHYS(kp))
                    CYCLE ! it is last word
                ENDIF
                
                DO j=i+2, neww-1
                    IF ( words(j)(1:1)==  'C') THEN 
                        valPHYS(kp,2) = evaluate(words(j + 1))
                        IF (valPHYS(kp,2) /=0.d0) THEN
                            size_PHYS_C = size_PHYS_C + 1
                            typPHYS(kp) = trim(typPHYS(kp))//'C'
                            PHYS_C(size_PHYS_C)%numdom = kp
                            PHYS_C(size_PHYS_C)%valdom = 2.d0*valPHYS(kp,2)/dt ! double for force calc dinamics
                            
                            PRINT '( "nsub= ",i3,  " valPHYS2= ",e10.3,  " typPHYS= ", a,   " namPHYS= ", a  )', &
                                  kp ,           valPHYS(kp,2),          trim(typPHYS(kp)),   trim(namePHYS(kp))
                        ENDIF
                    ELSEIF ( index( words(j),'VEX')  /=0) THEN 
                        valPHYS(kp,3) = evaluate(words(j + 1))
                        PRINT '(  " Vex=", g10.3 )', valPHYS(kp,3)
                    ELSEIF ( index( words(j),'VEY')  /=0) THEN 
                        valPHYS(kp,4) = evaluate(words(j + 1))
                        PRINT '(  " Vey=", g10.3 )', valPHYS(kp,4)
                    ELSEIF ( index( words(j),'VEZ')  /=0) THEN 
                        valPHYS(kp,5) = evaluate(words(j + 1))
                        PRINT '(  " Vez=", g10.3 )',   valPHYS(kp,5)
                    ENDIF
                ENDDO
                
                IF ( index( words(i+2),'SRC') /=0 ) THEN 
                    DO j=i+2, neww-1
                        IF     ( index( words(j),'SRCY') /=0 ) THEN 
                            CALL calcVmech('Y')
                            ! PRINT '(   " Vsx= ",e10.3,    " Vsy= ", e10.3,  " Vsz= ", e10.3  )', &
                                     ! fun_nod(numfun)%vel_Vmech(1), fun_nod(numfun)%vel_Vmech(2), fun_nod(numfun)%vel_Vmech(3)
                        ELSEIF ( index( words(j),'SRCX') /=0) THEN 
                            CALL calcVmech('X')
                            ! PRINT '(   " Vsx= ",e10.3,    " Vsy= ", e10.3,  " Vsz= ", e10.3  )', &
                                     ! fun_nod(numfun)%vel_Vmech(1), fun_nod(numfun)%vel_Vmech(2), fun_nod(numfun)%vel_Vmech(3)
                        ELSEIF ( index( words(j),'SRCZ') /=0) THEN 
                            CALL calcVmech('D')
                            ! PRINT '(   " Vsx= ",e10.3,    " Vsy= ", e10.3,  " Vsz= ", e10.3  )', &
                                     ! fun_nod(numfun)%vel_Vmech(1), fun_nod(numfun)%vel_Vmech(2), fun_nod(numfun)%vel_Vmech(3)
                        ENDIF
                    ENDDO
                    ! PRINT '(   " Vsx= ",e10.3,    " Vsy= ", e10.3,  " Vsz= ", e10.3  )', &
                            ! fun_nod(numfun)%vel_Vmech(1), fun_nod(numfun)%vel_Vmech(2), fun_nod(numfun)%vel_Vmech(3)
                ENDIF 
            ELSEIF ( index( words(i),'FUNC') /=0  )  THEN  !   Key word FUNC
                DO L=1,numfun
                    If (  trim(words(i+1)) == trim(nameFun(L)) ) THEN ! 
                        j=i+1  
                        Fun(L)%eqn = trim(words(j+1))
        
                        DO WHILE (0 ==0  )  ! j - eto namex
                            Fun(L)%args = Fun(L)%args + 1
                            j=j+2 ! next x
                            IF (j+1 > neww) EXIT
                        ENDDO
                        Fun(L)%args = Fun(L)%args - 1
                        
                        ALLOCATE (Fun(L)%velx(1:Fun(L)%args),  Fun(L)%namex(1:Fun(L)%args) )
    
                        k=0; j=i+1 
                        DO WHILE (0 ==0  )  ! j - it is namex
                            k=k+1
                            IF ( k > Fun(L)%args ) EXIT
                            Fun(L)%namex(k) = words(j+2)(1:8)
                            Fun(L)%velx(k) = evaluate(words(j + 3))
                            j=j+2             ! next x
                            IF (j+1 > neww) EXIT
                        ENDDO
                    ENDIF 
                ENDDO ! numfun
                
                DO L=1,numMech
                    If (  trim(words(i+1)) == trim(nameVmech(L)) ) THEN ! 
                        j=i+1  
                        Vmech(L)%eqn = trim(words(j+1))
        
                        DO WHILE ( 0 ==0  )  ! j - it is  namex
                            Vmech(L)%args = Vmech(L)%args + 1
                            j=j+2            ! next x
                            IF (j+1 > neww) EXIT
                        ENDDO
                        Vmech(L)%args = Vmech(L)%args - 1
                        
                        ALLOCATE (Vmech(L)%velx(1:Vmech(L)%args),  Vmech(L)%namex(1:Vmech(L)%args) )
                        
                        k=0; j=i+1 
                        DO WHILE (0==0)         ! j - it is  namex
                            k=k+1
                            IF ( k > Vmech(L)%args ) EXIT
                            Vmech(L)%namex(k) = words(j+2)(1:8)
                            Vmech(L)%velx(k) = evaluate(words(j + 3))
                            j=j+2                ! next x
                            IF (j+1 > neww) EXIT
                        ENDDO
                    ENDIF 
                ENDDO ! numMech
            ELSEIF ( index( words(i),'BOUNDARY') /=0  )  THEN
                DO j=i+1, neww-1,2
                    SELECT CASE ( words(j)(1:3) )
                        CASE ( 'BXM') 
                            BND(1,1) =  evaluate(words(j + 1))
                        CASE ( 'BXP') 
                            BND(1,2) =  evaluate(words(j + 1))
                        CASE ( 'BYM') 
                            BND(2,1) =  evaluate(words(j + 1))
                        CASE ( 'BYP') 
                            BND(2,2) =  evaluate(words(j + 1))
                        CASE ( 'BZM') 
                            BND(3,1) =  evaluate(words(j + 1))
                        CASE ( 'BZP') 
                            BND(3,2) =  evaluate(words(j + 1))
                        CASE ( 'ALL') 
                            BND =  evaluate(words(j + 1))
                        CASE DEFAULT
                            PRINT*,'not recognized BOUNDARY, required BXM, BXP ... etc. ',trim(words(j))
                            STOP
                    END SELECT
                ENDDO
            ELSEIF ( index( words(i),'ENVIRON') /=0  )  THEN  !   ! Key word  ENVIRON - here not utulise !
            ! for ENVIRON calculate D, C etc
                DO j=i+1, neww-1
                    IF     (  words(j)(1:1)=='D' ) THEN 
                        valPHYS(nsub_glob,1) = evaluate(words(j + 1))
                    ELSEIF (  words(j)(1:1)=='C') THEN 
                        valPHYS(nsub_glob,2) = evaluate(words(j + 1))
                        IF (valPHYS(nsub_glob,2) /=0.d0) THEN
                            size_PHYS_C = size_PHYS_C + 1
                            typPHYS(nsub_glob) = trim(typPHYS(nsub_glob))//'C'
                            PHYS_C(size_PHYS_C)%numdom = nsub_glob 
                            PHYS_C(size_PHYS_C)%valdom = 2.d0*valPHYS(nsub_glob,2)/dt ! double for force calc dinamics
                        ENDIF
                    ELSEIF ( index( words(j),'VEX')  /=0) THEN 
                        valPHYS(nsub_glob,3) = evaluate(words(j + 1))
                    ELSEIF ( index( words(j),'VEY')  /=0) THEN 
                        valPHYS(nsub_glob,4) = evaluate(words(j + 1))
                    ELSEIF ( index( words(j),'VEZ')  /=0) THEN 
                        valPHYS(nsub_glob,5) = evaluate(words(j + 1))
                    ENDIF
                ENDDO
                PRINT '( "parameter ENVIRON: nsub= ",i3,  " valPHYS1= ",e10.3,  " typPHYS= ", a,   " namPHYS= ", a  )', &
                        nsub_glob ,  valPHYS(nsub_glob,1),       trim(typPHYS(nsub_glob)),   trim(namePHYS(nsub_glob))
            ENDIF  ! key word FUNC
        ENDDO  ! words
    ENDIF ! <Name>

92  IF (ios/=0) EXIT

END DO convert_strings1
    
CLOSE(1)

ALLOCATE (geoPHYS(sdx,sdy,sdz),  source=0_1  )
ALLOCATE (geoPHYS_C(sdx,sdy,sdz),  source=0  )
geoPHYS = reshape(v, shape =(/sdx,sdy,sdz/) ,  order = (/ 1, 2, 3/) )

! forming an array of node numbers for each np from the list PHYS_C(m)%numdom
IF (size_PHYS_C /= 0) THEN
    ! exclusion of inertial regions from boundaries for matrix BCG method under boundary conditions 'A' or 'N'
    IF (solv == 'BCG' .and. ( index(bound, 'A') /=0  .or.  index(bound, 'N')/=0 ) )  THEN
        DO m=1,size_PHYS_C
            np = PHYS_C(m)%numdom
            n=0
            DO k = 1,sdz;  DO j = 1,sdy;  DO i = 1,sdx
                n = n+1
                IF (i==1 .or.j==1 .or.k==1 .or.i==sdx .or.j==sdy .or.k==sdz ) THEN
                    IF (v(n) == np ) v(n) = nsub_glob
                ENDIF
            ENDDO; ENDDO; ENDDO
        ENDDO
    ENDIF
    
    ! ordinal numbers of cells in the conductivity region
    m=0
    DO L=1,size_PHYS_C
        np = PHYS_C(L)%numdom
        n=0
        do k=1,sdz; do j=1,sdy; do i=1,sdx
            n = n + 1; 
            if (v(n) == np) THEN
                m = m + 1
                geoPHYS_C(i,j,k) = 3*Cells + m  !local number
            endif
        enddo;enddo;enddo
    enddo
    
! the global number of cell numbers for each np from the list PHYS_C(i)%numdom
    DO i=1,size_PHYS_C
        np = PHYS_C(i)%numdom
        n = count(v == np)
        PHYS_C(i)%siznod = n
        ALLOCATE (PHYS_C(i)%nod(n), source=0   )
        k=0
        DO j=1,Cells
            IF (v(j) == np) THEN
                k=k+1
                PHYS_C(i)%nod(k) = j
            ENDIF
        ENDDO
    ENDDO
ENDIF

DEALLOCATE (v)

IF (numfun /=0 ) THEN
    n=0; 
    DO k=1,sdz
        DO j=1,sdy
            DO i=1,sdx
                np=geoPHYS(i,j,k); 
                n=n+1
                DO m=1,numfun
                    IF (Fun(m)% nomsch == np) THEN
                        v_e(n) = np
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    
    ! init
    NULLIFY(  sp_F, next_F   )
    fun_nod(:)%numnod_Fx=0;  fun_nod(:)%numnod_Fy=0; fun_nod(:)%numnod_Fz=0; 
    
    nn = 0;
    DO k=1,sdz
      DO j=1,sdy
        DO i=1,sdx
            nn =  nn + 1                ! current cell
            if(v_e(nn)/=0) THEN         ! only those cells where there are independent sources
                DO m=1,numfun
                    IF (Fun(m)%nomsch == v_e(nn)) THEN
                        IF     (Fun(m)%ex == 'X' ) THEN
                            fun_nod(m)%numnod_Fx = fun_nod(m)%numnod_Fx + 1
                            ALLOCATE(sp_F)
                            sp_F=tmp(m,nn,'X',next_F) 
                            next_F => sp_F
                        ELSEIF (Fun(m)%ex == 'Y' ) THEN 
                            fun_nod(m)%numnod_Fy = fun_nod(m)%numnod_Fy + 1
                            ALLOCATE(sp_F)
                            sp_F=tmp(m,nn,'Y',next_F) 
                            next_F => sp_F
                        ELSEIF (Fun(m)%ex == 'Z' ) THEN
                            fun_nod(m)%numnod_Fz = fun_nod(m)%numnod_Fz + 1
                            ALLOCATE(sp_F)
                            sp_F=tmp(m,nn,'D',next_F) 
                            next_F => sp_F
                        ELSE
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
         ENDDO !i
      ENDDO ! j
    ENDDO ! k
    
    num_nodx = 0; num_nody = 0; num_nodz =0 ; 
    DO m=1,numfun                                          ! ALLOCATION OF MEMORY
        IF (fun_nod(m)%numnod_Fx /= 0) THEN
            n = fun_nod(m)%numnod_Fx 
            ALLOCATE(fun_nod(m)%nods_Fx(n), source=0)  
            num_nodx = num_nodx + n
        ENDIF
        IF (fun_nod(m)%numnod_Fy /= 0) THEN  
            n = fun_nod(m)%numnod_Fy 
            ALLOCATE(fun_nod(m)%nods_Fy(n), source=0)  
            num_nody = num_nody + n
        ENDIF
        IF ( fun_nod(m)%numnod_Fz /= 0) THEN 
            n = fun_nod(m)%numnod_Fz 
            ALLOCATE(fun_nod(m)%nods_Fz(n), source=0) 
            num_nodz = num_nodz + n
        ENDIF

    ENDDO
    
    num_nodV = num_nodx + num_nody + num_nodz
    ! WRITE(*, '(  *(a,i6) )')  ' num_nodV=', num_nodV
    IF (num_nodV /=0 ) fun_nod(:)%numnod_Fx=0; fun_nod(:)%numnod_Fy=0; fun_nod(:)%numnod_Fz=0; 

    DO m=1,num_nodV 
        k=sp_F%jm  
        n=sp_F%im  
        ch = sp_F%ch   
        IF       (trim(ch) == 'X') THEN 
            fun_nod(k)%numnod_Fx = fun_nod(k)%numnod_Fx + 1
            k0 = fun_nod(k)%numnod_Fx
            fun_nod(k)%nods_Fx(k0) = n
        ELSEIF (trim(ch) == 'Y') THEN
            fun_nod(k)%numnod_Fy = fun_nod(k)%numnod_Fy + 1
            k0 = fun_nod(k)%numnod_Fy
            fun_nod(k)%nods_Fy(k0) = n + Cells
        ELSEIF (trim(ch) == 'Z') THEN
            fun_nod(k)%numnod_Fz = fun_nod(k)%numnod_Fz + 1
            k0 = fun_nod(k)%numnod_Fz
            fun_nod(k)%nods_Fz(k0) = n + 2*Cells
        ELSE
        ENDIF
        sp_F => sp_F%prec
    ENDDO
ENDIF ! numfun

! open (1, file='outdata.txt')
! WRITE(1, WRITEdata)

! WRITE(1, *) '  '
! WRITE(1, *) 'Description of physical domains'
! WRITE(1, '( "namePHYS= ", *(a10)    )' ) namePHYS
! WRITE(1, '( "typPHYS = ", *(a10)    )' ) typPHYS
! WRITE(1, '( "valPHYS1 = ", *(g10.3)  )' ) valPHYS(:,1)
! WRITE(1, '( "valPHYS2 = ", *(g10.3)  )' ) valPHYS(:,2)

! WRITE(1, *) '  '
! WRITE(1, *) 'Descriptions of functions', numfun
! DO i=1,numfun
    ! WRITE(1 , '( "num func= ",i3)' )  i
    ! WRITE(1 , '( "name = " a, " num domain =",i3, " direction =", a )' ) nameFun(i), Fun(i)%nomsch, Fun(i)%ex
    ! WRITE(1 , '( "equation = ",  a )') trim(Fun(i)%eqn)
    ! WRITE(1 , '( "names argum=", *(a5) )') (' ',trim(Fun(i)%namex(j)) , j= 1,Fun(i)%args)
    ! WRITE(1 , '( "values argum=", *(a, g10.3) )') (' ', Fun(i)%velx(j) , j= 1,Fun(i)%args)
    ! WRITE(1 , '( "velocityX=", g10.3," velocityY=", g10.3," velocityZ=", g10.3 )' ) &
                 ! fun_nod(i)%vel_Vmech(1), fun_nod(i)%vel_Vmech(2), fun_nod(i)%vel_Vmech(3)
    ! WRITE(1 , '( "num_VmechX=", i5," num_VmechY=", i5," num_VmechZ=", i5 )' )       &
                 ! fun_nod(i)%num_Vmech(1), fun_nod(i)%num_Vmech(2), fun_nod(i)%num_Vmech(3)
    ! WRITE(1, *) '  '
! ENDDO

! WRITE(1, *) '  '
! WRITE(1, *) 'Descriptions of mech', numMech
! DO i=1,numMech
    ! WRITE(1 , '( "num func= ",i3)' )  i
    ! WRITE(1 , '( "name = " a, " num domain =",i3, " direction =", a )' ) nameVmech(i), Vmech(i)%nomsch, Vmech(i)%ex
    ! WRITE(1 , '( "equation = ",  a )') trim(Vmech(i)%eqn)
    ! WRITE(1 , '( "names argum=", *(a5) )') (' ',trim(Vmech(i)%namex(j)) , j= 1,Vmech(i)%args)
    ! WRITE(1 , '( "values argum=", *(a, g10.3) )') (' ', Vmech(i)%velx(j) , j= 1,Vmech(i)%args)
    ! WRITE(1, *) '  '
! ENDDO
! close(1)
! STOP

IF (numfun /=0) THEN
    NULLIFY(sp_F, next_F )
    DEALLOCATE(nameFun)
ENDIF
IF (numMech /=0 ) DEALLOCATE( nameVmech)
DEALLOCATE( v_e, st, ch)

WRITE(*,'( a)') 'input data complet'

!========================================
CONTAINS

SUBROUTINE readline(line,ch,ier, idev)
    CHARACTER(len=:),ALLOCATABLE,INTENT(out) :: line ,ch
    INTEGER,INTENT(out)                      :: ier
    INTEGER                                  :: idev
    CHARACTER :: buffer
    line=''
    ier=0
    do
        READ( idev, '(a)', advance='no', iostat=ier ) buffer
        IF ( ier /= 0 ) THEN
            EXIT
        ELSE 
            line = line // buffer
        ENDIF
    END DO
    ch = line
END SUBROUTINE readline

FUNCTION evaluate (string) 
    REAL(8)  :: evaluate
    CHARACTER(*) string

    SELECT CASE (string(1:1))
        CASE ('"',"'",'`')          ! function specified in quotes
            m=LEN_TRIM(string)
            CALL initf(1); CALL parsef(1,trim(string(2:m-1)),varconst(1:12));
            evaluate = evalf(1,valconst(1:12))
        CASE DEFAULT
            evaluate = numeric( string )
    END SELECT
END FUNCTION evaluate 

SUBROUTINE calcVmech(ch)
CHARACTER(*) ch
    numfun = numfun + 1
    Fun(numfun)% nomsch = kp
    Fun(numfun)%ex = ch          
    nameFun(numfun) = trim(words(j+1))  

    DO n=1,6  
        IF ( j+1 + n+1 <= neww) THEN                       !   number for Vsx or Vsy Vsz
            IF    ( index( words(j+1 + n ),'VSX') /=0 ) THEN 
                fun_nod(numfun)%move(1) = 1                 ! it is move X
                SELECT CASE (words(j+1 + n+1)(1:1) )        ! analysis of the character after the = sign
                    CASE('A':'Z')                           ! function name specified
                        numMech = numMech + 1                ! this is for the calculation of the mechanical function itself
                        Vmech(numMech)% nomsch = kp
                        Vmech(numMech)%ex = 'X'
                        nameVmech(numMech) = trim(words(j+1+ n+1) )
                        fun_nod(numfun)%num_Vmech(1) = numMech  ! the number of the mechanical func is recorded
                        fun_nod(numfun)%vel_Vmech(1) = 0.d0
                    CASE DEFAULT                                ! the rest are either a number or quotes
                        fun_nod(numfun)%num_Vmech(1) = 0         ! writing the number of the mechanical function = 0
                        fun_nod(numfun)%vel_Vmech(1) = evaluate( words(j+1 + n+1) )
                END SELECT
            ELSEIF ( index( words(j+1 + n ),'VSY') /=0 ) THEN
                fun_nod(numfun)%move(2) = 1                 ! it is move Y
                SELECT CASE (words(j+1 + n+1)(1:1) )        ! analysis of the character after the = sign
                    CASE('A':'Z') ! function name specified
                        numMech = numMech + 1               ! this is for the calculation of the mechanical function itself
                        Vmech(numMech)% nomsch = kp
                        Vmech(numMech)%ex = 'Y'
                        nameVmech(numMech) = trim(words(j+1+ n+1) )
                        fun_nod(numfun)%num_Vmech(2) = numMech  ! the number of the mechanical func is recorded
                        fun_nod(numfun)%vel_Vmech(2) = 0.d0 
                    CASE DEFAULT                                 ! the rest are either a number or quotes
                        fun_nod(numfun)%num_Vmech(2) = 0         ! çwriting the number of the mechanical function = 0
                        fun_nod(numfun)%vel_Vmech(2) = evaluate( words(j+1 + n+1) )
                END SELECT
            ELSEIF ( index( words(j+1 + n ),'VSZ') /=0 ) THEN
                fun_nod(numfun)%move(3) = 1 ! it is move X
                SELECT CASE (words(j+1 + n+1)(1:1) ) 
                    CASE('A':'Z')
                        numMech = numMech + 1 
                        Vmech(numMech)% nomsch = kp
                        Vmech(numMech)%ex = 'D'
                        nameVmech(numMech) = trim(words(j+1+ n+1) )
                        
                        fun_nod(numfun)%num_Vmech(3) = numMech 
                        fun_nod(numfun)%vel_Vmech(3) = 0.d0 
                    CASE DEFAULT 
                        fun_nod(numfun)%num_Vmech(3) = 0 
                        fun_nod(numfun)%vel_Vmech(3) = evaluate( words(j+1 + n+1) )
                END SELECT
            ENDIF
        ENDIF
    ENDDO
END SUBROUTINE calcVmech

END SUBROUTINE vxc2data
