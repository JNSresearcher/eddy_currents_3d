! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )

SUBROUTINE sprsBCGstabWR (valA, irow, jcol, n, b, x, tolerance, itmax, iter)
! Biconjugate Gradient Stabilized with restart method 
! This is the version if the right side (vector B) is not 0
! discussed in detail in my project https://github.com/JNSresearcher/SOLVERS_BCGSTAB_GMRES

INTEGER jcol(*), irow(n+1)
REAL(8) valA(*)
INTEGER iter, itmax,  n 
REAL(8):: alpha, beta, omega, tolerance, rr0, rr0_new, Bnorm, &
                    R(n), R0(n), P(n), AP(n), S(n), AS(n), B(n), X(n)
    iter=0 
    R = sprsAx(X)
    DO j=1,n
        R(j) = B(j) - R(j)
    ENDDO 
    R0 = R
    P = R

    Bnorm=norm2(b)

    IF (Bnorm == 0.d0) RETURN
    DO
        IF (iter > itmax) THEN
            PRINT*, norm2(R)
            EXIT
        ENDIF
        iter = iter + 1
        AP = sprsAx(P)
        rr0 = DOT_PRODUCT (R,R0)
        alpha = rr0/dot_product(AP,R0) 
        S = R - alpha*AP
        IF (norm2(S)/Bnorm < tolerance) THEN
        ! IF (norm2(S) < tolerance) THEN
            X = X + alpha*P;
            EXIT
        ENDIF
        AS = sprsAx(S)
        omega = dot_product(AS,S)/dot_product(AS,AS)
        X = X + alpha*P + omega*S
        R = S - omega*AS
        IF ( (norm2(R)/Bnorm) < tolerance) EXIT
        rr0_new = dot_product(R,R0)
        beta = (alpha/omega)*rr0_new/rr0
        P = R + beta*(P - omega*AP)
        IF ( (abs(rr0_new)/Bnorm) < tolerance) THEN
            R0 = R; P = R
        ENDIF
    ENDDO

CONTAINS

FUNCTION sprsAx  (V) 
    REAL(8)  :: sprsAx(n), V(n)
    INTEGER i1,i2
    DO  i=1,n
        i1=irow(i); i2=irow(i+1)-1
        sprsAx(i) =  dot_product(valA(i1:i2), V(jcol(i1:i2) ) )
    ENDDO 
END FUNCTION sprsAx 

END SUBROUTINE sprsBCGstabWR
