module SCF
    
    use IO
    use MAT
    
    implicit none
    contains
    
    subroutine buildX(NBASIS,SMAT,XMAT)
    
!   Variables
    integer NBASIS
    real(8), dimension(NBASIS,NBASIS) :: SMAT ,XMAT ,TXMAT
    real(8), dimension(NBASIS,NBASIS) :: AVAL, AVEC
!   Auxiliar variables
    real(8), dimension(NBASIS,NBASIS) :: saux, aux, aux2
    integer i,j
    
!   Start log on S diagonalization
    call writeLOG("S diagonalization matrices",3,"matrices.log")
    call writeLOG("",3,"matrices.log")
    
!   Diagonalization of the overlap matrix
    saux = SMAT
    call Jacobi(saux,NBASIS,AVAL,AVEC)
!   LOG
    call writeMATRIX(AVAL,NBASIS,NBASIS,"AutoValues S = ",3,"matrices.log")
    call writeMATRIX(AVEC,NBASIS,NBASIS,"AutoVectors S = ",3,"matrices.log")

!   Make s**-1/2 and 
    do i = 1,NBASIS
        AVAL(i,i) = AVAL(i,i)**(-0.5)
    enddo
!   LOG
    call writeMATRIX(AVAL,NBASIS,NBASIS,"AutoValues S**(-1/2) = ",3,"matrices.log")
    
!   LOG
    call writeMATRIX(transpose(AVEC),NBASIS,NBASIS,"AutoVector S Transposed = ",3,"matrices.log")
    
!   Make X = U(s**-1/2)Ut
    XMAT = 0.0
!   TXMAT = U*s**(-1/2)
    call multMATRIX(AVEC,NBASIS,NBASIS,AVAL,NBASIS,NBASIS,TXMAT)
    call writeMATRIX(TXMAT,NBASIS,NBASIS,"U * s**(-1/2) = ",3,"matrices.log")
!   XMAT = TXMAT*Ut
    call multMATRIX(TXMAT,NBASIS,NBASIS,transpose(AVEC),NBASIS,NBASIS,XMAT)
    call writeMATRIX(XMAT,NBASIS,NBASIS,"X (U*s-1/2*U) = ",3,"matrices.log")
    
!   Testing if the X matrix is corrected XSX = I
    call multMATRIX(XMAT,NBASIS,NBASIS,SMAT,NBASIS,NBASIS,AUX) !X*S
    call multMATRIX(AUX,NBASIS,NBASIS,transpose(XMAT),NBASIS,NBASIS,AUX2) !XS*X
    call writeMATRIX(AUX2,NBASIS,NBASIS,"XSX = ",3,"matrices.log")
    
    end subroutine buildX    
    
    subroutine buildFOCK(NBASIS,XMAT,FOCK,COEF)
     
    integer NBASIS
    real(8), dimension(NBASIS,NBASIS) :: XMAT,FOCK,FOCKL,FLAVAL,FLAVEC
    real(8), dimension(NBASIS,NBASIS) :: COEF
!   Auxiliar variables
    real(8), dimension(NBASIS,NBASIS) :: faux
    integer i
    character*100 TEXT
    
!   Generating the initial F' from F and X -> F'=X*FX
    FOCKL = 0.0
    FOCKL = matmul(matmul(transpose(XMAT),FOCK),XMAT)
    !LOG Fock'
    call writeMATRIX(FOCKL,NBASIS,NBASIS,"[F'] =",3,"matrices.log")

!   Diagonalize the F' matrix to obtain the C' matrix -> F'C' = C'e
!   C' = FLAVEC , e = FLAVAL
    faux = FOCKL
    call jacobi(faux,NBASIS,FLAVAL,FLAVEC)
    !LOG C'
    call writeMATRIX(FLAVEC,NBASIS,NBASIS,"[C'] =",3,"matrices.log")

    
!   Restoring the initial coeficients
    COEF = matmul(XMAT,FLAVEC)
    
!   Print the autovalues
    do i = 1,NBASIS
        write(TEXT,11) "Orbital ", i, " Energy = ", FLAVAL(i,i)
        call writeLOG(TEXT,3,"matrices.log")
    enddo
    call writeLOG("",3,"matrices.log")
    
    11  format(A8,1X,I3,1X,A10,F12.8)    
    
    end subroutine buildFOCK
    
    subroutine buildPDENS(NELEC,NBASIS,COEF,PDENS)
    
!   Variables
    integer NELEC,NBASIS
    real(8), intent(out), dimension(NBASIS,NBASIS) :: PDENS
    real(8), intent(in), dimension(NBASIS,NBASIS) :: COEF
    real(8), dimension(NBASIS,NBASIS) :: COEFT
!   Auxiliar variables
    integer i,j,k
    
    COEFT = COEF!transpose(COEF)
    
!   Build the density matrix
    ! If one electron
    PDENS = 0.0
    if (NELEC .EQ. 1) then
        do i = 1,NBASIS
            do j = 1,NBASIS
                    PDENS(i,j) = PDENS(i,j) + COEF(i,1)*COEF(j,1)
            enddo
        enddo
    ! If many electron closed-shell
    else
        do i = 1,NBASIS
            do j = 1,NBASIS
                do k = 1,NELEC/2
                    PDENS(i,j) = PDENS(i,j) + 2.0*COEF(i,k)*COEF(j,k) 
                enddo
            enddo
        enddo
    endif
       
    end subroutine buildPDENS
    
    subroutine buildENERGYDET(NBASIS,PDENS,TMAT,VNUC,JMAT,EKIN,ENUC,EJ) !Specify all the energy terms (detailed).
    
    !   Variables
    integer NBASIS
    real(8), dimension(NBASIS,NBASIS) :: TMAT, VNUC, JMAT, PDENS
    real(8) EKIN, ENUC, EJ
!   Auxiliar variables
    integer i,j
    character*100 TEXT
    
    EKIN = 0.0
    EJ = 0.0
    ENUC = 0.0
    do i = 1,NBASIS
        do j = 1,NBASIS
            EKIN = EKIN + PDENS(j,i)*(TMAT(i,j))
            ENUC = ENUC + PDENS(j,i)*(VNUC(i,j))
            EJ = EJ + PDENS(j,i)*((0.5)*JMAT(i,j))
        enddo
    enddo
   
    end subroutine buildENERGYDET

end module 