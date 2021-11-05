module INTEGRALS
    
!   Modules to use
    use IO
    
    implicit none
    contains
    
    !!!SUBROUTINES!!!
    
    subroutine overlap(NBASIS,ALPHA,SMAT)

!   Variables
    integer NBASIS
    real(8) SMAT(NBASIS,NBASIS),ALPHA(NBASIS)
    real(8), parameter :: pi = 4.0*atan(1.0)
!   Auxiliary variables
    integer i,j

!   Create the Matrix doing the "integrals".  
    do i = 1,NBASIS
        do j=1,NBASIS
            SMAT(i,j) = intOVERLAP(ALPHA(i),ALPHA(j))
        enddo
    enddo
    
    end subroutine overlap
    
    subroutine nuclear(ZATOM,NBASIS,ALPHA,VNUC)

!   Variables    
    integer NBASIS,ZATOM
    real(8) ZATOM1,VNUC(NBASIS,NBASIS),ALPHA(NBASIS)
    real(8), parameter :: pi = 4.0*atan(1.0)
!   Auxiliary variables
    integer i,j
    
    ZATOM1 = real(ZATOM)
!   The integrals itself to create the Matrix
    do i = 1,NBASIS
        do j=1,NBASIS
            VNUC(i,j) = intNUCLEAR(ZATOM1,ALPHA(i),ALPHA(j))
        enddo
    enddo
    
    end subroutine nuclear
    
    subroutine kinetic(NBASIS,ALPHA,TMAT)

!   Variables    
    integer NBASIS
    real(8) TMAT(NBASIS,NBASIS),ALPHA(NBASIS)
    real(8), parameter :: pi = 4.0*atan(1.0)
!   Auxiliary variables
    integer i,j

!   The integrals itself to create the Matrix
    do i = 1,NBASIS
        do j=1,NBASIS
            TMAT(i,j) = intKINETIC(ALPHA(i),ALPHA(j))
        enddo
    enddo
    
    end subroutine
    
    subroutine coulomb (NBASIS,ALPHA,PDENS,JMAT)
    
!   Variables    
    integer NBASIS
    real(8), dimension(NBASIS,NBASIS) :: PDENS, JMAT
    real(8), dimension(NBASIS) :: ALPHA
    real(8) intc, soma
    real(8), parameter :: pi = 4.0*atan(1.0)
!   Auxiliary variables
    integer u,n,s,l
    
    !  Four center integrals
    JMAT = 0.0
        do u=1,nbasis
            do n=1,nbasis
                do s=1,nbasis
                    do l=1,nbasis
                        JMAT(u,n) = JMAT(u,n) + PDENS(l,s)*intCOULOMB(ALPHA(u),ALPHA(n),ALPHA(s),ALPHA(l))
                    enddo
                enddo
            enddo
        enddo

    end subroutine coulomb
    
    !!!FUNCTIONS!!!
    
    real*8 function intOVERLAP(mi,ni)
    
    real(8) mi,ni
        
    intOVERLAP = (4.0*mi*ni/((mi+ni)**2.0))**(0.75)
    
    end function intOVERLAP
     
    real*8 function intNUCLEAR(ZATOM,mi,ni)
    
    real(8) ZATOM,mi,ni
    real(8), parameter :: pi = 4.0*atan(1.0)
        
    intNUCLEAR = -ZATOM * 2.0**(2.5)*pi**(-0.5) * (mi*ni)**(0.75)/(mi+ni)
    
    end function intNUCLEAR
    
    real*8 function intKINETIC(mi,ni)
    
    real(8) mi,ni
    
    intKINETIC = (6.0*2.0**(0.5)) * (mi**(0.75))*(ni**(1.75))/((mi+ni)**(1.5))&
    & *(1-(ni/(mi+ni)) )
    
    end function intKINETIC
    
    real*8 function intCOULOMB(mi,ni,sigma,lambda)
    
    real(8)  mi, ni, sigma, lambda
    real(8), parameter :: pi = 4.0*atan(1.0)
      
    intCOULOMB = (16.0/(pi**0.5))*(((mi*ni*sigma*lambda)**0.75)&
    &/((mi+ni)*(lambda+sigma)*(mi+ni+sigma+lambda)**0.5))
    
    end function intCOULOMB
    
end module