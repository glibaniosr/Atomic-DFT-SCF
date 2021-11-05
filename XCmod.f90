module XC
    
    use IO
    use MAT
     
    implicit none
    
    contains
    
    !!! SUBROUTINES !!!
    
    subroutine calcRHO(PDENS,ALPHA,NBASIS,STEPSIZE,INTPREC,NPOINTS,XRHO)
    
    !Variables
    integer NPOINTS,NBASIS
    real*8 INTPREC, STEPSIZE, integral
    real*8, dimension(NBASIS,NBASIS) :: PDENS
    real*8, dimension(NBASIS) :: ALPHA
    real*8, dimension(1000) :: XRHO, RXRHO
    real(8), parameter :: pi = 4.0*atan(1.0)
    !Auxiliar
    integer k,l,i
    real*8 rho,radrho,rr,xalpha,xbeta
    character*100 TEXT
    
    rho = 10.0
    rr = 0 + STEPSIZE
    NPOINTS = 1
    
1   format(F8.4,3X,F8.4,3X,F8.4,3X,F8.4)
2   format(F8.4,3X,F8.4,3X,F8.4)
    
!   Start plot log
    TEXT = "  r(a.u.)    Density   radDensity"
    call writeLOG(TEXT,4,"density.log")
        
!   Calculate Npoints       
    do while (rho .GT. INTPREC .AND. NPOINTS .LT. 1000)
        rho = 0.0
        radrho = 0.0
        do l=1,NBASIS
            xalpha = phiVALUE(ALPHA(l),rr)
            do k=1,NBASIS
                xbeta = phiVALUE(ALPHA(k),rr)
                rho = rho + PDENS(l,k)*xalpha*xbeta
                radrho = rho*4.0*pi*rr**2
            enddo
        enddo
        XRHO(NPOINTS) = rho
        RXRHO(NPOINTS) = radrho
        write(TEXT,2) rr, XRHO(NPOINTS), radrho
        call writeLOG(TEXT,4,"density.log")
        rr = rr + stepsize
        NPOINTS = NPOINTS + 1
    enddo
    
!   LOG write
    write(TEXT,03) "NPOINTS = ", NPOINTS
    call writeLOG(TEXT,4,"density.log")
    
!   Integrating RHO to give electron number    
    call simpson(RXRHO,NPOINTS,STEPSIZE,integral)
    write(TEXT,04) "Nelec (numeric) = ", integral
    call writeLOG(TEXT,4,"density.log")
    call writeLOG("",4,"density.log")
    
    03 format (A17,3x,I5)     
    04 format (A17,1X,F12.8)
    
    end subroutine calcRHO
    
    subroutine calcEXC(PDENS,ALPHA,NBASIS,STEPSIZE,INTPREC,NPOINTS,XRHO,VXC,EXC)
    
    !Variables
    integer NPOINTS, NBASIS
    real*8 STEPSIZE, INTPREC, EXC, integral, integral_result
    real*8, dimension(NBASIS,NBASIS) :: PDENS, Vxc
    real*8, dimension(NBASIS) :: ALPHA
    real*8, dimension(1000) :: intVXC, intEXC, XRHO
    real(8), parameter :: pi = 4.0*atan(1.0)
    !Auxiliar
    integer i,j,k,m
    real*8 rr,xalpha,xbeta
    character*100 TEXT   

    !Calculate the exchange potential & build the VXC matrix
    integral = 0.0
    do i = 1,NBASIS        
        do j = 1,NBASIS
            do k = 1,NPOINTS
                rr = k*STEPSIZE
                xalpha = phiVALUE(ALPHA(i),rr)
                xbeta = phiVALUE(ALPHA(j),rr)
                intVXC(k) = valueVXC(XRHO(k),rr,xalpha,xbeta)
            enddo
        call SIMPSON(intVXC,NPOINTS,STEPSIZE,integral)
        VXC(i,j) = integral
        enddo
    enddo
   
    !Calculate the Exchange energy
    integral = 0.0
    do k = 1,NPOINTS
        rr = k*STEPSIZE
        intEXC(k) = valueEXC(XRHO(k),rr)
    enddo
    call SIMPSON(intEXC,NPOINTS,STEPSIZE,integral)
    EXC = integral
        
    end subroutine
    
    !!!! FUNCTIONS !!!!
    
    real*8 function phiVALUE(alpha,rr)
    
    !Variables
    real(8), intent(in) :: alpha,rr
    !Auxiliar
    real(8), parameter :: pi = 4.0*atan(1.0) 
    real(8) Cte
    
    Cte = (2.0*alpha/pi)**(0.75)
    phiVALUE = Cte*exp(-alpha*(rr**2.0))
    
    end function phiVALUE
    
    real*8 function valueVXC(rho,rr,xalpha,xbeta)
    
    real(8), intent(in) :: rho,rr,xalpha,xbeta
    real(8) Cx
    real(8), parameter :: pi = 4.0*atan(1.0)
    
    Cx = (1.5)*(3.0/(4.0*pi))**(1.0/3.0) !Cx from Slater article
    !Cx = (3.0)*(3.0/(4.0*pi))**(1.0/3.0) !Cx from Becke's article on J. Chem. Phys. 2004
    
    valueVXC = xalpha*(-Cx*rho**(1.0/3.0))*xbeta*(4*pi)*(rr**2)
    
    end function valueVXC
    
    real*8 function valueEXC(rho,rr)
    
    real(8), intent(in) :: rho,rr
    real(8) Cx
    real(8), parameter :: pi = 4.0*atan(1.0)
    
    Cx = (3.0/4.0)*((3.0/pi)**(1.0/3.0)) !Cx from Slater article
    !Cx = (1.5)*(3.0/(4.0*pi))**(1.0/3.0) !Cx from Becke's article on J. Chem. Phys. 2004
    
    valueEXC = (-Cx*rho**(4.0/3.0))*(4*pi*rr**2)
    
    end function valueEXC
        
end module