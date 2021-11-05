program DFT_ATOM
    
!   Program DFT_ATOM to calculate the eletronic energy of an atom using 
!   SCF calculations with DFT level of theory
!   Authors: Gabriel Libânio
!   Date: 1th semester of 2016

!   Modules to use    
    use IO
    use INTEGRALS
    use MAT
    use SCF
    use XC
   
!   Starting the program
    implicit none

!   Main Variables   
    integer ZATOM,CHARGE,NBASIS,NELEC
    real(8), allocatable, dimension(:) :: ALPHA
    real(8), allocatable, dimension(:,:) :: SMAT, VNUC, TMAT, JMAT, VXC
    real(8), allocatable, dimension(:,:) :: PDENS, HCORE, FOCK, COEF
!   For diagonalization
    real(8), allocatable, dimension(:,:) :: XMAT, AVEC !overlap diagonal
    real(8), allocatable, dimension(:,:) :: HCORED !HCORE diagonal
!   SCF variables
    integer nstep,maxstep
    real(8) crit, erro, ietot, etot, EFOCK, EXC, EHARTREE, EJ
    !For detailed energies
    real(8) EKIN,ENUC
!   For numerical integration
    integer NPOINTS
    real(8) STEPSIZE, INTPREC
    real(8) XRHO(1000)
!   Log writing variables and files
    integer UNITS
    character*100 TEXT
    character*30 FILES(3)
!   Date and Time variables
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: VALUES
!   Auxiliar Variables
    integer i,j,k
    real(8) ZATOM1, soma
    real(8), allocatable, dimension(:,:) :: aux, aux2 !Auxiliar
    
!!!!!   Create or reset all necessary output files
    UNITS = 3
    FILES = (/"output.log","matrices.log","density.log"/)
    do j = 1, UNITS
        call buildFILES(UNITS+1,FILES(j))
    end do
    
!!!!   Read the input information
    call inREAD(ZATOM,CHARGE,NELEC,MAXSTEP,CRIT,STEPSIZE,INTPREC,NBASIS,ALPHA)
    
!!! Start building CORE matrices
    call writeLOG("Core matrices",3,"matrices.log")
    call writeLOG("",3,"matrices.log")
 
    
!!!!!!!!   Build the overlap matrix (SMAT) !!!!!!!
!   Log write (TEXT,UNIT,FILE)
    TEXT = "Building the overlap matrix..."
    call writeLOG(TEXT,2,"output.log")
    allocate (SMAT(NBASIS,NBASIS))
    SMAT = 0.0
    call overlap(NBASIS,ALPHA,SMAT)
!   Log write    
    TEXT = "Done."
    call writeLOG(TEXT,2,"output.log")
!   Write the matrix on output writeMATRIX(MAT,M,N,TEXT,UNIT,FILE)    
    call writeMATRIX(SMAT,NBASIS,NBASIS,"Overlap [S] = ",3,"matrices.log")
!    pause
!!!!!!   Done building overlap matrix !!!!!!
 
    
!!!!!!   Build the nuclear potential matrix (VNUC) !!!!!
!   Log write (TEXT,UNIT,FILE)
    TEXT = "Building nuclear potential matrix..."
    call writeLOG(TEXT,2,"output.log")
    allocate(VNUC(NBASIS,NBASIS))
    VNUC = 0.0  
    call nuclear(ZATOM,NBASIS,ALPHA,VNUC)
!   Log write
    TEXT = "Done."
    call writeLOG(TEXT,2,"output.log")
    call writeMATRIX(VNUC,NBASIS,NBASIS,"Nuclear [Vnuc] = ",3,"matrices.log")
!!!!!!   Done building nuclear matrix !!!!!!!
 
    
!!!!!!!   Build the kinetic energy matrix (TMAT) !!!!!!
!   Log write (TEXT,UNIT,FILE)
    TEXT = "Building kinetic energy matrix..."
    call writeLOG(TEXT,2,"output.log")
    
    allocate(TMAT(NBASIS,NBASIS))
    TMAT = 0.0
    call kinetic(NBASIS,ALPHA,TMAT)   
!   Log write
    TEXT = "Done."
    call writeLOG(TEXT,2,"output.log")
!   Write the matrix on output writeMATRIX(MAT,M,N,TEXT,UNIT,FILE)    
    call writeMATRIX(TMAT,NBASIS,NBASIS,"Kinetic [T] = ",3,"matrices.log")
!!!!!   Done building kinetic energy matrix !!!!!!!

    
!!!!!!   Build the HCORE  !!!!!!
!   Log write
    call writeLOG("Building HCORE...",2,"output.log")   
    allocate(HCORE(NBASIS,NBASIS))
    HCORE = TMAT + VNUC   
!   Log write
    call writeMATRIX(HCORE,NBASIS,NBASIS,"HCORE [HC] = ",3,"matrices.log")
    call writeLOG("Done.",2,"output.log")
!!!!!!   Done building HCORE !!!!!!
 
    
!!!!!!   Diagonalize the overlap matrix and construct X matrix (X = S**-1/2) !!!!!!!
    ! LOG
    call writeLOG("Diagonalizing the overlap matrix",2,"output.log")        
    allocate(XMAT(NBASIS,NBASIS),AVEC(NBASIS,NBASIS))
!   Calculate the X matrix
    call buildX(NBASIS,SMAT,XMAT)  
    ! LOG
    call writeMATRIX(XMAT,NBASIS,NBASIS,"[X] = ",3,"matrices.log")
    call writeLOG("Done.",2,"output.log")
!!!!!!   Done overlap matrix diagonalization  !!!!!!!

    
!!!!!   Start of initial guess  !!!!!!!!
    call writeLOG("Initial Guess matrices",3,"matrices.log")
    call writeLOG("",3,"matrices.log")
!   Initial guess with P = 0 -> equivalent to neglect all electron-electron interactions
!   Fock aprox HCORE
    ! LOG
    call writeLOG("Building the initial guess to Fock and density matrices...",2,"output.log")    
    allocate(FOCK(NBASIS,NBASIS),COEF(NBASIS,NBASIS),PDENS(NBASIS,NBASIS),JMAT(NBASIS,NBASIS),VXC(NBASIS,NBASIS))
    FOCK = HCORE
    call buildFOCK(NBASIS,XMAT,FOCK,COEF)
    ! LOG
    call writeMATRIX(FOCK,NBASIS,NBASIS,"Initial FOCK = ",3,"matrices.log")
    call writeMATRIX(COEF,NBASIS,NBASIS,"Initial COEF = ",3,"matrices.log")
    call writeLOG("Initial guess done.",2,"output.log")

!!!!!!   End of initial guess !!!!!!    
    
!!!!!!!!!!!   STARTING THE SCF CYCLE !!!!!!!!!!!
    call writeLOG("",2,"output.log")    
    call writeLOG("Starting the SCF cycle...",2,"output.log")
    call writeLOG("SCF matrices",3,"matrices.log")
    call writeLOG("",3,"matrices.log")
    
!   Defining auxiliar parameters
    ietot = 100
    erro = 1
    nstep = 1

    do while ((nstep .LE. MAXSTEP) .AND. (erro .GT. CRIT))
        
        ! LOG SCF Matrices
        write(TEXT,01) "SCF Step", nstep
        call writeLOG(TEXT,3,"matrices.log")
        
        !!! Build density matrix
        call buildPDENS(NELEC,NBASIS,COEF,PDENS)
        ! LOG density
        call writeMATRIX(PDENS,NBASIS,NBASIS,"[P] =",3,"matrices.log")           
        !   Calculating the number of numerical integration steps and plot the density on plot file
        write(TEXT,01) "SCF Step", nstep
        call writeLOG(TEXT,4,"density.log")
        call calcRHO(PDENS,ALPHA,NBASIS,STEPSIZE,INTPREC,NPOINTS,XRHO)
        
        !!! Build Coulomb
        call coulomb(NBASIS,ALPHA,PDENS,JMAT)
        !LOG coulomb
        call writeMATRIX(JMAT,NBASIS,NBASIS,"[J]=",3,"matrices.log")
        
        !!! Calculate the Exchange Correlation potential and energy
        call calcEXC(PDENS,ALPHA,NBASIS,STEPSIZE,INTPREC,NPOINTS,XRHO,VXC,EXC)
        !LOG Exc & Vxc
        call writeMATRIX(VXC,NBASIS,NBASIS,"[Vxc] =",3,"matrices.log")
        
        !!! Build Fock
        FOCK = HCORE + JMAT + VXC
        call buildFOCK(NBASIS,XMAT,FOCK,COEF)
        !LOG Fock
        call writeMATRIX(FOCK,NBASIS,NBASIS,"[F] =",3,"matrices.log")
        !LOG COEF 
        call writeMATRIX(COEF,NBASIS,NBASIS,"[C] =",3,"matrices.log")
        !Calculate the Hartree and Coulomb Energies
        !call buildENERGY(NBASIS,PDENS,HCORE,JMAT,EHARTREE,EJ,EFOCK)
        
        !!! Calculate the Energy !!!!
        call buildENERGYDET(NBASIS,PDENS,TMAT,VNUC,JMAT,EKIN,ENUC,EJ)
        !!! Total energy
        ETOT = EKIN + ENUC + EJ + EXC
              
        !LOG SCF ENERGY
        call writeLOG("",2,"output.log")
        write(TEXT,01) "SCF cycle ", nstep
        call writeLOG(TEXT,2,"output.log")       
        write(TEXT,02) "EKINETIC =",EKIN    !Print electron kinetic energy
        call writeLOG(TEXT,2,"output.log")
        write(TEXT,02) "ENUCLEAR =",ENUC    !Print Coulomb atraction potential energy
        call writeLOG(TEXT,2,"output.log")
        write(TEXT,02) "ECOULOMB =",EJ      !Print Coulomb repulsion potential energy
        call writeLOG(TEXT,2,"output.log")
        write(TEXT,02) "EXCHANGE =",EXC     !Print exchange potential energy
        call writeLOG(TEXT,2,"output.log")
        write(TEXT,02) "    ETOTAL = ", ETOT !Print total electronic energy
        call writeLOG(TEXT,2,"output.log")
        
        !Develop SCF variables
        erro = abs(ietot - ETOT)
        ietot = ETOT
        nstep = nstep + 1    
    enddo
    
    call writeLOG("",2,"output.log")
    call writeLOG("The SCF has converged!",2,"output.log")
    call writeLOG("",2,"output.log")
    
!!!!!!!!!!! END OF SCF !!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!! Formats
01 format (A10,1X,I4)
02 format (A12,1X,F12.8)
03 format (A22,1X,F10.6)
     
!!!!!   End of the program
!   Date and Time
    call DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
10  format (A57,I2,"/",I2,"/",I4," at ",I2,":",I2,".",I2)   
    write (TEXT,10) "Ended at ", &
    & VALUES(3), VALUES(2), VALUES(1), VALUES(5), VALUES(6), VALUES(7)
!   Write log
    call writeLOG("The calculation ended smoothly bro!!!",2,"output.log")
    call writeLOG(TEXT,2,"output.log")
      
    deallocate (ALPHA,SMAT,VNUC,TMAT,HCORE,XMAT,AVEC,FOCK,PDENS,JMAT,VXC,COEF)
    
end program DFT_ATOM