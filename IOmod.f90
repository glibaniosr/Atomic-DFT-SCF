module IO
    
    implicit none
    contains
    
!   Subroutine to create necessary blank files
    subroutine buildFILES(UNIT,FILE)
    
    !   Variables
    integer UNIT
    character(len=*) :: FILE
    logical :: exist
    
    !   Checking if the file exists
    inquire(file=FILE, exist=exist) 
    if (exist) then
        open(unit=UNIT, file=FILE, status="replace", action="write")
    else
        open(unit=UNIT, file=FILE)
    endif
    
    close(UNIT)
    
    end subroutine buildFILES
    
!   Reading routines (I)
    subroutine inREAD(ZATOM,CHARGE,NELEC,MAXSTEP,CRIT,STEPSIZE,INTPREC,NBASIS,ALPHA)
!   This is the subroutine to read all the input information, including the atomic number (ZATOM),
!   the atomic charge (CHARGE), the number of electrons (NELEC), the number of gaussian primitives
!   (NBASIS) and its expoents (ALPHA).     
    
    integer ZATOM,CHARGE,NBASIS,NELEC,MAXSTEP !Atom data
    real(8) CRIT,STEPSIZE,INTPREC !Numerical Integrals
    real(8), allocatable :: ALPHA(:)
    character*100 :: TITLE,TEXT
    logical :: exist
!   For Date and Time
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: VALUES
!   Auxiliary variables
    integer i,j,k,l,m,n

!   Check if the input file is avaiable
    inquire(file="input.txt", exist=exist)
    if (exist==.false.) then
        print *, "There is no input.txt file in the current directory"
        print *, "Press ENTER to exit the program"
        pause
        stop
    else
        open(unit=1, file="input.txt")
    endif

!   Date and Time
    call DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
10  format (A57,I2,"/",I2,"/",I4," at ",I2,":",I2,".",I2)   
    write (TEXT,10) "Initializing DFT_ATOM routine made by Gabriel Libanio. ", &
    & VALUES(3), VALUES(2), VALUES(1), VALUES(5), VALUES(6), VALUES(7)
!   Write log
    call writeLOG(TEXT,2,"output.log")
    call writeLOG("Reading the input",2,"output.log")
    
!   STARTING THE REAL DEAL           
!   Read the integers
    read (1,'(A)') TITLE
    read (1,*) ZATOM, CHARGE
    ZATOM = abs(ZATOM)
    NELEC = ZATOM - CHARGE
    
!   Avoiding UHF calculations with NELEC = odd number different of 1
    if (mod(NELEC,2) .EQ. 1 .AND. NELEC .NE. 1 ) then
        TEXT = "Your number of electrons must be even (RHF calculation) or must be 1"
        print *, TEXT
        call writeLOG(TEXT,2,"output.log")
        print *, "Press enter to exit"
        pause
        stop
    else
        continue
    endif
    
!   Read the SCF parameters
    read (1,*) MAXSTEP, CRIT
    read (1,*) STEPSIZE, INTPREC

!   Reading the Basis parameters
    read(1,*) NBASIS
    !Read the exponents
    allocate (ALPHA(1:NBASIS))
    do i=1,NBASIS
        read (1,*) ALPHA(i)
    enddo
    
    call wlogHEADER(ZATOM,CHARGE,NELEC,NBASIS,ALPHA,TITLE)    
                  
    close(1)
    
    end subroutine inREAD
    
!   Writing routines (O)
    
    subroutine wlogHEADER(ZATOM,CHARGE,NELEC,NBASIS,ALPHA,TITLE)
!   Subroutine to write the HEADER of the log file    
    
    integer ZATOM,CHARGE,NBASIS,NELEC
    real(8), allocatable :: ALPHA(:)
    character*100 TITLE
    character*50 FILE
    logical :: exist
!   For Date and Time
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: VALUES
!   Auxiliary variables
    integer i
    
    open(unit=2, file="output.log", status="old", position="append", action="write")
    open(unit=3, file="matrices.log", status="old", position="append", action="write")

!   Formats
1   format (A8,1X,I3)
2   format (A25,1X,I3)
3   format (A11,1X,A25)
    write (2,*) " "
    write (2,*) "Input information"
    write (2,3) "Calculation ", TITLE
    write (2,2) "Atomic Number (Z) = ", ZATOM
    write (2,2) "Charge =", CHARGE
    write (2,2) "Number of Electrons =", NELEC
    write (2,2) "Number of Basis =", NBASIS
    write (2,*) "The basis S-type Gaussian exponents are: "
    do i=1,NBASIS
        write (2,*) ALPHA(i)
    enddo
    write (2,*) " "
!   End of input reading log
    
!   matrices.log file header
    call writeLOG("File containing the calculation matrices",3,"matrices.log")
    call writeLOG("",3,"matrices.log")
    
    close(2)
    close(3)
    
    end subroutine wlogHEADER
    
    subroutine writeLOG(TEXT,UNIT,FILE)
    
!   Variables
    integer UNIT
    character(len=*) :: TEXT,FILE
    logical :: exist
    
!   Checking if the file exists
    inquire(file=FILE, exist=exist) 
    if (exist) then
        open(unit=UNIT, file=FILE, status="old", position="append", action="write")
    else
        open(unit=UNIT, file=FILE)
    endif

!    write(UNIT,*) ""
    write(UNIT,"(A)",advance="yes"), TEXT
    
    close(UNIT)
    
    end subroutine writeLOG
    
    subroutine writeMATRIX(MAT,M,N,TEXT,UNIT,FILE)
    
    integer i,j
    integer, intent(in) :: M,N,UNIT  
    real(8), intent(in) :: MAT(M,N)
    character(len=*), intent(in) :: TEXT,FILE
    
    open(unit=UNIT,file=FILE,position="append")
    
    write(UNIT,"(A)"), TEXT
    
!    Writing the elements of MAT one by one)
    
    do i=1,M
        write(UNIT,"(100g15.8)") ( MAT(i,j), j=1,N )
    enddo
    
     write(UNIT,*) ""
    
    close(UNIT)
    
    end subroutine writeMATRIX
 
end module