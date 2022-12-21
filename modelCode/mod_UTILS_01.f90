module UTILS_1
    use GLOBAL
    implicit none

contains

    subroutine READ_MODEL_CONSTANTS(MODEL_CONSTANTS, INPUT_NO)

        real(kind = DBL), dimension(:), intent(out) :: MODEL_CONSTANTS

        integer :: INPUT_NO

        integer :: CONSTANT_NO
        character(len = 100) :: CONSTANT_NAME
        real(kind = DBL) :: CONSTANT_VALUE

        !print *, 'PELAGIC MODEL CONSTANTS:'
        !'i6, a35, f18.5'
        do
            read(unit = INPUT_NO, fmt = *, end = 200) &
                 CONSTANT_NO, CONSTANT_NAME, CONSTANT_VALUE
            !print *, CONSTANT_NO, CONSTANT_VALUE
            MODEL_CONSTANTS(CONSTANT_NO) = CONSTANT_VALUE
        end do

        200 continue

    end subroutine READ_MODEL_CONSTANTS


    subroutine READ_MODEL_CONSTANTS_PRICE(MODEL_CONSTANTS, INPUT_NO)

        real(kind = DBL), dimension(:), intent(out) :: MODEL_CONSTANTS

        integer :: INPUT_NO
        integer :: NUM_CONSTANTS
        integer :: i

        read(unit = INPUT_NO, fmt = *) NUM_CONSTANTS

        do i = 1, NUM_CONSTANTS
            read(unit = INPUT_NO, fmt = *) MODEL_CONSTANTS(i)
        end do

    end subroutine READ_MODEL_CONSTANTS_PRICE


    integer function STRANGER(SUSPECTED_VALUE)

        !Checks for NaN and Inf
        
        real(kind = DBL) ::  SUSPECTED_VALUE, BIGNUMBER, SMALLVALUE, RATIO
        
        BIGNUMBER  = 1.0D300
        SMALLVALUE = 1.0D-30
        
        STRANGER=0
        
        if (dabs(SUSPECTED_VALUE) < SMALLVALUE) then
            return
        end if
        
        if (.not.(SUSPECTED_VALUE .lt. 0.D0).and..not.(SUSPECTED_VALUE .ge. 0.D0)) then
            STRANGER=1
        end if
        
        if ((SUSPECTED_VALUE .lt. 0.D0).and.(SUSPECTED_VALUE .ge. 0.D0)) then
            STRANGER=1
        end if

        if (SUSPECTED_VALUE > BIGNUMBER) then
            RATIO = BIGNUMBER / SUSPECTED_VALUE

            if(RATIO .eq. 0.D0) then
                STRANGER=1
            endif
        endif

        return
    end function stranger


    integer function STRANGERSD(SUSPECTED_VALUE)
          
        ! Cheks for NaN and Inf
        ! Input is real(kind = DBL)!
    
          real(kind = DBL) :: SUSPECTED_VALUE, BIGNUMBER, RATIO
          integer :: STRANGERD 
          
          BIGNUMBER=1.0D30
          STRANGERD=0
          
          if (.not.(SUSPECTED_VALUE .lt. 0.).and..not.(SUSPECTED_VALUE .ge. 0.)) then
              STRANGERD=1
          end if
          
          if ((SUSPECTED_VALUE .lt. 0.).and.(SUSPECTED_VALUE .ge. 0.D0)) then
              STRANGERD=1
          end if
          
          RATIO = BIGNUMBER/SUSPECTED_VALUE
          
          if(RATIO .eq. 0.) then
            STRANGERSD=1
          endif

          return
    end function STRANGERSD             


    character function INTEGER_TO_STRING(I_VALUE)
        integer, intent(in) :: I_VALUE
        character(len = 2048) :: INTERNAL_FILE
        write(INTERNAL_FILE, *) I_VALUE
        INTEGER_TO_STRING = trim(adjustl(INTERNAL_FILE))
    end function INTEGER_TO_STRING


    integer function CALCULATE_DAY_OF_YEAR(DAY_OF_SIMULATION, BASE_YEAR) result (DAY_OF_YEAR)

        double precision, intent(in) :: DAY_OF_SIMULATION
        integer         , intent(in) :: BASE_YEAR
        
        integer :: i
        integer :: YEAR
        integer :: NUM_DAYS_IN_BASE_YEAR
        integer :: NUM_DAYS_IN_YEAR
        
        integer :: SIM_DAY
        
        SIM_DAY               = int(DAY_OF_SIMULATION)
        NUM_DAYS_IN_BASE_YEAR = CALCULATE_NUM_DAYS_IN_YEAR(BASE_YEAR)

        if (SIM_DAY .le. NUM_DAYS_IN_BASE_YEAR) then
            DAY_OF_YEAR = DAY_OF_SIMULATION
        else
            YEAR    = BASE_YEAR + 1
            SIM_DAY = SIM_DAY - NUM_DAYS_IN_BASE_YEAR
            NUM_DAYS_IN_YEAR = CALCULATE_NUM_DAYS_IN_YEAR(YEAR)

            do while (SIM_DAY .gt. NUM_DAYS_IN_YEAR)
                NUM_DAYS_IN_YEAR = CALCULATE_NUM_DAYS_IN_YEAR(YEAR) 
                SIM_DAY          = SIM_DAY - NUM_DAYS_IN_YEAR
                YEAR             = YEAR + 1
            end do
        end if
        
        DAY_OF_YEAR = SIM_DAY
    end function CALCULATE_DAY_OF_YEAR


    integer function CALCULATE_NUM_DAYS_IN_YEAR(YEAR) result(NUM_DAYS_IN_YEAR)
        integer, intent(in) :: YEAR
            
        if (mod(YEAR, 4) > 0) then
            NUM_DAYS_IN_YEAR = 365
        else
            NUM_DAYS_IN_YEAR = 366
        end if

    end function CALCULATE_NUM_DAYS_IN_YEAR

end module UTILS_1