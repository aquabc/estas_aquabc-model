module AQUATIC_MODEL
    use GLOBAL
    use PELAGIC_BOX_MODEL
    use RESUSPENSION
    use BOTTOM_SEDIMENTS
    use INITIALIZE_PELAGIC_BOX_MODEL
    use TIME_SERIES

    implicit none

    type AQUATIC_MODEL_DS
        integer          :: BASE_YEAR
        real(kind = DBL) :: SIMULATION_START
        real(kind = DBL) :: SIMULATION_END
        integer          :: TIME_STEPS_PER_DAY
        real(kind = DBL) :: TIME_STEP
        integer          :: PRINT_INTERVAL
        integer          :: NUM_REPEATS
        integer          :: DAY_OF_YEAR

        type(PELAGIC_BOX_MODEL_DS) :: PELAGIC_BOX_MODEL_DATA
    end type AQUATIC_MODEL_DS

contains

    subroutine READ_AQUATIC_MODEL_INPUTS(AQUATIC_MODEL_DATA, IN_FILE, OUT_FILE)

        type(AQUATIC_MODEL_DS), intent(inout) :: AQUATIC_MODEL_DATA
        integer, intent(in) :: IN_FILE
        integer, intent(in) :: OUT_FILE

        integer :: i, j
        character(len = 2048) :: FILE_NAME
        real(kind = DBL)      :: TIME

        !READ DESCRIPTION LINES
        do i = 1, 5
            read(unit = IN_FILE, fmt = *)
        end do

        !READ BASIC MODEL SETUP INFORMATION
        read (unit = IN_FILE, fmt = *)
        read (unit = IN_FILE, fmt = *) AQUATIC_MODEL_DATA % BASE_YEAR

        read (unit = IN_FILE, fmt = *)
        read (unit = IN_FILE, fmt = *) AQUATIC_MODEL_DATA % SIMULATION_START

        write(unit = *, fmt = *) &
              'SIMULATION START  (Julian days) : ', AQUATIC_MODEL_DATA % SIMULATION_START

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) AQUATIC_MODEL_DATA % SIMULATION_END

        write(unit = *, fmt = *) &
              'SIMULATION END     (Julian days) : ', AQUATIC_MODEL_DATA % SIMULATION_END

        !Read the number of repeating simulations
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) AQUATIC_MODEL_DATA % NUM_REPEATS

        write(unit = *, fmt = *) &
              'NUMBER OF REPEATS (Julian days) : ', AQUATIC_MODEL_DATA % SIMULATION_END

        !Read time step
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) AQUATIC_MODEL_DATA % TIME_STEPS_PER_DAY

        AQUATIC_MODEL_DATA % TIME_STEP = 1.0D0 / float(AQUATIC_MODEL_DATA % TIME_STEPS_PER_DAY)

        write(unit = *, fmt = *) &
              'TIME STEP               (days) : ', AQUATIC_MODEL_DATA % TIME_STEP

        !Read the print interval
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) AQUATIC_MODEL_DATA % PRINT_INTERVAL

        write(unit = *, fmt = *) &
              'PRINT INTERVAL   (time steps) : ', AQUATIC_MODEL_DATA % PRINT_INTERVAL

        !Read the pelagic model output folder
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) PELAGIC_INPUT_FOLDER

        ! This is for linux
        !DEC$ IF DEFINED(_WIN32)
        print *,'Windows'
        !DEC$ ELSEIF DEFINED(__linux)
        print *,'Operating system is Linux'
        PELAGIC_INPUT_FOLDER = trim(PELAGIC_INPUT_FOLDER) // '/'
        ! End linux
        !DEC$ ELSE
        print *, 'Oops'
        !DEC$ ENDIF


        write(unit = *, fmt = *) &
              'PELAGIC_INPUT_FOLDER : ', trim(adjustl(PELAGIC_INPUT_FOLDER))

        !Read the pelagic model input file name
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) FILE_NAME



        open(unit   = IN_FILE + 2, &
             file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // trim(adjustl(FILE_NAME)), &
             status = 'OLD')

        !Read the pelagic model output folder
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) PELAGIC_OUTPUT_FOLDER

        ! This is for linux
        !DEC$ IF DEFINED(_WIN32)
        print *,'Windows'
        !DEC$ ELSEIF DEFINED(__linux)
        print *,'Operating system is Linux'
        PELAGIC_OUTPUT_FOLDER = trim(PELAGIC_OUTPUT_FOLDER) // '/'
        ! End linux
        !DEC$ ELSE
        print *, 'Oops'
        !DEC$ ENDIF

        write(unit = *, fmt = *) &
              'PELAGIC_OUTPUT_FOLDER : ', trim(adjustl(PELAGIC_OUTPUT_FOLDER))

        !Read the pelagic model inputs
        call READ_PELAGIC_BOX_MODEL_INPUTS &
             (AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA, IN_FILE + 2, IN_FILE + 4)

        close(IN_FILE + 2)
        close(IN_FILE + 4)

        TIME      = AQUATIC_MODEL_DATA % SIMULATION_START
        INIT_TIME = TIME
        call INIT_AQUATIC_MODEL(AQUATIC_MODEL_DATA, TIME)

        nkn = AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES

        if (AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS.ne. &
            nstate) then

            write(*,*) &
                 'The number of pelagic state variables are not compitable', &
                 ' with AQUABC pelagic module'

            write(*,*) &
                'The number of pelagic state variables                     : ', &
                AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

            write (*,*) 'The number of pelagic state variables requested by AQUABC : ', nstate
            stop("")
        end if

        if (AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % NUM_MODEL_CONSTANTS.ne. &
            nconst) then

            write(*,*) &
                 'The number of pelagic model constants are not compitable', &
                 ' with AQUABC pelagic module'

            write(*,*) &
                'The number of pelagic model constants                             : ', &
                AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % NUM_MODEL_CONSTANTS

            write (*,*) 'The number of pelagic model constants requested by AQUABC : ', nconst
            stop("")
        end if

        allocate(node_active          (nkn)                    )
        allocate(STATE_VARIABLES      (nkn,nstate)             )
        allocate(DERIVATIVES          (nkn,nstate)             )
        allocate(MODEL_CONSTANTS      (nconst)                 )
        allocate(DRIVING_FUNCTIONS    (nkn,n_driving_functions))
        allocate(FLAGS                (nflags)                 )
        allocate(PROCESS_RATES        (nkn,nstate, NDIAGVAR)   )
        allocate(SAVED_OUTPUTS        (nkn,n_saved_outputs)    )
        allocate(pH                   (nkn)                    )
        allocate(CHLA                 (nkn)                    )
        allocate(SURFACE_BOXES(nkn))

        allocate(FLUXES_TO_WATER_COLUMN       (nkn,nstate))
        allocate(FLUXES_OUTPUT_TO_WATER_COLUMN(nkn,nstate))

        do i = 1, nkn
            SURFACE_BOXES(i) = &
                AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % SURFACE_BOX
        end do

        MODEL_CONSTANTS = AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % MODEL_CONSTANTS(:, 1)

        allocate(DISSOLVED_FRACTIONS           (nkn, nstate      ))
        allocate(FRACTION_OF_DEPOSITION        (nkn, nstate      ))
        allocate(SETTLING_RATES                (nkn, nstate      ))
        allocate(NOT_DEPOSITED_FLUXES          (nkn, nstate      ))
        allocate(FLUXES                        (nkn, NUM_SED_VARS))
        allocate(SETTLING_VELOCITIES_OUTPUT    (nkn, nstate      ))
        allocate(EFFECTIVE_DISSLOVED_FRACTIONS (nkn, nstate      ))
        allocate(EFFECTIVE_DEPOSITION_FRACTIONS(nkn, nstate      ))
        allocate(DEPOSITION_AREA_RATIOS        (nkn, nstate      ))

        call INIT_PELAGIC_MODEL_CONSTANTS()

        ! Print out the water levels
        if (PRODUCE_ONLY_BINARY_PELAGIC_OUTPUT == 0) then

            open (unit   = 1000, &
                  file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // 'WATER_LEVELS.out', &
                  status = 'UNKNOWN')

            write(unit = 1000, fmt = '(2A10, 3A30)') &
                  '      TIME', '    BOX NO',  &
                  '         SURFACE ELEVATION (m)', '                  VOLUME (m^3)', &
                  '                     DEPTH (m)'

            open (unit   = 1001, &
                  file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // 'MASS_BALANCES.out', &
                  status = 'UNKNOWN')

            write(unit = 1001, fmt = '(3A10, 7A30)') &
                  '      TIME', '    BOX NO',  '    VAR_NO', &
                  '        ADVECTION (g/m^3/days)', '        DIFFUSION (g/m^3/days)', &
                  '         SETTLING (g/m^3/days)', '       MASS LOADS (g/m^3/days)', &
                  ' MASS WITHDRAWALS (g/m^3/days)', '         KINETICS (g/m^3/days)', &
                  '  SEDIMENT FLUXES (g/m^3/days)'
        else
            open(unit   = 1001, &
                 file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                          trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" // &
                          'MASS_BALANCES.bin', &
                 status = 'UNKNOWN', access = 'STREAM')
        end if

        if (PRODUCE_COCOA_OUTPUTS > 0) then
            open(unit   = 2020, &
                 file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                          trim(adjustl(COCOA_PELAGIC_OUTPUTS_FILENAME)), &
                 status = 'UNKNOWN')

            write(unit = 2020, fmt = '(3A10, 8A30)') &
                  '      TIME', '    BOX NO', '  LAYER NO', &
                  '         N_PEL_DENITRIFICATION', '                    N_FIXATION', &
                  '                   N_PEL_ASSIM', '                  N_PEL_EXCRET', &
                  '           N_PEL_DECOMP_OF_DET', '                   P_PEL_ASSIM', &
                  '                  P_PEL_EXCRET', '           P_PEL_DECOMP_OF_DET'
        end if
        ! -----------------------------------------------------------------------------------
        ! END OF INITIALIZATION OF THE WATER COLIMN MODEL
        ! -----------------------------------------------------------------------------------


        ! Read the resuspension option
        !
        ! Option 1 : Fully prescribed resuspension. The user supplies the time series
        !            for the resuspension velocities and resuspended concentrations.
        !
        ! Option 2 : Semi-prescribed resuspension. ESTAS calculates the time series
        !            for the resuspension velocities, but the user still supplies
        !            the resuspended concentrations. (TO BE IMPLEMENTED)
        !
        ! Option 3 : Full water column - sediment coupling

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) RESUSPENSION_OPTION
        write(unit = *, fmt = *) 'RESUSPENSION_OPTION : ', RESUSPENSION_OPTION

        SHUT_DOWN_SETTLING = 0

        if (RESUSPENSION_OPTION < 1) then
            CONSIDER_RESUSPENSION = 0
            write(unit = *, fmt = *) &
                  'RESUSPENSION_OPTION : ', RESUSPENSION_OPTION

            write(unit = *, fmt = *) 'Resuspension will not be considered'
        else
            !Read the resuspension model output folder
            read(unit = IN_FILE, fmt = *)
            read(unit = IN_FILE, fmt = *) RESUSPENSION_INPUT_FOLDER

            write(unit = *, fmt = *) &
                  'RESUSPENSION INPUT FOLDER : ', trim(adjustl(RESUSPENSION_INPUT_FOLDER))

            !Read the resuspension model input file name
            read(unit = IN_FILE, fmt = *)
            read(unit = IN_FILE, fmt = *) FILE_NAME

            !Read the resuspension model output folder
            read(unit = IN_FILE, fmt = *)
            read(unit = IN_FILE, fmt = *) RESUSPENSION_OUTPUT_FOLDER

            write(unit = *, fmt = *) &
                  'RESUSPENSION OUTPUT FOLDER : ', &
                  trim(adjustl(RESUSPENSION_OUTPUT_FOLDER))

            select case (RESUSPENSION_OPTION)

                case (1)
                    CONSIDER_RESUSPENSION = 1

                    write(unit = *, fmt = *) &
                        'Resuspension will be considered as in Option 1'

                    !Open the resuspension model input file
                    open(unit   = IN_FILE + 2, &
                         file   = trim(adjustl(RESUSPENSION_INPUT_FOLDER)) // &
			                      trim(adjustl(FILE_NAME)), &
                         status = 'OLD')

                    !Read the resuspension model inputs
                    call READ_RESUSPENSION_FILE_OPTION_1(IN_FILE + 2)
                    close(IN_FILE + 2)

                case (2)
                    CONSIDER_RESUSPENSION = 1
                    write(unit = *, fmt = *) &
                          'RESUSPENSION_OPTION : ', RESUSPENSION_OPTION

                    write(unit = *, fmt = *) &
                        'Resuspension will be considered as in Option 2'

                    !Open the resuspension model input file
                    open(unit   = IN_FILE + 2, &
                         file   = trim(adjustl(RESUSPENSION_INPUT_FOLDER)) // &
			                      trim(adjustl(FILE_NAME)), &
                         status = 'OLD')

                    call READ_RESUSPENSION_FILE_OPTION_2(IN_FILE + 2)
                    close(IN_FILE + 2)

                case (3)
                    CONSIDER_RESUSPENSION = 0

                    write(unit = *, fmt = *) &
                        'RESUSPENSION_OPTION : ', RESUSPENSION_OPTION

                    write(unit = *, fmt = *) &
                         'Resuspension will not be considered'

                    write(unit = *, fmt = *) &
                         'becuse resuspension option 3 is not implemented yet.'

                    close(IN_FILE + 2)

            end select
        end if

        ! Read the bottom sediments inputs
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) MODEL_BOTTOM_SEDIMENTS

        if (MODEL_BOTTOM_SEDIMENTS == 1) then
            read(unit = IN_FILE, fmt = *)
            read(unit = IN_FILE, fmt = *) NUM_PRESCRIBED_SEDIMENT_FLUX_SETS
            read(unit = IN_FILE, fmt = *)

            allocate(AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % &
                         SEDIMENT_FLUX_TS_NOS &
                             (NUM_PRESCRIBED_SEDIMENT_FLUX_SETS, nkn, nstate))

            allocate(AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % &
                         SEDIMENT_FLUX_TS_VAR_NOS &
                             (NUM_PRESCRIBED_SEDIMENT_FLUX_SETS, nkn, nstate))

            do i = 1, NUM_PRESCRIBED_SEDIMENT_FLUX_SETS
                read(unit = IN_FILE, fmt = *) FILE_NAME

                open(unit   = IN_FILE + 2, &
                     file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // &
                              trim(adjustl(FILE_NAME)), &
                     status = 'OLD')

                call READ_BOTTOM_SEDIMENTS_FLUXES_INPUTS &
                         (AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA, IN_FILE + 2, i)

                close(IN_FILE + 2)
            end do
        end if

        if (MODEL_BOTTOM_SEDIMENTS > 1) then
            if (CONSIDER_RESUSPENSION > 0) then
                write(unit = *, fmt = *) &
                      'Bottom sediments are not coupled with resuspension ' // &
                      'in this version of ESTAS-AQUABC. Program halted.'
                stop
            end if

            read(unit = IN_FILE, fmt = *)
            read(unit = IN_FILE, fmt = *) FILE_NAME

            open(unit   = IN_FILE + 2, &
                 file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // &
                          trim(adjustl(FILE_NAME)), &
                 status = 'OLD')

            call READ_BOTTOM_SEDIMENTS_MODEL_INPUTS(IN_FILE + 2)

            if (AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % &
                    ADVANCED_REDOX_SIMULATION.ne. &
                BOTTOM_SED_ADVANCED_REDOX_SIMULATION) then

                write(*, *) 'Pelagic model and bottom sediment model cannot have different'
                write(*, *) 'options for advanced redox simulation'

                if (AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA % &
                        ADVANCED_REDOX_SIMULATION > 0) then

                    write(*, *) 'Advanced redox simulation option is ON by the pelagic model'
                else
                    write(*, *) 'Advanced redox simulation option is OFF by the pelagic model'
                end if

                if (BOTTOM_SED_ADVANCED_REDOX_SIMULATION > 0) then

                    write(*, *) 'Advanced redox simulation option is ON by the bottom sediment model'
                else
                    write(*, *) 'Advanced redox simulation option is OFF by the bottom sediment model'
                end if
            end if

            call INIT_BSED_MODEL_CONSTANTS()

            !call INIT_BS(INIT_SED_STATE_VARS, nkn, NUM_SED_LAYERS, NUM_SED_VARS)


            open(unit   = 1021, &
                 file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                          trim(adjustl(BOTTOM_SEDIMENT_CONCENTRATIONS_FILENAME)), &
                 status = 'UNKNOWN')

            open(unit   = 1023, &
                 file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                          trim(adjustl(BOTTOM_SEDIMENT_FLUXES_FILENAME)), &
                 status = 'UNKNOWN')

            if (PRODUCE_COCOA_OUTPUTS > 0) then

                open(unit   = 2021, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl(COCOA_SEDIMENT_PROCESS_RATES_FILENAME)), &
                     status = 'UNKNOWN')

                write(unit  = 2021, fmt = '(3A10, 8A30)') &
                      '      TIME', '    BOX NO', '  LAYER NO', &
                      '         N_SED_DENITRIFICATION', '        N_SED_REMINERALIZATION', &
                      '        P_SED_REMINERALIZATION'


                open(unit   = 2022, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl(COCOA_SEDIMENT_BURIAL_RATES_FILENAME)) , &
                     status = 'UNKNOWN')

                write(unit  = 2022, fmt = '(3A10, 2A30)') &
                      '      TIME', '    BOX NO', '  LAYER NO', &
                      '                     N_BURRIAL', '                     P_BURRIAL'


                open(unit   = 2031, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl(COCOA_FLUXES_FROM_SEDIMENTS_FILENAME)) , &
                     status = 'UNKNOWN')

                write(unit  = 2031, fmt = '(3A10, 2A30)') &
                      '      TIME', '    BOX NO', '  LAYER NO', &
                      '                N_OUT_FROM_SED', '                P_OUT_FROM_SED'


                open(unit    = 2032, &
                     file    = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                               trim(adjustl(COCOA_FLUXES_TO_SEDIMENTS_FILENAME))   , &
                      status = 'UNKNOWN')

                write(unit  = 2032, fmt = '(3A10, 3A40)') &
                      '      TIME', '    BOX NO', '  LAYER NO', &
                      '    N_FLX_FROM_WATER_TO_SED_DUE_TO_DENIT', &
                      '              N_PARTICULATE_ORG_INTO_SED', &
                      '              P_PARTICULATE_ORG_INTO_SED'
            end if

            do i = 1, nkn
                do j = 1, NUM_SED_LAYERS
                    write(unit = 1021, fmt = '(F10.4,2I10,24F20.10)') &
                          TIME, i, j, INIT_SED_STATE_VARS(i, j, :)
                end do
            end do

            close(IN_FILE + 2)
        else
            FLUXES_TO_WATER_COLUMN        = 0.0D0
            FLUXES_OUTPUT_TO_WATER_COLUMN = 0.0D0
        end if
    end subroutine READ_AQUATIC_MODEL_INPUTS


    subroutine INIT_AQUATIC_MODEL(AQUATIC_MODEL_DATA, TIME)
        type(AQUATIC_MODEL_DS), intent(inout) :: AQUATIC_MODEL_DATA
        real(kind = DBL), intent(in) :: TIME

        integer :: i

        !INITIALIZE THE PELAGIC BOX MODEL
        call INIT_PELAGIC_BOX_MODEL(AQUATIC_MODEL_DATA % PELAGIC_BOX_MODEL_DATA, TIME)

    end subroutine INIT_AQUATIC_MODEL

end module AQUATIC_MODEL
