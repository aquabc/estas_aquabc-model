module PELAGIC_BOX_MODEL

    use GLOBAL
    use TIME_SERIES
    use PELAGIC_BOX
    use PELAGIC_LINK
    use EXTERNAL_FORCING
    use INITIAL_CONDITIONS
    use BASIN
    use UTILS_1
    use COST_FUNCTION

    implicit none

    type PELAGIC_BOX_MODEL_DS
        integer :: NUM_PELAGIC_STATE_VARS
        integer :: NUM_BASINS
        integer :: NUM_BATHYMETRIES
        integer :: NUM_PELAGIC_BOXES
        integer :: NUM_PELAGIC_INIT_COND_SETS
        integer :: NUM_PELAGIC_ADVECTIVE_LINKS
        integer :: NUM_PELAGIC_DISPERSIVE_LINKS
        integer :: NUM_SETTLING_VELOCITIES
        integer :: NUM_OPEN_BOUNDARIES
        integer :: NUM_MASS_LOADS
        integer :: NUM_MASS_WITHDRAWALS

        integer :: NUM_FLOW_TS
        integer :: NUM_MIXING_TS
        integer :: NUM_FORCING_TS

        integer :: NUM_MODEL_CONSTANTS
        integer :: NUM_DRIVING_FUNCTIONS
        integer :: NUM_FLAGS
        integer :: NUM_PROCESS_RATES
        integer :: NUM_SAVED_OUTPUTS

        integer :: NUM_EXERGY_COMPONENTS
        integer :: PRODUCE_COST_FUNC

        character(len = 30), pointer, dimension(:) :: PELAGIC_STATE_VAR_NAMES

        real(kind = DBL), pointer, dimension(:) :: MEASUREMENT_ERRORS

        type(PELAGIC_BASIN_DS), pointer, dimension(:) :: BASINS
        type(BATHYMETRY_DS)   , pointer, dimension(:) :: BATHYMETRIES
        type(PELAGIC_BOX_DS)  , pointer, dimension(:) :: PELAGIC_BOXES

        type(PELAGIC_INIT_COND_SET_DS)  , pointer, dimension(:) :: PELAGIC_INIT_COND_SETS
        type(PELAGIC_ADVECTIVE_LINK_DS) , pointer, dimension(:) :: PELAGIC_ADVECTIVE_LINKS
        type(PELAGIC_DISPERSIVE_LINK_DS), pointer, dimension(:) :: PELAGIC_DISPERSIVE_LINKS

        type(OPEN_BOUNDARY_DS)  , pointer, dimension(:) :: OPEN_BOUNDARIES
        type(MASS_LOAD_DS)      , pointer, dimension(:) :: MASS_LOADS
        type(MASS_WITHDRAWAL_DS), pointer, dimension(:) :: MASS_WITHDRAWALS

        type(TIME_SERIE), pointer, dimension(:) :: FORCING_TS
        type(TIME_SERIE), pointer, dimension(:) :: FLOWS
        type(TIME_SERIE), pointer, dimension(:) :: MIXING_EXCHANGES
        type(TIME_SERIE), pointer, dimension(:) :: SETTLING_VELOCITIES

        real(kind = DBL), pointer, dimension(:,:) :: FORCINGS

        integer :: DAY_OF_YEAR
        !Exergy
        integer :: CALCULATE_PELAGIC_EXERGY
        integer :: CREATE_PELAGIC_EXERGY_OUTPUTS
        integer :: CREATE_PELAGIC_SAVED_OUTPUT
        integer :: START_REPEAT_NO_PEL_EX_OUTS
        integer :: CREATE_PELAGIC_ECOL_OUTPUT

        real(kind = DBL)    , pointer, dimension(:) :: EXERGY_COMPONENTS
        real(kind = DBL)    , pointer, dimension(:) :: AVERAGED_EXERGY_COMPONENTS
        character(len = 128), pointer, dimension(:) :: EXERGY_COMPONENT_NAMES
        real(kind = DBL) :: TOTAL_EXERGY

        !Internal data storage
        real(kind = DBL), pointer, dimension(:,:)   :: VOLUME_DERIVS
        real(kind = DBL), pointer, dimension(:,:)   :: SALT_ADVECTION_DERIVS
        real(kind = DBL), pointer, dimension(:,:)   :: TEMP_ADVECTION_DERIVS
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_ADVECTION_DERIVS
        real(kind = DBL), pointer, dimension(:,:)   :: SALT_DISPERSION_DERIVS
        real(kind = DBL), pointer, dimension(:,:)   :: TEMP_DISPERSION_DERIVS
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_DISPERSION_DERIVS
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_SETTLING_DERIVS
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_MASS_LOAD_DERIVS
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_MASS_WITHDRAWAL_DERIVS
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_KINETIC_DERIVS
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_PRESCRIBED_SEDIMENT_FLUX_DERIVS
        real(kind = DBL), pointer, dimension(:,:)   :: MODEL_CONSTANTS

        !OUTPUT STORAGE
        real(kind = DBL), pointer, dimension(:,:,:) :: ECOL_RESULTS
        real(kind = DBL), pointer, dimension(:,:,:) :: AREA_BASED_ECOL_RESULTS
        real(kind = DBL), pointer, dimension(:,:,:) :: SAVED_OUTPUT_RESULTS
        real(kind = DBL), pointer, dimension(:,:,:) :: EXERGY_RESULTS

        !MEASURED VALUES
        type(MEASURED_VALUE_DS), pointer, dimension(:) :: MEASURED_VALUES
        integer :: PRODUCE_COST_FUNCTION

        !SAVED OUTPUTS
        character(len = 128), pointer, dimension(:) :: SAVED_OUTPUT_NAMES

        !PROCESSES
        real(kind = DBL), pointer, dimension(:,:,:) :: PROCESS_RATE_RESULTS

        !PELAGIC OUTPUT INFORMATION
        integer, pointer, dimension(:) :: PRODUCE_PEL_STATE_VAR_OUTPUTS
        integer, pointer, dimension(:) :: PRODUCE_PEL_PROCESS_RATE_OUTPUTS
        integer, pointer, dimension(:) :: PRODUCE_PEL_MASS_BALANCE_OUTPUTS
        integer :: CREATE_STATE_VARIABLE_OUTPUT

        integer, allocatable, dimension(:, :, :) :: SEDIMENT_FLUX_TS_NOS
        integer, allocatable, dimension(:, :, :) :: SEDIMENT_FLUX_TS_VAR_NOS

        integer, pointer, dimension(:) :: ADVECTION_ON
        integer, pointer, dimension(:) :: DIFFUSION_ON
        integer, pointer, dimension(:) :: SETTLING_ON
        integer, pointer, dimension(:) :: STATE_VAR_OUTPUT_TYPES
        integer, pointer, dimension(:) :: STATE_VAR_INITIAL_CONDITION_TYPES

        !PELAGIC OPTIONS
        integer :: ZOOPLANKTON_OPTION
        integer :: ADVANCED_REDOX_SIMULATION
        integer :: LIGHT_EXTINCTION_OPTION
        integer :: CYANO_BOUYANT_STATE_SIMULATION
        integer :: CONSIDER_NON_OBLIGATORY_FIXERS
        integer :: CONSIDER_NOSTOCALES
    end type PELAGIC_BOX_MODEL_DS


contains


    subroutine ALLOC_PELAGIC_BOX_MODEL_DATA &
               (PELAGIC_BOX_MODEL_DATA      ,  &
                NUM_PELAGIC_STATE_VARS      , NUM_MODEL_CONSTANTS         , &
                NUM_BASINS                  , NUM_BATHYMETRIES            , &
                NUM_PELAGIC_BOXES           , NUM_PELAGIC_INIT_COND_SETS  , &
                NUM_PELAGIC_ADVECTIVE_LINKS , NUM_PELAGIC_DISPERSIVE_LINKS, &
                NUM_SETTLING_VELOCITIES     , NUM_OPEN_BOUNDARIES         , &
                NUM_FLOW_TS                 , NUM_MIXING_TS               , &
                NUM_MASS_LOADS              , NUM_MASS_WITHDRAWALS        , &
                NUM_FORCING_TS)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: NUM_PELAGIC_STATE_VARS
        integer, intent(in) :: NUM_MODEL_CONSTANTS
        integer, intent(in) :: NUM_BASINS
        integer, intent(in) :: NUM_BATHYMETRIES
        integer, intent(in) :: NUM_PELAGIC_BOXES
        integer, intent(in) :: NUM_PELAGIC_INIT_COND_SETS
        integer, intent(in) :: NUM_PELAGIC_ADVECTIVE_LINKS
        integer, intent(in) :: NUM_PELAGIC_DISPERSIVE_LINKS
        integer, intent(in) :: NUM_SETTLING_VELOCITIES
        integer, intent(in) :: NUM_OPEN_BOUNDARIES
        integer, intent(in) :: NUM_FLOW_TS
        integer, intent(in) :: NUM_MIXING_TS
        integer, intent(in) :: NUM_MASS_LOADS
        integer, intent(in) :: NUM_MASS_WITHDRAWALS
        integer, intent(in) :: NUM_FORCING_TS

        integer :: i

        allocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_STATE_VAR_NAMES(NUM_PELAGIC_STATE_VARS))
        allocate(PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES(NUM_PELAGIC_STATE_VARS))
        allocate(PELAGIC_BOX_MODEL_DATA % MEASUREMENT_ERRORS(NUM_PELAGIC_STATE_VARS))

        allocate(PELAGIC_BOX_MODEL_DATA % MODEL_CONSTANTS(NUM_MODEL_CONSTANTS, NUM_PELAGIC_BOXES))
        allocate(PELAGIC_BOX_MODEL_DATA % BASINS(NUM_BASINS))
        allocate(PELAGIC_BOX_MODEL_DATA % BATHYMETRIES(NUM_BATHYMETRIES))
        allocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(NUM_PELAGIC_BOXES))
        allocate(PELAGIC_BOX_MODEL_DATA % VOLUME_DERIVS(NUM_PELAGIC_BOXES, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % SALT_ADVECTION_DERIVS(NUM_PELAGIC_BOXES, 1))
        allocate(PELAGIC_BOX_MODEL_DATA % TEMP_ADVECTION_DERIVS(NUM_PELAGIC_BOXES, 1))
        allocate(PELAGIC_BOX_MODEL_DATA % ECOL_ADVECTION_DERIVS(NUM_PELAGIC_BOXES, NUM_PELAGIC_STATE_VARS, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % SALT_DISPERSION_DERIVS(NUM_PELAGIC_BOXES, 1))
        allocate(PELAGIC_BOX_MODEL_DATA % TEMP_DISPERSION_DERIVS(NUM_PELAGIC_BOXES, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % ECOL_DISPERSION_DERIVS               &
                     (NUM_PELAGIC_BOXES, NUM_PELAGIC_STATE_VARS, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % ECOL_SETTLING_DERIVS                 &
                     (NUM_PELAGIC_BOXES, NUM_PELAGIC_STATE_VARS, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % ECOL_MASS_LOAD_DERIVS                &
                     (NUM_PELAGIC_BOXES, NUM_PELAGIC_STATE_VARS, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % ECOL_MASS_WITHDRAWAL_DERIVS          &
                     (NUM_PELAGIC_BOXES, NUM_PELAGIC_STATE_VARS, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % ECOL_KINETIC_DERIVS                  &
                     (NUM_PELAGIC_BOXES, NUM_PELAGIC_STATE_VARS, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % ECOL_PRESCRIBED_SEDIMENT_FLUX_DERIVS &
                     (NUM_PELAGIC_BOXES, NUM_PELAGIC_STATE_VARS, 1))

        allocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_INIT_COND_SETS               &
                     (NUM_PELAGIC_INIT_COND_SETS))

        do i = 1, NUM_PELAGIC_INIT_COND_SETS
            call ALLOC_PELAGIC_INIT_COND_SET&
                     (PELAGIC_BOX_MODEL_DATA % PELAGIC_INIT_COND_SETS(i), &
                      NUM_PELAGIC_STATE_VARS)
        end do

        allocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_ADVECTIVE_LINKS  &
                     (NUM_PELAGIC_ADVECTIVE_LINKS))

        allocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_DISPERSIVE_LINKS &
                    (NUM_PELAGIC_DISPERSIVE_LINKS))

        allocate(PELAGIC_BOX_MODEL_DATA % FLOWS(NUM_FLOW_TS))

        allocate(PELAGIC_BOX_MODEL_DATA % MIXING_EXCHANGES(NUM_MIXING_TS))

        allocate(PELAGIC_BOX_MODEL_DATA % SETTLING_VELOCITIES(NUM_SETTLING_VELOCITIES))

        allocate(PELAGIC_BOX_MODEL_DATA % OPEN_BOUNDARIES(NUM_OPEN_BOUNDARIES))

        do i = 1, NUM_OPEN_BOUNDARIES
            call ALLOC_OPEN_BOUNDARY_DATA &
                     (PELAGIC_BOX_MODEL_DATA % OPEN_BOUNDARIES(i), NUM_PELAGIC_STATE_VARS)
        end do

        allocate(PELAGIC_BOX_MODEL_DATA % MASS_LOADS(NUM_MASS_LOADS))

        do i = 1, NUM_MASS_LOADS
            call ALLOC_MASS_LOAD_DATA &
                     (PELAGIC_BOX_MODEL_DATA % MASS_LOADS(i), NUM_PELAGIC_STATE_VARS)
        end do

        allocate(PELAGIC_BOX_MODEL_DATA % MASS_WITHDRAWALS(NUM_MASS_WITHDRAWALS))

        do i = 1, NUM_MASS_WITHDRAWALS
            call ALLOC_MASS_WITHDRAWAL_DATA &
                     (PELAGIC_BOX_MODEL_DATA % MASS_WITHDRAWALS(i), NUM_PELAGIC_STATE_VARS)
        end do

        allocate(PELAGIC_BOX_MODEL_DATA % FORCING_TS(NUM_FORCING_TS))

        allocate(PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS   (NUM_PELAGIC_BOXES))
        allocate(PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_PROCESS_RATE_OUTPUTS(NUM_PELAGIC_BOXES))
        allocate(PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_MASS_BALANCE_OUTPUTS(NUM_PELAGIC_BOXES))

        allocate(PELAGIC_BOX_MODEL_DATA % ADVECTION_ON(NUM_PELAGIC_STATE_VARS))
        allocate(PELAGIC_BOX_MODEL_DATA % DIFFUSION_ON(NUM_PELAGIC_STATE_VARS))
        allocate(PELAGIC_BOX_MODEL_DATA % SETTLING_ON (NUM_PELAGIC_STATE_VARS))

        allocate(PELAGIC_BOX_MODEL_DATA % &
                 STATE_VAR_OUTPUT_TYPES           (NUM_PELAGIC_STATE_VARS))

        allocate(PELAGIC_BOX_MODEL_DATA % &
                 STATE_VAR_INITIAL_CONDITION_TYPES(NUM_PELAGIC_STATE_VARS))

        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS       = NUM_PELAGIC_STATE_VARS
        PELAGIC_BOX_MODEL_DATA % NUM_MODEL_CONSTANTS          = NUM_MODEL_CONSTANTS
        PELAGIC_BOX_MODEL_DATA % NUM_BASINS                   = NUM_BASINS
        PELAGIC_BOX_MODEL_DATA % NUM_BATHYMETRIES             = NUM_BATHYMETRIES
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES            = NUM_PELAGIC_BOXES
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_INIT_COND_SETS   = NUM_PELAGIC_INIT_COND_SETS
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_ADVECTIVE_LINKS  = NUM_PELAGIC_ADVECTIVE_LINKS
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_DISPERSIVE_LINKS = NUM_PELAGIC_DISPERSIVE_LINKS
        PELAGIC_BOX_MODEL_DATA % NUM_FLOW_TS                  = NUM_FLOW_TS
        PELAGIC_BOX_MODEL_DATA % NUM_MIXING_TS                = NUM_MIXING_TS
        PELAGIC_BOX_MODEL_DATA % NUM_SETTLING_VELOCITIES      = NUM_SETTLING_VELOCITIES
        PELAGIC_BOX_MODEL_DATA % NUM_OPEN_BOUNDARIES          = NUM_OPEN_BOUNDARIES
        PELAGIC_BOX_MODEL_DATA % NUM_MASS_LOADS               = NUM_MASS_LOADS
        PELAGIC_BOX_MODEL_DATA % NUM_MASS_WITHDRAWALS         = NUM_MASS_WITHDRAWALS
        PELAGIC_BOX_MODEL_DATA % NUM_FORCING_TS               = NUM_FORCING_TS

        do i = 1, NUM_PELAGIC_BOXES
            call ALLOC_PELAGIC_BOX_DATA &
                    (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i), NUM_PELAGIC_STATE_VARS)
        end do

    end subroutine ALLOC_PELAGIC_BOX_MODEL_DATA


    subroutine DEALLOC_PELAGIC_BOX_MODEL_DATA(PELAGIC_BOX_MODEL_DATA)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA

        integer :: i

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_BASINS
            call DEALLOCATE_BASIN_DATA(PELAGIC_BOX_MODEL_DATA % BASINS(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % BASINS)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_BATHYMETRIES
            call DEALLOCATE_BATHYMETRY_DATA(PELAGIC_BOX_MODEL_DATA % BATHYMETRIES(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % BATHYMETRIES)
        deallocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_STATE_VAR_NAMES)
        deallocate(PELAGIC_BOX_MODEL_DATA % MODEL_CONSTANTS)
        deallocate(PELAGIC_BOX_MODEL_DATA % VOLUME_DERIVS)

        deallocate(PELAGIC_BOX_MODEL_DATA % SALT_ADVECTION_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % TEMP_ADVECTION_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % ECOL_ADVECTION_DERIVS)

        deallocate(PELAGIC_BOX_MODEL_DATA % SALT_DISPERSION_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % TEMP_DISPERSION_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % ECOL_DISPERSION_DERIVS)

        deallocate(PELAGIC_BOX_MODEL_DATA % ECOL_SETTLING_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % ECOL_MASS_LOAD_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % ECOL_MASS_WITHDRAWAL_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % ECOL_KINETIC_DERIVS)
        deallocate(PELAGIC_BOX_MODEL_DATA % ECOL_PRESCRIBED_SEDIMENT_FLUX_DERIVS)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            call DEALLOC_PELAGIC_BOX_DATA(PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_INIT_COND_SETS
            call DEALLOC_PELAGIC_INIT_COND_SET(PELAGIC_BOX_MODEL_DATA % PELAGIC_INIT_COND_SETS(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_INIT_COND_SETS)
        deallocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_ADVECTIVE_LINKS)
        deallocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_DISPERSIVE_LINKS)
        deallocate(PELAGIC_BOX_MODEL_DATA % FLOWS)
        deallocate(PELAGIC_BOX_MODEL_DATA % MIXING_EXCHANGES)
        deallocate(PELAGIC_BOX_MODEL_DATA % SETTLING_VELOCITIES)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_OPEN_BOUNDARIES
            call DEALLOC_OPEN_BOUNDARY_DATA(PELAGIC_BOX_MODEL_DATA % OPEN_BOUNDARIES(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % OPEN_BOUNDARIES)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_MASS_LOADS
            call DEALLOC_MASS_LOAD_DATA(PELAGIC_BOX_MODEL_DATA % MASS_LOADS(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % MASS_LOADS)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_MASS_WITHDRAWALS
            call DEALLOC_MASS_WITHDRAWAL_DATA(PELAGIC_BOX_MODEL_DATA % MASS_WITHDRAWALS(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % MASS_WITHDRAWALS)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
            call DEALLOC_MEASURED_VALUE(PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES(i))
        end do

        deallocate(PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES)
        deallocate(PELAGIC_BOX_MODEL_DATA % MEASUREMENT_ERRORS)

        deallocate(PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS)
        deallocate(PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_PROCESS_RATE_OUTPUTS)
        deallocate(PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_MASS_BALANCE_OUTPUTS)

        if (allocated(PELAGIC_BOX_MODEL_DATA % SEDIMENT_FLUX_TS_NOS)) then
            deallocate(PELAGIC_BOX_MODEL_DATA % SEDIMENT_FLUX_TS_NOS)
        end if

        if (allocated(PELAGIC_BOX_MODEL_DATA % SEDIMENT_FLUX_TS_VAR_NOS)) then
            deallocate(PELAGIC_BOX_MODEL_DATA % SEDIMENT_FLUX_TS_VAR_NOS)
        end if

        deallocate(PELAGIC_BOX_MODEL_DATA % ADVECTION_ON)
        deallocate(PELAGIC_BOX_MODEL_DATA % DIFFUSION_ON)
        deallocate(PELAGIC_BOX_MODEL_DATA % SETTLING_ON )

        deallocate(PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES           )
        deallocate(PELAGIC_BOX_MODEL_DATA % STATE_VAR_INITIAL_CONDITION_TYPES)

        PELAGIC_BOX_MODEL_DATA % MODEL_CONSTANTS              = 0
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES            = 0
        PELAGIC_BOX_MODEL_DATA % NUM_BASINS                   = 0
        PELAGIC_BOX_MODEL_DATA % NUM_BATHYMETRIES             = 0
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_INIT_COND_SETS   = 0
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_ADVECTIVE_LINKS  = 0
        PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_DISPERSIVE_LINKS = 0
        PELAGIC_BOX_MODEL_DATA % NUM_FLOW_TS                  = 0
        PELAGIC_BOX_MODEL_DATA % NUM_MIXING_TS                = 0
        PELAGIC_BOX_MODEL_DATA % NUM_SETTLING_VELOCITIES      = 0
        PELAGIC_BOX_MODEL_DATA % NUM_OPEN_BOUNDARIES          = 0
        PELAGIC_BOX_MODEL_DATA % NUM_MASS_LOADS               = 0
        PELAGIC_BOX_MODEL_DATA % NUM_MASS_WITHDRAWALS         = 0
        PELAGIC_BOX_MODEL_DATA % NUM_FORCING_TS               = 0
    end subroutine DEALLOC_PELAGIC_BOX_MODEL_DATA


    subroutine READ_PELAGIC_BOX_MODEL_INPUTS(PELAGIC_BOX_MODEL_DATA, IN_FILE, OUT_FILE)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: IN_FILE
        integer, intent(in) :: OUT_FILE

        integer :: i
        integer :: j

        character(len = 2048) :: FILE_NAME

        integer :: NUM_PELAGIC_STATE_VARS
        integer :: NUM_MODEL_CONSTANTS
        integer :: NUM_BASINS
        integer :: NUM_BATHYMETRIES
        integer :: NUM_PELAGIC_BOXES
        integer :: NUM_PELAGIC_INIT_COND_SETS
        integer :: NUM_PELAGIC_ADVECTIVE_LINKS
        integer :: NUM_PELAGIC_DISPERSIVE_LINKS
        integer :: NUM_FLOW_TS
        integer :: NUM_MIXING_TS
        integer :: NUM_SETTLING_VELOCITIES
        integer :: NUM_OPEN_BOUNDARIES
        integer :: NUM_FORCING_TS

        integer :: NUM_MASS_LOADS
        integer :: NUM_MASS_WITHDRAWALS
        integer :: NUM_MASS_LOADS_OF_BOX
        integer :: NUM_MASS_WITHDRAWALS_OF_BOX

        integer :: BASIN_NO
        integer :: BATHYMETRY_NO
        integer :: NUM_BOXES
        integer :: BOX_NO
        integer :: PELAGIC_INIT_COND_SET_NO
        integer :: ADVECTIVE_LINK_NO
        integer :: UPSTREAM_BOX_NO
        integer :: DOWNSTREAM_BOX_NO
        integer :: FLOW_TS_NO
        integer :: FLOW_TS_VAR_NO

        integer :: DISPERSIVE_LINK_NO
        integer :: FIRST_BOX_NO
        integer :: SECOND_BOX_NO
        integer :: MIXING_TS_NO

        real(kind = SELECTED_REAL_KIND(15,307)) :: MIXING_LENGTH
        real(kind = SELECTED_REAL_KIND(15,307)) :: VOLUME
        real(kind = SELECTED_REAL_KIND(15,307)) :: SURFACE_ELEVATION
        real(kind = SELECTED_REAL_KIND(15,307)) :: BOTTOM_ELEVATION
        real(kind = SELECTED_REAL_KIND(15,307)) :: DISSOLVED_FRACTION
        real(kind = SELECTED_REAL_KIND(15,307)) :: DEPOSITED_FRACTION
        integer :: CHLA_SUPRESSION_OF_SETTLING

        integer :: SETTLING_VELOCITY_NO

        integer :: OPEN_BOUNDARY_NO
        integer :: MASS_LOAD_NO
        integer :: MASS_WITHDRAWAL_NO
        integer :: PELAGIC_STATE_VAR_NO
        integer :: FORCING_TS_NO
        integer :: FORCING_TS_VAR_NO

        integer :: STATE_VAR_NO
        character(len = 30) :: STATE_VAR_NAME
        integer :: NUM_BOXES_WITH_MEASURED_VALS
        integer :: AUX_INT

        integer :: PRODUCE_PEL_STATE_VAR_OUTPUT
        integer :: PRODUCE_PEL_PROCESS_RATE_OUTPUT
        integer :: PRODUCE_PEL_MASS_BALANCE_OUTPUT

        character(len = 20) :: FORMAT_STRING
        character(len = 5)  :: INTEGER_STRING
        character(len = 5)  :: BOX_NO_STRING

        integer          :: AUX_INTEGER_1, AUX_INTEGER_2, AUX_INTEGER_3, AUX_INTEGER_4
        real(kind = DBL) :: AUX_DOUBLE_1 , AUX_DOUBLE_2 , AUX_DOUBLE_3 , AUX_DOUBLE_4

        real(kind = DBL) :: TIME

        !READ DESCRIPTION LINES
        do i = 1, 5
            read(IN_FILE, *)
        end do

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_PELAGIC_STATE_VARS
        write(*,*) NUM_PELAGIC_STATE_VARS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_MODEL_CONSTANTS
        write(*,*) NUM_MODEL_CONSTANTS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_BASINS
        write(*,*) NUM_BASINS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_BATHYMETRIES
        write(*,*) NUM_BATHYMETRIES

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_PELAGIC_BOXES
        write(*,*) NUM_PELAGIC_BOXES

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_PELAGIC_INIT_COND_SETS
        write(*,*) NUM_PELAGIC_INIT_COND_SETS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_PELAGIC_ADVECTIVE_LINKS
        write(*,*) NUM_PELAGIC_ADVECTIVE_LINKS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_PELAGIC_DISPERSIVE_LINKS
        write(*,*) NUM_PELAGIC_DISPERSIVE_LINKS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_FLOW_TS
        write(*,*) NUM_FLOW_TS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_MIXING_TS
        write(*,*) NUM_MIXING_TS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_SETTLING_VELOCITIES
        write(*,*) NUM_SETTLING_VELOCITIES

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_OPEN_BOUNDARIES
        write(*,*) NUM_OPEN_BOUNDARIES

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_MASS_LOADS
        write(*,*) NUM_MASS_LOADS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_MASS_WITHDRAWALS
        write(*,*) NUM_MASS_WITHDRAWALS

        read(IN_FILE, *)
        read(IN_FILE, *) NUM_FORCING_TS
        write(*,*) NUM_FORCING_TS

        call ALLOC_PELAGIC_BOX_MODEL_DATA(PELAGIC_BOX_MODEL_DATA,  &
                    NUM_PELAGIC_STATE_VARS      , NUM_MODEL_CONSTANTS         , &
                    NUM_BASINS                  , NUM_BATHYMETRIES            , &
                    NUM_PELAGIC_BOXES           , NUM_PELAGIC_INIT_COND_SETS  , &
                    NUM_PELAGIC_ADVECTIVE_LINKS , NUM_PELAGIC_DISPERSIVE_LINKS, &
                    NUM_SETTLING_VELOCITIES     , NUM_OPEN_BOUNDARIES         , &
                    NUM_FLOW_TS                 , NUM_MIXING_TS               , &
                    NUM_MASS_LOADS              , NUM_MASS_WITHDRAWALS        , &
                    NUM_FORCING_TS)

        read(IN_FILE, *)
        read(IN_FILE, *) FILE_NAME

        open(unit   = IN_FILE + 1, &
             file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // trim(adjustl(FILE_NAME)), &
             status = 'OLD')

        read(IN_FILE + 1, *)
        read(IN_FILE + 1, *) PELAGIC_BOX_MODEL_DATA % ZOOPLANKTON_OPTION

        read(IN_FILE + 1, *)
        read(IN_FILE + 1, *) PELAGIC_BOX_MODEL_DATA % ADVANCED_REDOX_SIMULATION

        if (PELAGIC_BOX_MODEL_DATA % ADVANCED_REDOX_SIMULATION > 0) then
            write(*,*) 'Advanced redox simulation is on.'
        else
            write(*,*) 'Advanced redox simulation is off.'
        end if

        read(IN_FILE + 1, *)
        read(IN_FILE + 1, *) PELAGIC_BOX_MODEL_DATA % LIGHT_EXTINCTION_OPTION

        if (PELAGIC_BOX_MODEL_DATA % LIGHT_EXTINCTION_OPTION > 0) then
            write(*,*) 'Light extinction modelling option as given in "KD".'
        else
            write(*,*) 'Standard light extinction sub model'
        end if

        read(IN_FILE + 1, *)
        read(IN_FILE + 1, *) PELAGIC_BOX_MODEL_DATA % CYANO_BOUYANT_STATE_SIMULATION

        if (PELAGIC_BOX_MODEL_DATA % CYANO_BOUYANT_STATE_SIMULATION > 0) then
            write(*,*) 'Bouyant cyanobacteria option is on.'
        else
            write(*,*) 'Bouyant cyanobacteria option is off. ' // &
                       'Nostocales (if simulated) are allways considered bouyant.'
        end if

        read(IN_FILE + 1, *)
        read(IN_FILE + 1, *) PELAGIC_BOX_MODEL_DATA % CONSIDER_NON_OBLIGATORY_FIXERS

        if (PELAGIC_BOX_MODEL_DATA % CONSIDER_NON_OBLIGATORY_FIXERS > 0) then
            write(*,*) 'Non-abligatory fixers will be modelled'
        else
            write(*,*) 'Non-abligatory fixers will not be modelled.'
        end if

        read(IN_FILE + 1, *)
        read(IN_FILE + 1, *) PELAGIC_BOX_MODEL_DATA % CONSIDER_NOSTOCALES

        if (PELAGIC_BOX_MODEL_DATA % CONSIDER_NOSTOCALES > 0) then
            write(*,*) 'Heterocyst forming fixers and akinetes will be modelled.'
        else
            write(*,*) 'Heterocyst forming fixers and akinetes will not be modelled.'
        end if

        close(IN_FILE + 1)

        ! READ THE MODEL OUTPUT OPTION
        read(IN_FILE, *)
        read(IN_FILE, *) FILE_NAME

        open(unit   = IN_FILE + 1, &
             file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // trim(adjustl(FILE_NAME)), &
             status = 'OLD')

        read(IN_FILE + 1, *)

        do i = 1,  NUM_PELAGIC_BOXES
            read(IN_FILE + 1, *) &
                BOX_NO, PRODUCE_PEL_STATE_VAR_OUTPUT, PRODUCE_PEL_PROCESS_RATE_OUTPUT, &
                PRODUCE_PEL_MASS_BALANCE_OUTPUT

            PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS   (BOX_NO) = &
                PRODUCE_PEL_STATE_VAR_OUTPUT

            PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_PROCESS_RATE_OUTPUTS(BOX_NO) = &
                PRODUCE_PEL_PROCESS_RATE_OUTPUT

            PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_MASS_BALANCE_OUTPUTS(BOX_NO) = &
                PRODUCE_PEL_MASS_BALANCE_OUTPUT
        end do

        close(IN_FILE + 1)

        read(IN_FILE, *)
        read(IN_FILE, *) PEL_PROCESS_RATE_OUTPUT_OPTION

        if ((PEL_PROCESS_RATE_OUTPUT_OPTION < 1).or. &
            (PEL_PROCESS_RATE_OUTPUT_OPTION > 2)) then
            write(*,*) 'Wrong option was entered for the pelagic process output option'
            write(*,*) 'Enter 1 for the volume based output - g/m^3/day'
            write(*,*) 'Enter 2 for the area   based output - g/m^2/day'
            stop
        else
            if (PEL_PROCESS_RATE_OUTPUT_OPTION.eq.1) then
                write(*,*) 'Process rates will be written in g/m^3/day'
            else
                write(*,*) 'Process rates will be written in g/m^2/day'
            end if
        end if

        !READ PELAGIC STATE VAR INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do  i = 1, NUM_PELAGIC_STATE_VARS
            read(IN_FILE, *) &
                 STATE_VAR_NO, STATE_VAR_NAME, NUM_BOXES_WITH_MEASURED_VALS, &
                 AUX_DOUBLE_1, AUX_INTEGER_1 , AUX_INTEGER_2

            PELAGIC_BOX_MODEL_DATA % PELAGIC_STATE_VAR_NAMES(STATE_VAR_NO) = &
                trim(adjustl(STATE_VAR_NAME))

            PELAGIC_BOX_MODEL_DATA % MEASUREMENT_ERRORS               &
                (STATE_VAR_NO) = AUX_DOUBLE_1

            PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES           &
                (STATE_VAR_NO) = AUX_INTEGER_1

            PELAGIC_BOX_MODEL_DATA % STATE_VAR_INITIAL_CONDITION_TYPES&
                (STATE_VAR_NO) = AUX_INTEGER_2

            call ALLOC_MEASURED_VALUE &
                 (PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES(STATE_VAR_NO), &
                  NUM_BOXES_WITH_MEASURED_VALS)

        end do

        !READ PELAGIC MODEL CONSTANTS INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_BOXES
            read(IN_FILE, *) AUX_INT, FILE_NAME

            if (USE_PELAGIC_CONSTANTS_FILE_NAME > 0) then
                open(unit   = IN_FILE + 1, &
                     file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // &
                              trim(adjustl(PELAGIC_CONSTANTS_FILE_NAME)), &
                     status = 'OLD')
            else
                open(unit   = IN_FILE + 1, &
                     file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // trim(adjustl(FILE_NAME)), &
                     status = 'OLD')
            end if

            call READ_MODEL_CONSTANTS &
                 (PELAGIC_BOX_MODEL_DATA % MODEL_CONSTANTS(:, i), IN_FILE + 1)

            close(IN_FILE + 1)
        end do

        ! Read the extra model constant information
        read(IN_FILE, *)
        read(IN_FILE, *) FILE_NAME

        open(unit   = IN_FILE + 1, &
             file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // trim(adjustl(FILE_NAME)), &
             status = 'OLD')

        !  Read the user entered fraction of available DON for cyanobacteria uptake
        read(IN_FILE + 1, *)
        read(IN_FILE + 1, fmt = *) USER_ENTERED_frac_avail_DON

        !  Read the user entered fraction of available DON for cyanobacteria uptake
        read(IN_FILE + 1, *)
        read(IN_FILE + 1, fmt = *) USER_ENTERED_K_B_E

        close(IN_FILE + 1)
        ! End of read the extra model constant information

        !READ PELAGIC BOX INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_BOXES
            read(IN_FILE, *) &
                BOX_NO, NUM_MASS_LOADS_OF_BOX, NUM_MASS_WITHDRAWALS_OF_BOX

            call ALLOC_TOPLOGY_AND_LOADS &
                 (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO), &
                  NUM_MASS_LOADS_OF_BOX, NUM_MASS_WITHDRAWALS_OF_BOX)
        end do

        !READ PELAGIC BASIN INFORMATION
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_BASINS

            read(IN_FILE, *)
            read(IN_FILE, *) BASIN_NO, NUM_BOXES, BATHYMETRY_NO

            call ALLOCATE_BASIN_DATA &
                     (PELAGIC_BOX_MODEL_DATA % BASINS(BASIN_NO), NUM_BOXES)

            read(IN_FILE, *)

            do j = 1, NUM_BOXES
                read(IN_FILE, *) BOX_NO

                PELAGIC_BOX_MODEL_DATA % BASINS(BASIN_NO) % &
                    PELAGIC_BOXES(j)   = BOX_NO

                PELAGIC_BOX_MODEL_DATA % BASINS(BASIN_NO) % &
                    BATHYMETRY_NO      = BATHYMETRY_NO

                PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                    BASIN_NO      = BASIN_NO

                PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                    BATHYMETRY_NO = BATHYMETRY_NO

                PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % BASIN_BOX_NO  = j
            end do

        end do

        !READ PELAGIC BATHYMETRY INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_BATHYMETRIES

            read(IN_FILE, *) BATHYMETRY_NO, FILE_NAME

            open(unit   = IN_FILE + 1, &
                 file   = trim(adjustl(PELAGIC_INPUT_FOLDER)) // trim(adjustl(FILE_NAME)), &
                 status = 'OLD')

            call READ_BATHYMETRY_DATA_FROM_FILE &
                 (PELAGIC_BOX_MODEL_DATA % BATHYMETRIES(BATHYMETRY_NO), IN_FILE + 1)

            PELAGIC_BOX_MODEL_DATA % BATHYMETRIES(BATHYMETRY_NO) % ID_NUM = BATHYMETRY_NO

            close(IN_FILE + 1)
        end do

        !READ INITIAL CONDITIONS
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_BOXES

            read(IN_FILE, *) &
                 BOX_NO, PELAGIC_INIT_COND_SET_NO, SURFACE_ELEVATION, BOTTOM_ELEVATION

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                PELAGIC_INIT_COND_SET_NO = PELAGIC_INIT_COND_SET_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                SURFACE_ELEVATION        = SURFACE_ELEVATION

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                BOTTOM_ELEVATION         = BOTTOM_ELEVATION

            BATHYMETRY_NO = &
                PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % BATHYMETRY_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % VOLUME = &
                CALCULATE_VOLUME_BETWEEN_LEVELS &
                    (PELAGIC_BOX_MODEL_DATA % BATHYMETRIES(BATHYMETRY_NO), &
                     SURFACE_ELEVATION, BOTTOM_ELEVATION)

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % DEPTH = &
                SURFACE_ELEVATION - BOTTOM_ELEVATION

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % SURFACE_AREA = &
                    CALCULATE_SURFACE_AREA &
                       (PELAGIC_BOX_MODEL_DATA % BATHYMETRIES(BATHYMETRY_NO), &
                        SURFACE_ELEVATION)
        end do

        !READ LOAD INFORMATION : MASS LOADS
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_BOXES
            read(IN_FILE, *)

            NUM_MASS_LOADS_OF_BOX = &
                PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % NUM_MASS_LOADS

            do j = 1, NUM_MASS_LOADS_OF_BOX
                read(unit = IN_FILE, fmt = *) &
                    PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % MASS_LOAD_NOS(j)
            end do
        end do

        !READ LOAD INFORMATION : MASS WITHDRAWALS
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_BOXES
            read(IN_FILE, *)

            NUM_MASS_WITHDRAWALS_OF_BOX = &
                PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % NUM_MASS_WITHDRAWALS

            do j = 1, NUM_MASS_WITHDRAWALS_OF_BOX
                read(unit = IN_FILE, fmt = *) &
                     PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % MASS_WITHDRAWAL_NOS(j)
            end do
        end do

        !READ ADVECTIVE LINK INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_ADVECTIVE_LINKS

            read(IN_FILE, *) ADVECTIVE_LINK_NO, UPSTREAM_BOX_NO, &
                             DOWNSTREAM_BOX_NO, FLOW_TS_NO, FLOW_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_ADVECTIVE_LINKS(i) % &
                ADVECTIVE_LINK_NO = ADVECTIVE_LINK_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_ADVECTIVE_LINKS(i) % &
                UPSTREAM_BOX_NO   = UPSTREAM_BOX_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_ADVECTIVE_LINKS(i) % &
                DOWNSTREAM_BOX_NO = DOWNSTREAM_BOX_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_ADVECTIVE_LINKS(i) % &
                FLOW_TS_NO        = FLOW_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_ADVECTIVE_LINKS(i) % &
                FLOW_TS_VAR_NO    = FLOW_TS_VAR_NO
        end do

        !READ DISPERSIVE LINK INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_DISPERSIVE_LINKS

            read(IN_FILE, *) DISPERSIVE_LINK_NO, FIRST_BOX_NO, SECOND_BOX_NO, &
                             MIXING_TS_NO, MIXING_LENGTH

            PELAGIC_BOX_MODEL_DATA % PELAGIC_DISPERSIVE_LINKS(i) % &
                DISPERSIVE_LINK_NO = DISPERSIVE_LINK_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_DISPERSIVE_LINKS(i) % &
                FIRST_BOX_NO       = FIRST_BOX_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_DISPERSIVE_LINKS(i) % &
                SECOND_BOX_NO      = SECOND_BOX_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_DISPERSIVE_LINKS(i) % &
                MIXING_TS_NO       = MIXING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_DISPERSIVE_LINKS(i) % &
                MIXING_LENGTH      = MIXING_LENGTH
        end do

        !READ SETTLING INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_BOXES * NUM_PELAGIC_STATE_VARS

            read(IN_FILE, *) &
                BOX_NO, PELAGIC_STATE_VAR_NO, DISSOLVED_FRACTION, &
                SETTLING_VELOCITY_NO, DEPOSITED_FRACTION, CHLA_SUPRESSION_OF_SETTLING

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                DISSOLVED_FRACTIONS(PELAGIC_STATE_VAR_NO) = DISSOLVED_FRACTION

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                SETTLING_TS_NOS    (PELAGIC_STATE_VAR_NO) = SETTLING_VELOCITY_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                DEPOSITED_FRACTIONS(PELAGIC_STATE_VAR_NO) = DEPOSITED_FRACTION

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                CHLA_SUPRESSION_OF_SETTLING(PELAGIC_STATE_VAR_NO) = &
            CHLA_SUPRESSION_OF_SETTLING
        end do

        !READ OPEN BOUNDARY INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_OPEN_BOUNDARIES * NUM_PELAGIC_STATE_VARS

            read(IN_FILE, *) &
                OPEN_BOUNDARY_NO, PELAGIC_STATE_VAR_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % OPEN_BOUNDARIES(OPEN_BOUNDARY_NO) % &
                FORCING_TS_NOS    (PELAGIC_STATE_VAR_NO) = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % OPEN_BOUNDARIES(OPEN_BOUNDARY_NO) % &
                FORCING_TS_VAR_NOS(PELAGIC_STATE_VAR_NO) = FORCING_TS_VAR_NO
        end do

        !READ MASS LOAD INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_MASS_LOADS * NUM_PELAGIC_STATE_VARS

            read(IN_FILE, *) &
                MASS_LOAD_NO , PELAGIC_STATE_VAR_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % MASS_LOADS(MASS_LOAD_NO) % &
                FORCING_TS_NOS    (PELAGIC_STATE_VAR_NO)     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % MASS_LOADS(MASS_LOAD_NO) % &
                FORCING_TS_VAR_NOS(PELAGIC_STATE_VAR_NO) = FORCING_TS_VAR_NO
        end do

        !READ MASS WITHDRAWAL INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_MASS_WITHDRAWALS * NUM_PELAGIC_STATE_VARS

            read(IN_FILE, *) MASS_WITHDRAWAL_NO, PELAGIC_STATE_VAR_NO, &
                             FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % MASS_WITHDRAWALS(MASS_WITHDRAWAL_NO) % &
                    FORCING_TS_NOS(PELAGIC_STATE_VAR_NO)     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % MASS_WITHDRAWALS(MASS_WITHDRAWAL_NO) % &
                    FORCING_TS_VAR_NOS(PELAGIC_STATE_VAR_NO) = FORCING_TS_VAR_NO
        end do

        !READ TEMPERATURE INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                TEMPERATURE_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                TEMPERATURE_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ SALINITY INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                SALINITY_TS_NO    = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                SALINITY_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ SOLAR RADIATION INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                SOLAR_RADIATION_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                SOLAR_RADIATION_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ FRACTION OF DAY INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                FRACTION_OF_DAY_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                FRACTION_OF_DAY_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ AIR TEMPERATURE INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                AIR_TEMPERATURE_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                AIR_TEMPERATURE_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ WIND SPEED INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                WIND_SPEED_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                WIND_SPEED_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ PRECIPITATION INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                PRECIPITATION_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                PRECIPITATION_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ EVAPORATION INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                EVAPORATION_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                EVAPORATION_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ ICE FRACTION INFORMATION
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            read(IN_FILE, *) BOX_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                ICE_FRACTION_TS_NO     = FORCING_TS_NO

            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(BOX_NO) % &
                ICE_FRACTION_TS_VAR_NO = FORCING_TS_VAR_NO
        end do

        !READ INITIAL CONDITIONS
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_PELAGIC_INIT_COND_SETS
            read(unit = IN_FILE, fmt = *) PELAGIC_INIT_COND_SET_NO, FILE_NAME

            open(unit   = IN_FILE + 1, &
                 file   = trim(adjustl(PELAGIC_INPUT_FOLDER))//trim(adjustl(FILE_NAME)),&
                 status = 'OLD')

            call READ_INITIAL_CONCENTRATIONS &
                 (PELAGIC_BOX_MODEL_DATA % &
                      PELAGIC_INIT_COND_SETS(PELAGIC_INIT_COND_SET_NO), IN_FILE + 1)

            close(IN_FILE + 1)
        end do

        !READ FLOW TIME SERIES
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_FLOW_TS
            read(unit = IN_FILE, fmt = *) FLOW_TS_NO, FILE_NAME

            open(unit   = IN_FILE + 1, &
                 file   = trim(adjustl(PELAGIC_INPUT_FOLDER))//trim(adjustl(FILE_NAME)),&
                 status = 'OLD')

            call INITIALIZE_TIME_SERIE(PELAGIC_BOX_MODEL_DATA % FLOWS(FLOW_TS_NO))

            call READ_TIME_SERIE_FROM_FILE &
                     (PELAGIC_BOX_MODEL_DATA % FLOWS(FLOW_TS_NO), IN_FILE + 1)

            close(IN_FILE + 1)
        end do

        !READ MIXING TIME SERIES
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_MIXING_TS
            read(unit = IN_FILE, fmt = *) MIXING_TS_NO, FILE_NAME

            open(unit   = IN_FILE + 1, &
                 file   = trim(adjustl(PELAGIC_INPUT_FOLDER))//trim(adjustl(FILE_NAME)),&
                 status = 'OLD')

            call INITIALIZE_TIME_SERIE &
                     (PELAGIC_BOX_MODEL_DATA % MIXING_EXCHANGES(MIXING_TS_NO))

            call READ_TIME_SERIE_FROM_FILE &
                     (PELAGIC_BOX_MODEL_DATA % MIXING_EXCHANGES(MIXING_TS_NO), &
                      IN_FILE + 1)

            close(IN_FILE + 1)
        end do

        !READ SETTLING VELOCITY TIME SERIES
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_SETTLING_VELOCITIES
            read(unit = IN_FILE, fmt = *) SETTLING_VELOCITY_NO, FILE_NAME

            open(unit   = IN_FILE + 1, &
                 file   = trim(adjustl(PELAGIC_INPUT_FOLDER))//trim(adjustl(FILE_NAME)),&
                 status = 'OLD')

            call INITIALIZE_TIME_SERIE &
                     (PELAGIC_BOX_MODEL_DATA % SETTLING_VELOCITIES(SETTLING_VELOCITY_NO))

            call READ_TIME_SERIE_FROM_FILE &
                     (PELAGIC_BOX_MODEL_DATA % &
                         SETTLING_VELOCITIES(SETTLING_VELOCITY_NO), IN_FILE + 1)

            close(IN_FILE + 1)
        end do

        !READ FORCING TIME SERIES
        read(IN_FILE, *)
        read(IN_FILE, *)

        do i = 1, NUM_FORCING_TS
            read(unit = IN_FILE, fmt = *) FORCING_TS_NO, FILE_NAME

            open(unit   = IN_FILE + 1, &
                 file   = trim(adjustl(PELAGIC_INPUT_FOLDER))//trim(adjustl(FILE_NAME)),&
                 status = 'OLD')

            call INITIALIZE_TIME_SERIE &
                 (PELAGIC_BOX_MODEL_DATA % FORCING_TS(FORCING_TS_NO))

            call READ_TIME_SERIE_FROM_FILE &
                 (PELAGIC_BOX_MODEL_DATA % FORCING_TS(FORCING_TS_NO), &
                  IN_FILE + 1)

            close(IN_FILE + 1)
        end do

        !READ PELAGIC OUTPUT INFORMARION
        PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_ECOL_OUTPUT = 0

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *)

        read(unit = IN_FILE, fmt = *) &
             PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_ECOL_OUTPUT

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *)

        read(unit = IN_FILE, fmt = *) &
             PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_SAVED_OUTPUT

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *)

        read(unit = IN_FILE, fmt = *) &
             PELAGIC_BOX_MODEL_DATA % CREATE_STATE_VARIABLE_OUTPUT

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *) PRODUCE_COCOA_OUTPUTS

        if (PRODUCE_COCOA_OUTPUTS > 0) then
            read(unit = IN_FILE, fmt = *)
            read(unit = IN_FILE, fmt = *) COCOA_PELAGIC_OUTPUTS_FILENAME
        end if

        !READ PELAGIC EXERGY OUTPUT INFORMATION
        PELAGIC_BOX_MODEL_DATA % CALCULATE_PELAGIC_EXERGY      = 0
        PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_EXERGY_OUTPUTS = 0
        PELAGIC_BOX_MODEL_DATA % START_REPEAT_NO_PEL_EX_OUTS   = 0

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *)

        read(unit = IN_FILE, fmt = *) &
             PELAGIC_BOX_MODEL_DATA % CALCULATE_PELAGIC_EXERGY

        if (PELAGIC_BOX_MODEL_DATA % CALCULATE_PELAGIC_EXERGY > 0) then
            read(unit = IN_FILE, fmt = *)

            read(unit = IN_FILE, fmt = *) &
                 PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_EXERGY_OUTPUTS

            if (PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_EXERGY_OUTPUTS > 0) then
                read(unit = IN_FILE, fmt = *)

                read(unit = IN_FILE, fmt = *) &
                     PELAGIC_BOX_MODEL_DATA % START_REPEAT_NO_PEL_EX_OUTS
            end if

        end if

        !READ THE COST FUNCTION INFORMATION
        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *)

        read(unit = IN_FILE, fmt = *) &
             PELAGIC_BOX_MODEL_DATA % PRODUCE_COST_FUNCTION

        if (PELAGIC_BOX_MODEL_DATA % PRODUCE_COST_FUNCTION > 0) then

            read(unit = IN_FILE, fmt = *)

            do i = 1, NUM_PELAGIC_STATE_VARS
                read(unit = IN_FILE, fmt = *)
                read(unit = IN_FILE, fmt = *)

                do j = 1, PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES(i) % NUM_BOXES

                    read(unit = IN_FILE, fmt = *) &
                         PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES(i) % &
                             BOX_NOS           (j), &
                         PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES(i) % &
                             FORCING_TS_NOS    (j), &
                         PELAGIC_BOX_MODEL_DATA % MEASURED_VALUES(i) % &
                            FORCING_TS_VAR_NOS(j)
                end do
            end do
        end if

        !INITIALIZE THE OUTPUT FILES
        TIME                   = INIT_TIME
        NUM_BOXES              = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
        NUM_PELAGIC_STATE_VARS = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

        if (PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_ECOL_OUTPUT > 0) then

            if (PRODUCE_ONLY_BINARY_PELAGIC_OUTPUT == 0) then

                !Pelagic state variables output files
                if (PELAGIC_BOX_MODEL_DATA % CREATE_STATE_VARIABLE_OUTPUT > 0) then
                    do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

                        open(unit = AUX_OUTPUT_UNIT, &
                             file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                    trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                                 PELAGIC_STATE_VAR_NAMES(i))) &
                                    // '.out', status = 'UNKNOWN')

                        FORMAT_STRING = '(a20, '
                        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES
                        FORMAT_STRING = '(a20, ' // INTEGER_STRING // 'i20)'

                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) &
                              '         TIME_DAYS', (j, j = 1, NUM_BOXES)

                        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES + 1
                        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                               CONCENTRATIONS(i), j = 1, NUM_BOXES)

                        close(AUX_OUTPUT_UNIT)
                    end do
                end if

                !Pelagic derived variables output files
                if (PELAGIC_BOX_MODEL_DATA % CREATE_STATE_VARIABLE_OUTPUT > 0) then
                    do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

                        open(unit = AUX_OUTPUT_UNIT, &
                             file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                    trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                                 PELAGIC_STATE_VAR_NAMES(i))) &
                                    // '.out', status = 'UNKNOWN')

                        FORMAT_STRING = '(a20, '
                        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES
                        FORMAT_STRING = '(a20, ' // INTEGER_STRING // 'i20)'

                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) &
                              '         TIME_DAYS', (j, j = 1, NUM_BOXES)

                        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES + 1
                        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                               CONCENTRATIONS(i), j = 1, NUM_BOXES)

                        close(AUX_OUTPUT_UNIT)
                    end do
                end if

                !Pelagic box files
                do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
                    if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS(i) > 0) then

                        PELAGIC_INIT_COND_SET_NO = &
                            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % PELAGIC_INIT_COND_SET_NO

                        do j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

                            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % CONCENTRATIONS(j) = &
                                PELAGIC_BOX_MODEL_DATA % &
                                    PELAGIC_INIT_COND_SETS(PELAGIC_INIT_COND_SET_NO) % &
                                    PELAGIC_CONCENTRATIONS(j)
                        end do

                        write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                        open(unit   = AUX_OUTPUT_UNIT, &
                             file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                      trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // '.out')), &
                             status = 'UNKNOWN')

                        FORMAT_STRING = '(a20, '
                        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_PELAGIC_STATE_VARS
                        FORMAT_STRING = '(a20, ' // INTEGER_STRING // 'a20)'

                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) '         TIME_DAYS', &
                             (trim(adjustl(PELAGIC_BOX_MODEL_DATA % PELAGIC_STATE_VAR_NAMES(j))), &
                                           j = 1, NUM_PELAGIC_STATE_VARS)

                        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_PELAGIC_STATE_VARS + 1
                        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) INIT_TIME, &
                              (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % CONCENTRATIONS(j), &
                              j = 1, NUM_PELAGIC_STATE_VARS)

                        close(AUX_OUTPUT_UNIT)

                        write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                        open(unit = AUX_OUTPUT_UNIT, &
                             file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                    trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // &
                                                 '.mtrx')), status = 'UNKNOWN')

                        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_PELAGIC_STATE_VARS + 1
                        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) INIT_TIME, &
                              (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % CONCENTRATIONS(j), &
                               j = 1, NUM_PELAGIC_STATE_VARS)

                        close(AUX_OUTPUT_UNIT)
                    end if
                end do
            end if
        end if

        FORMAT_STRING = '(a20, '
        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES
        FORMAT_STRING = '(a20, ' // INTEGER_STRING // 'i20)'

        !Pelagic saved output files
        if (PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_SAVED_OUTPUT == 1) then
            if (PRODUCE_ONLY_BINARY_PELAGIC_OUTPUT == 0) then
                do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_SAVED_OUTPUTS

                    open(unit = AUX_OUTPUT_UNIT, &
                         file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                                 SAVED_OUTPUT_NAMES(i))) // '.out', &
                         status = 'UNKNOWN')

                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) &
                          '         TIME_DAYS', (j, j = 1, NUM_BOXES)

                    close(AUX_OUTPUT_UNIT)
                end do
            end if
        end if

        !Pelagic exergy output files
        if (PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_EXERGY_OUTPUTS == 1) then
            if (PRODUCE_ONLY_BINARY_PELAGIC_OUTPUT == 0) then
                do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_EXERGY_COMPONENTS

                    open(unit  = AUX_OUTPUT_UNIT, &
                         file  = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                 trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                                 EXERGY_COMPONENT_NAMES(i)))// '.out',  &
                         status = 'UNKNOWN')

                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) &
                          '         TIME_DAYS', (j, j = 1, NUM_BOXES)

                    close(AUX_OUTPUT_UNIT)
                end do
            end if
        end if
    end subroutine READ_PELAGIC_BOX_MODEL_INPUTS


    subroutine READ_BOTTOM_SEDIMENTS_FLUXES_INPUTS&
               (PELAGIC_BOX_MODEL_DATA, IN_FILE, FLUX_SET_NO)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: IN_FILE
        integer, intent(in) :: FLUX_SET_NO

        integer :: BOX_NO
        integer :: STATE_VAR_NO
        integer :: FORCING_TS_NO
        integer :: FORCING_TS_VAR_NO
        integer :: i
        integer :: j

        read(unit = IN_FILE, fmt = *)
        read(unit = IN_FILE, fmt = *)

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            do j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
                read(unit = IN_FILE, fmt = *) &
                    BOX_NO, STATE_VAR_NO, FORCING_TS_NO, FORCING_TS_VAR_NO

                PELAGIC_BOX_MODEL_DATA % SEDIMENT_FLUX_TS_NOS     &
                    (FLUX_SET_NO, BOX_NO, STATE_VAR_NO) = FORCING_TS_NO

                PELAGIC_BOX_MODEL_DATA % SEDIMENT_FLUX_TS_VAR_NOS &
                    (FLUX_SET_NO, BOX_NO, STATE_VAR_NO) = FORCING_TS_VAR_NO
            end do
        end do

    end subroutine READ_BOTTOM_SEDIMENTS_FLUXES_INPUTS


    subroutine WRITE_PELAGIC_OUTPUT &
               (PELAGIC_BOX_MODEL_DATA, TIME, AUX_OUTPUT_UNIT, &
                WRITE_SAVED_OUTPUT, WRITE_PELAGIC_EXERGY_OUTPUT, &
                MODEL_AT_INIT)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA
        real(kind = DBL), intent(in) :: TIME
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: WRITE_SAVED_OUTPUT
        integer, intent(in) :: WRITE_PELAGIC_EXERGY_OUTPUT
        integer, intent(in) :: MODEL_AT_INIT

        character(len = 20) :: FORMAT_STRING
        character(len = 5)  :: INTEGER_STRING
        character(len = 5)  :: BOX_NO_STRING

        integer i
        integer j
        integer NUM_BOXES

        double precision, dimension(nkn)    :: MEAN_BOX_DEPTHS
        double precision, dimension(nstate) :: WRITTEN_CONC

        do i = 1, nkn

            MEAN_BOX_DEPTHS(i) = &
                 PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % VOLUME       /   &
                 PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % SURFACE_AREA
        end do

        BOX_NO_STRING = '     '

        NUM_BOXES = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES + 1
        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

        if (PELAGIC_BOX_MODEL_DATA % CREATE_STATE_VARIABLE_OUTPUT > 0) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                         PELAGIC_STATE_VAR_NAMES(i))) &
                            // '.out', status = 'UNKNOWN', position = 'APPEND')

                ! If the desired output is in g/m^2 than convert the out concentration
                ! from g/m^3 to g/m^2
                if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(i) > 1) then
                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) &
                          TIME, &
                          ((PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                               CONCENTRATIONS(i) * MEAN_BOX_DEPTHS(j)), j = 1, NUM_BOXES)
                else
                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) &
                          TIME, &
                          (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                               CONCENTRATIONS(i), j = 1, NUM_BOXES)
                end if

                close(AUX_OUTPUT_UNIT)
            end do
        end if

        ! -------------------------------------------------------------------------------
        ! CREATE PELAGIC BOX OUTPUT FILES
        ! -------------------------------------------------------------------------------
        do i = 1, NUM_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                if (MODEL_AT_INIT > 0) then
                    open(unit   = AUX_OUTPUT_UNIT, &
                         file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                  trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // '.out')), &
                         status = 'UNKNOWN')

                    FORMAT_STRING = '(a20, '

                    write(unit = INTEGER_STRING, fmt = '(i5)') &
                          PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

                    FORMAT_STRING = '(a20, ' // INTEGER_STRING // 'a20)'

                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) '         TIME_DAYS', &
                             (trim(adjustl(PELAGIC_BOX_MODEL_DATA % PELAGIC_STATE_VAR_NAMES(j))), &
                             j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS)

                    close(AUX_OUTPUT_UNIT)
                end if

                open(unit   = AUX_OUTPUT_UNIT, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // '.out')), &
                     status = 'OLD', position = 'APPEND')

                write(unit = INTEGER_STRING, fmt = '(i5)') &
                PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS + 1

                FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'


                ! If the desired output is in g/m^2 than convert the out concentration
                ! from g/m^3 to g/m^2
                do j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
                   if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(j) > 1) then
                       WRITTEN_CONC(j) = &
                           PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % CONCENTRATIONS(j) * &
                           MEAN_BOX_DEPTHS(i)
                   else
                       WRITTEN_CONC(j) = &
                           PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % CONCENTRATIONS(j)
                   end if
                end do

                write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                          ((WRITTEN_CONC(j)), &
                          j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS)

                close(AUX_OUTPUT_UNIT)

                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                             trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // &
                                          '.mtrx')), status = 'UNKNOWN')

                write(unit = INTEGER_STRING, fmt = '(i5)') &
                      PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS + 1

                FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

                write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                          ((WRITTEN_CONC(j)), &
                          j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS)

                close(AUX_OUTPUT_UNIT)
            end if
        end do
        ! -------------------------------------------------------------------------------
        ! END OF CREATE PELAGIC BOX OUTPUT FILES
        ! -------------------------------------------------------------------------------


        if (WRITE_SAVED_OUTPUT == 1) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_SAVED_OUTPUTS

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                         SAVED_OUTPUT_NAMES(i))) &
                            // '.out', status = 'UNKNOWN', position = 'APPEND')

                write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                      (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                           SAVED_OUTPUTS(i), j = 1, NUM_BOXES)

                close(AUX_OUTPUT_UNIT)
            end do
        end if

        if (WRITE_PELAGIC_EXERGY_OUTPUT == 1) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_EXERGY_COMPONENTS

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                         EXERGY_COMPONENT_NAMES(i))) &
                            // '.out', status = 'UNKNOWN', position = 'APPEND')

                write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                      (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                           EXERGY_COMPONENTS(i), j = 1, NUM_BOXES)

                close(AUX_OUTPUT_UNIT)
            end do
        end if
    end subroutine WRITE_PELAGIC_OUTPUT


    subroutine WRITE_PELAGIC_BINARY_OUTPUT &
               (PELAGIC_BOX_MODEL_DATA, TIME, AUX_OUTPUT_UNIT, &
                WRITE_SAVED_OUTPUT, WRITE_PELAGIC_EXERGY_OUTPUT, &
                MODEL_AT_INIT)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA
        real(kind = DBL), intent(in) :: TIME
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: WRITE_SAVED_OUTPUT
        integer, intent(in) :: WRITE_PELAGIC_EXERGY_OUTPUT
        integer, intent(in) :: MODEL_AT_INIT

        integer i
        integer j
        integer NUM_BOXES
        integer NUM_PROCESS_RATES

        double precision, dimension(nstate) :: WRITTEN_CONC
        double precision, dimension(nkn, (nstate * NDIAGVAR)) :: WRITTEN_PROCESS_RATE
        double precision, dimension(nkn)    :: MEAN_BOX_DEPTHS
        character(len = 5)  :: BOX_NO_STRING

        NUM_PROCESS_RATES = nstate * NDIAGVAR

        do i = 1, nkn

            MEAN_BOX_DEPTHS(i) = &
                 PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % VOLUME       /   &
                 PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % SURFACE_AREA
        end do

        NUM_BOXES = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES


        ! -------------------------------------------------------------------------------
        ! CREATE PELAGIC BOX OUTPUT FILES
        ! -------------------------------------------------------------------------------
        do i = 1, NUM_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                ! If the desired output is in g/m^2 than convert the out concentration
                ! from g/m^3 to g/m^2
                do j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
                    if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(j) > 1) then

                        WRITTEN_CONC(j) = &
                            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % &
                            CONCENTRATIONS(j) * MEAN_BOX_DEPTHS(i)
                    else
                        WRITTEN_CONC(j) = &
                            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % CONCENTRATIONS(j)
                    end if
                end do

                if (MODEL_AT_INIT > 0) then
                    open(unit   = AUX_OUTPUT_UNIT, &
                         file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                  trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) &
                                  // "_" // &
                                  trim(adjustl('PELAGIC_BOX_' // &
                                                BOX_NO_STRING // '.bin')), &
                         status = 'REPLACE', access = 'STREAM')

                    write(unit = AUX_OUTPUT_UNIT)  &
                          TIME, ((WRITTEN_CONC(j)), &
                              j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS)

                    close(AUX_OUTPUT_UNIT)
                else
                    open(unit   = AUX_OUTPUT_UNIT, &
                         file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                  trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME))  &
                                  // "_" // &
                                  trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // &
                                  '.bin')), &
                         status = 'OLD', position = 'APPEND', access = 'STREAM')

                    write(unit = AUX_OUTPUT_UNIT)  &
                          TIME, ((WRITTEN_CONC(j)), &
                              j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS)

                    close(AUX_OUTPUT_UNIT)
                end if
            end if
        end do
        ! -------------------------------------------------------------------------------
        ! END OF CREATE PELAGIC BOX OUTPUT FILES
        ! -------------------------------------------------------------------------------


        ! -------------------------------------------------------------------------------
        ! CREATE PELAGIC BOX PROCESS RATES OUTPUT FILES
        ! -------------------------------------------------------------------------------
        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_PROCESS_RATE_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                if (MODEL_AT_INIT > 0) then
                    open(unit   = AUX_OUTPUT_UNIT, &
                         file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                                  trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) &
                                  // "_" // &
                                  trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // &
                                               '_PROCESS_RATES.bin')), &
                         status ='REPLACE')

                    close(AUX_OUTPUT_UNIT)
                end if

                open(unit   = AUX_OUTPUT_UNIT, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) &
                              // "_" // &
                              trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // &
                                           '_PROCESS_RATES.bin')), &
                     status = 'OLD', position = 'APPEND', access = 'STREAM')

                ! If the desired output is in g/m^2 than convert the out concentration
                ! from g/m^3 to g/m^2
                do j = 1, NUM_PROCESS_RATES
                    if (PEL_PROCESS_RATE_OUTPUT_OPTION == 1) then

                        WRITTEN_PROCESS_RATE(i, j) = &
                            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % PROCESS_RATES(j)
                    else
                        WRITTEN_PROCESS_RATE(i, j) = &
                            PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % PROCESS_RATES(j) * &
                            (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % VOLUME / &
                             PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % SURFACE_AREA)
                    end if
                end do

                write(unit = AUX_OUTPUT_UNIT)  &
                      TIME, (WRITTEN_PROCESS_RATE(i, j), j = 1, NUM_PROCESS_RATES)

                close(AUX_OUTPUT_UNIT)
            end if
        end do
        ! -------------------------------------------------------------------------------
        ! END OF CREATE PELAGIC BOX PROCESS RATES OUTPUT FILES
        ! -------------------------------------------------------------------------------


        if (PELAGIC_BOX_MODEL_DATA % CREATE_STATE_VARIABLE_OUTPUT > 0) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" &
				            // trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                            PELAGIC_STATE_VAR_NAMES(i))) // '.bin', &
						    status = 'UNKNOWN', position = 'APPEND', access = 'STREAM')

                ! If the desired output is in g/m^2 than convert the out concentration
                ! from g/m^3 to g/m^2
                if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(i) > 1) then

                       write(unit = AUX_OUTPUT_UNIT ) TIME, &
                          ((PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                              CONCENTRATIONS(i) *  MEAN_BOX_DEPTHS(j)), j = 1, NUM_BOXES)
                else
                       write(unit = AUX_OUTPUT_UNIT ) TIME, &
                          (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                              CONCENTRATIONS(i), j = 1, NUM_BOXES)
                end if

                close(AUX_OUTPUT_UNIT)
            end do
        end if

        if (WRITE_SAVED_OUTPUT == 1) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_SAVED_OUTPUTS

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" &
				            // trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                            SAVED_OUTPUT_NAMES(i))) // '.bin', &
							status = 'UNKNOWN', position = 'APPEND', access = 'STREAM')

                write(unit = AUX_OUTPUT_UNIT) TIME, &
                      (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                           SAVED_OUTPUTS(i), j = 1, NUM_BOXES)

                close(AUX_OUTPUT_UNIT)
            end do
        end if


        if (WRITE_PELAGIC_EXERGY_OUTPUT == 1) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_EXERGY_COMPONENTS

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" &
				            // trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                             EXERGY_COMPONENT_NAMES(i))) // '.bin', &
                            status = 'UNKNOWN', position = 'APPEND', access = 'STREAM')

                write(unit = AUX_OUTPUT_UNIT) TIME, &
                      (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(j) % &
                           EXERGY_COMPONENTS(i), j = 1, NUM_BOXES)

                close(AUX_OUTPUT_UNIT)
            end do
        end if
    end subroutine WRITE_PELAGIC_BINARY_OUTPUT


    subroutine WRITE_PELAGIC_MEM_OUTPUT &
                  (PELAGIC_BOX_MODEL_DATA, AUX_OUTPUT_UNIT, NUM_LINES)

        type(PELAGIC_BOX_MODEL_DS), intent(in) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: NUM_LINES

        real(kind = DBL) :: TIME

        character(len = 20) :: FORMAT_STRING
        character(len = 20) :: FORMAT_STRING_PROCESS
        character(len = 5)  :: INTEGER_STRING
        character(len = 5)  :: BOX_NO_STRING

        integer i
        integer j
        integer LINE_NO
        integer NUM_BOXES
        integer NUM_STATE_VARS
        integer NUM_PROCESS_RATES

        NUM_BOXES = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES + 1
        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

        NUM_PROCESS_RATES = PELAGIC_BOX_MODEL_DATA % NUM_PROCESS_RATES
        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_PROCESS_RATES + 1
        FORMAT_STRING_PROCESS = '(' // INTEGER_STRING // 'f20.6)'

        if (PELAGIC_BOX_MODEL_DATA % CREATE_STATE_VARIABLE_OUTPUT > 0) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                         PELAGIC_STATE_VAR_NAMES(i))) &
                            // '.out', status = 'UNKNOWN', position = 'APPEND')

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, 0, i)

                    ! If the desired output is in g/m^2 than convert the out concentration
                    ! from g/m^3 to g/m^2
                    if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(i) > 1) then
                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % &
                                  AREA_BASED_ECOL_RESULTS(LINE_NO, j, i), &
                              j = 1, NUM_BOXES)
                    else
                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, j, i), &
                                   j = 1, NUM_BOXES)
                    end if
                end do

                close(AUX_OUTPUT_UNIT)
            end do
        end if

        NUM_STATE_VARS = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_STATE_VARS + 1
        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // '.out')), &
                     status = 'UNKNOWN', position = 'APPEND')

                do j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

                    ! If the desired output is in g/m^2 than convert the out concentration
                    ! from g/m^3 to g/m^2
                    if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(j) > 1) then
                        PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(:, i, j) = &
                            PELAGIC_BOX_MODEL_DATA % AREA_BASED_ECOL_RESULTS(:, i, j)
                    end if
                end do

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, 0, 1)

                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                          (PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, i, j), &
                           j = 1, NUM_STATE_VARS)
                end do

                close(AUX_OUTPUT_UNIT)
            end if
        end do


        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                open(unit   = AUX_OUTPUT_UNIT, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // '.mtrx')), &
                     status = 'UNKNOWN', position = 'APPEND')

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, 0, 1)

                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                          (PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, i, j), &
                           j = 1, NUM_STATE_VARS)
                end do

                close(AUX_OUTPUT_UNIT)
            end if
        end do


        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_PROCESS_RATE_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                open(unit   = AUX_OUTPUT_UNIT, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // &
                                            '_PROCESS_RATES.mtrx')), status = 'UNKNOWN')

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % PROCESS_RATE_RESULTS(LINE_NO, 0, 1)

                    if (PEL_PROCESS_RATE_OUTPUT_OPTION == 1) then
                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING_PROCESS) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % PROCESS_RATE_RESULTS(LINE_NO, i, j), &
                               j = 1, NUM_PROCESS_RATES)
                    else
                        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING_PROCESS) TIME, &
                              ((PELAGIC_BOX_MODEL_DATA % PROCESS_RATE_RESULTS(LINE_NO, i, j) * &
                                ((PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % VOLUME) / &
                                 (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % SURFACE_AREA))), &
                               j = 1, NUM_PROCESS_RATES)
                    end if
                end do

                close(AUX_OUTPUT_UNIT)
            end if
        end do
    end subroutine WRITE_PELAGIC_MEM_OUTPUT


    subroutine WRITE_PELAGIC_MEM_BINARY_OUTPUT &
                  (PELAGIC_BOX_MODEL_DATA, AUX_OUTPUT_UNIT, NUM_LINES)

        type(PELAGIC_BOX_MODEL_DS), intent(in) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: NUM_LINES

        real(kind = DBL) :: TIME

        integer i
        integer j
        integer LINE_NO
        integer NUM_BOXES
        integer NUM_STATE_VARS
        integer NUM_PROCESS_RATES

        character(len = 5)  :: BOX_NO_STRING

        NUM_BOXES         = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
        NUM_PROCESS_RATES = PELAGIC_BOX_MODEL_DATA % NUM_PROCESS_RATES
        NUM_STATE_VARS    = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

        if (PELAGIC_BOX_MODEL_DATA % CREATE_STATE_VARIABLE_OUTPUT > 0) then
            do i = 1, NUM_STATE_VARS
                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" // &
				            trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                          PELAGIC_STATE_VAR_NAMES(i))) // '.bin', &
                     status = 'UNKNOWN', access = 'STREAM')

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, 0, i)

                    ! If the desired output is in g/m^2 than convert the out concentration
                    ! from g/m^3 to g/m^2
                    if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(i) > 1) then
                        write(unit = AUX_OUTPUT_UNIT) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % &
                                  AREA_BASED_ECOL_RESULTS(LINE_NO, j, i), &
                              j = 1, NUM_BOXES)
                    else
                        write(unit = AUX_OUTPUT_UNIT) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, j, i), &
                                   j = 1, NUM_BOXES)
                    end if
                end do

                close(AUX_OUTPUT_UNIT)
            end do
        end if

        do i = 1, NUM_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_STATE_VAR_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" // &
                            trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // '.bin')), &
                     status = 'UNKNOWN', access = 'STREAM')

                do j = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS

                    ! If the desired output is in g/m^2 than convert the out concentration
                    ! from g/m^3 to g/m^2
                    if (PELAGIC_BOX_MODEL_DATA % STATE_VAR_OUTPUT_TYPES(j) > 1) then
                        PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(:, i, j) = &
                            PELAGIC_BOX_MODEL_DATA % AREA_BASED_ECOL_RESULTS(:, i, j)
                    end if
                end do

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, 0, 1)

                    write(unit = AUX_OUTPUT_UNIT) TIME, &
                          (PELAGIC_BOX_MODEL_DATA % ECOL_RESULTS(LINE_NO, i, j), &
                           j = 1, NUM_STATE_VARS)
                end do

                close(AUX_OUTPUT_UNIT)
            end if
        end do

        do i = 1, NUM_BOXES
            if (PELAGIC_BOX_MODEL_DATA % PRODUCE_PEL_PROCESS_RATE_OUTPUTS(i) > 0) then
                write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

                open(unit   = AUX_OUTPUT_UNIT, &
                     file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                              trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" // &
              	              trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // &
                             '_PROCESS_RATES.bin')), status = 'UNKNOWN', access = 'STREAM')

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % PROCESS_RATE_RESULTS(LINE_NO, 0, 1)

                    if (PEL_PROCESS_RATE_OUTPUT_OPTION == 1) then
                        write(unit = AUX_OUTPUT_UNIT) TIME, &
                              (PELAGIC_BOX_MODEL_DATA % PROCESS_RATE_RESULTS(LINE_NO, i, j), &
                               j = 1, NUM_PROCESS_RATES)
                    else
                        write(unit = AUX_OUTPUT_UNIT) TIME, &
                              ((PELAGIC_BOX_MODEL_DATA % PROCESS_RATE_RESULTS(LINE_NO, i, j) * &
                                (PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % VOLUME       /   &
                                 PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % SURFACE_AREA)), &
                               j = 1, NUM_PROCESS_RATES)
                    end if
                end do

                close(AUX_OUTPUT_UNIT)
            end if
        end do
    end subroutine WRITE_PELAGIC_MEM_BINARY_OUTPUT


    subroutine WRITE_PELAGIC_MEM_SAVED_OUTPUT &
               (PELAGIC_BOX_MODEL_DATA, AUX_OUTPUT_UNIT, NUM_LINES)

        type(PELAGIC_BOX_MODEL_DS), intent(in) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: NUM_LINES

        real(kind = DBL) :: TIME

        character(len = 20) :: FORMAT_STRING
        character(len = 5)  :: INTEGER_STRING
        character(len = 5)  :: BOX_NO_STRING

        integer i
        integer j
        integer LINE_NO
        integer NUM_BOXES
        integer NUM_SAVED_OUTPUTS

        NUM_BOXES  = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES + 1
        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_SAVED_OUTPUTS

            open(unit   = AUX_OUTPUT_UNIT, &
                 file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                          trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                         SAVED_OUTPUT_NAMES(i))) // '.out', &
                 status = 'UNKNOWN', position = 'APPEND')

            do LINE_NO = 1, NUM_LINES
                TIME = PELAGIC_BOX_MODEL_DATA % SAVED_OUTPUT_RESULTS(LINE_NO, 0, i)

                write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                      (PELAGIC_BOX_MODEL_DATA % SAVED_OUTPUT_RESULTS(LINE_NO, j, i), &
                       j = 1, NUM_BOXES)
            end do

            close(AUX_OUTPUT_UNIT)
        end do

        !NUM_SAVED_OUTPUTS = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_STATE_VARS
        !write(unit = INTEGER_STRING, fmt = '(i5)') NUM_STATE_VARS + 1
        !FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'


        !do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES

        !    write(unit = BOX_NO_STRING, fmt = '(i5.5)') i

        !    open(unit = AUX_OUTPUT_UNIT, &
        !    &    file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
        !                trim(adjustl('PELAGIC_BOX_' // BOX_NO_STRING // 'SAVED_OUTPUTS.out')), &
        !    &    status = 'UNKNOWN', position = 'APPEND')


        !    do LINE_NO = 1, NUM_LINES

        !        TIME = PELAGIC_BOX_MODEL_DATA % SAVED_OUTPUT_RESULTS(LINE_NO, 0, 1)

        !        write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
        !        &     (PELAGIC_BOX_MODEL_DATA % SAVED_OUTPUT_RESULTS(LINE_NO, i, j), &
        !        &      j = 1, NUM_STATE_VARS)

        !    end do


        !    close(AUX_OUTPUT_UNIT)

        !end do

    end subroutine WRITE_PELAGIC_MEM_SAVED_OUTPUT


    subroutine WRITE_PELAGIC_EX_MEM_OUTPUT &
               (PELAGIC_BOX_MODEL_DATA, AUX_OUTPUT_UNIT, NUM_LINES)

        type(PELAGIC_BOX_MODEL_DS), intent(in) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: NUM_LINES

        real(kind = DBL) :: TIME

        character(len = 20) :: FORMAT_STRING
        character(len = 5)  :: INTEGER_STRING

        integer i
        integer j
        integer LINE_NO
        integer NUM_BOXES

        NUM_BOXES = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
        write(unit = INTEGER_STRING, fmt = '(i5)') NUM_BOXES + 1
        FORMAT_STRING = '(' // INTEGER_STRING // 'f20.6)'

        if (PELAGIC_BOX_MODEL_DATA % CREATE_PELAGIC_EXERGY_OUTPUTS > 0) then
            do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_EXERGY_COMPONENTS

                open(unit = AUX_OUTPUT_UNIT, &
                     file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                            trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                         EXERGY_COMPONENT_NAMES(i))) &
                            // '.out', status = 'UNKNOWN', position = 'APPEND')

                do LINE_NO = 1, NUM_LINES
                    TIME = PELAGIC_BOX_MODEL_DATA % EXERGY_RESULTS(LINE_NO, 0, i)

                    write(unit = AUX_OUTPUT_UNIT, fmt = FORMAT_STRING) TIME, &
                          (PELAGIC_BOX_MODEL_DATA % EXERGY_RESULTS(LINE_NO, j, i), &
                           j = 1, NUM_BOXES)
                end do

                close(AUX_OUTPUT_UNIT)
            end do
        end if
    end subroutine WRITE_PELAGIC_EX_MEM_OUTPUT


    subroutine WRITE_PELAGIC_MEM_SAVED_BINARY_OUTPUT &
               (PELAGIC_BOX_MODEL_DATA, AUX_OUTPUT_UNIT, NUM_LINES)

        type(PELAGIC_BOX_MODEL_DS), intent(in) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: NUM_LINES

        real(kind = DBL) :: TIME

        integer i
        integer j
        integer LINE_NO
        integer NUM_BOXES
        integer NUM_SAVED_OUTPUTS

        NUM_BOXES = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_SAVED_OUTPUTS

            open(unit = AUX_OUTPUT_UNIT, &
                 file = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                        trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" // &
				        trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                            SAVED_OUTPUT_NAMES(i))) // '.bin', &
                        status = 'UNKNOWN', access = 'STREAM')

            do LINE_NO = 1, NUM_LINES
                TIME = PELAGIC_BOX_MODEL_DATA % SAVED_OUTPUT_RESULTS(LINE_NO, 0, i)

                write(unit = AUX_OUTPUT_UNIT) TIME, &
                      (PELAGIC_BOX_MODEL_DATA % SAVED_OUTPUT_RESULTS(LINE_NO, j, i), &
                       j = 1, NUM_BOXES)
            end do

            close(AUX_OUTPUT_UNIT)
        end do
    end subroutine WRITE_PELAGIC_MEM_SAVED_BINARY_OUTPUT


    subroutine WRITE_PELAGIC_EX_MEM_BINARY_OUTPUT &
               (PELAGIC_BOX_MODEL_DATA, AUX_OUTPUT_UNIT, NUM_LINES)

        type(PELAGIC_BOX_MODEL_DS), intent(in) :: PELAGIC_BOX_MODEL_DATA
        integer, intent(in) :: AUX_OUTPUT_UNIT
        integer, intent(in) :: NUM_LINES

        real(kind = DBL) :: TIME

        integer i
        integer j
        integer LINE_NO
        integer NUM_BOXES

        NUM_BOXES = PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES

        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_EXERGY_COMPONENTS

            open(unit   = AUX_OUTPUT_UNIT, &
                 file   = trim(adjustl(PELAGIC_OUTPUT_FOLDER)) // &
                          trim(adjustl(BINARY_PELAGIC_OUTPUT_FILE_NAME)) // "_" // &
				          trim(adjustl(PELAGIC_BOX_MODEL_DATA % &
                                            EXERGY_COMPONENT_NAMES(i))) // '.bin', &
                 status = 'UNKNOWN', access = 'STREAM')

            do LINE_NO = 1, NUM_LINES
                TIME = PELAGIC_BOX_MODEL_DATA % EXERGY_RESULTS(LINE_NO, 0, i)

                write(unit = AUX_OUTPUT_UNIT) TIME, &
                      (PELAGIC_BOX_MODEL_DATA % EXERGY_RESULTS(LINE_NO, j, i), &
                       j = 1, NUM_BOXES)
            end do

            close(AUX_OUTPUT_UNIT)
        end do
    end subroutine WRITE_PELAGIC_EX_MEM_BINARY_OUTPUT

end module PELAGIC_BOX_MODEL
