module PELAGIC_EXERGY

    use PELAGIC_BOX_MODEL
    use aquabc_pel_state_var_indexes
    implicit none

contains


    subroutine CALCULATE_PELAGIC_EXERGY(STATE_VARIABLES, TEMP, EXERGY_COMPONENTS, TIME, BOX_NO)
        real(kind = SELECTED_REAL_KIND(15, 307)), dimension(:), intent(in)  :: STATE_VARIABLES
        real(kind = SELECTED_REAL_KIND(15, 307)), dimension(:), intent(out) :: EXERGY_COMPONENTS
        real(kind = SELECTED_REAL_KIND(15, 307)), intent(in) :: TEMP
        real(kind = SELECTED_REAL_KIND(15, 307)), intent(in) :: TIME 

        integer, intent(in) :: BOX_NO

        !State variables
        real(kind = SELECTED_REAL_KIND(15, 307)) :: NH4_N
        real(kind = SELECTED_REAL_KIND(15, 307)) :: NO3_N
        real(kind = SELECTED_REAL_KIND(15, 307)) :: PO4_P
        real(kind = SELECTED_REAL_KIND(15, 307)) :: DISS_OXYGEN
        real(kind = SELECTED_REAL_KIND(15, 307)) :: PHYTO_C
        real(kind = SELECTED_REAL_KIND(15, 307)) :: ZOO_C

        !Exergy variables
        real(kind = SELECTED_REAL_KIND(15,307)) :: EX_DEAD_ORGANIC_CARBON
        real(kind = SELECTED_REAL_KIND(15,307)) :: EX_BACTERIA_CARBON
        real(kind = SELECTED_REAL_KIND(15,307)) :: EX_PHYTOPLANKTON_CARBON
        real(kind = SELECTED_REAL_KIND(15,307)) :: EX_ZOOPLANKTON_CARBON

        !Auxillary variables
        real(kind = SELECTED_REAL_KIND(15,307)) :: TEMP_KELVIN
        real(kind = SELECTED_REAL_KIND(15,307)) :: DEAD_ORGANIC_CARBON
        real(kind = SELECTED_REAL_KIND(15,307)) :: BACTERIA_CARBON
        real(kind = SELECTED_REAL_KIND(15,307)) :: PHYTOPLANKTON_CARBON
        real(kind = SELECTED_REAL_KIND(15,307)) :: ZOOPLANKTON_CARBON


        !Initialize state variables
        NH4_N          = STATE_VARIABLES(NH4_N_INDEX)
        NO3_N          = STATE_VARIABLES(NO3_N_INDEX)
        PO4_P          = STATE_VARIABLES(PO4_P_INDEX)  
        DISS_OXYGEN    = STATE_VARIABLES(DISS_OXYGEN_INDEX)
        
        PHYTO_C        = &
            STATE_VARIABLES(DIA_C_INDEX) + STATE_VARIABLES(CYN_C_INDEX) + &
            STATE_VARIABLES(OPA_C_INDEX) + STATE_VARIABLES(FIX_CYN_C_INDEX)
        
        ZOO_C          = STATE_VARIABLES(ZOO_C_INDEX)

        TEMP_KELVIN = TEMP + 2.7315D2

        DEAD_ORGANIC_CARBON  = STATE_VARIABLES(DET_PART_ORG_C_INDEX)
        BACTERIA_CARBON      = 0.0D0
        PHYTOPLANKTON_CARBON = PHYTO_C
        ZOOPLANKTON_CARBON   = ZOO_C

        EX_DEAD_ORGANIC_CARBON     = TEMP_KELVIN * (6.23D1 + (1.74D-1 * 0.0D0 )) * DEAD_ORGANIC_CARBON
        EX_BACTERIA_CARBON         = TEMP_KELVIN * (6.23D1 + (1.74D-1 * 6.0D0 )) * BACTERIA_CARBON
        EX_PHYTOPLANKTON_CARBON    = TEMP_KELVIN * (6.23D1 + (1.74D-1 * 8.5D2 )) * PHYTOPLANKTON_CARBON
        EX_ZOOPLANKTON_CARBON      = TEMP_KELVIN * (6.23D1 + (1.74D-1 * 1.25D4)) * ZOOPLANKTON_CARBON
        EXERGY_COMPONENTS(1) = EX_DEAD_ORGANIC_CARBON
        EXERGY_COMPONENTS(2) = EX_BACTERIA_CARBON
        EXERGY_COMPONENTS(3) = EX_PHYTOPLANKTON_CARBON
        EXERGY_COMPONENTS(4) = EX_ZOOPLANKTON_CARBON
    end subroutine CALCULATE_PELAGIC_EXERGY


    subroutine ALLOCATE_PELAGIC_EXERGY(PELAGIC_BOX_MODEL_DATA)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA

        integer :: i

        integer :: NUM_EXERGY_COMPONENTS

        NUM_EXERGY_COMPONENTS = 4


        allocate(PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENTS(NUM_EXERGY_COMPONENTS))
        allocate(PELAGIC_BOX_MODEL_DATA % AVERAGED_EXERGY_COMPONENTS(NUM_EXERGY_COMPONENTS))
        allocate(PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENT_NAMES(NUM_EXERGY_COMPONENTS))

        PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENT_NAMES(1) = 'EX_DEAD_ORGANIC_CARBON'
        PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENT_NAMES(2) = 'EX_BACTERIA_CARBON'
        PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENT_NAMES(3) = 'EX_PHYTOPLANKTON_CARBON'
        PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENT_NAMES(4) = 'EX_ZOOPLANKTON_CARBON'


        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            allocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % EXERGY_COMPONENTS(NUM_EXERGY_COMPONENTS))
            allocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % AVERAGED_EXERGY_COMPONENTS(NUM_EXERGY_COMPONENTS))
            PELAGIC_BOX_MODEL_DATA % NUM_EXERGY_COMPONENTS = NUM_EXERGY_COMPONENTS
        end do

    end subroutine ALLOCATE_PELAGIC_EXERGY



    subroutine DEALLOCATE_PELAGIC_EXERGY(PELAGIC_BOX_MODEL_DATA)

        type(PELAGIC_BOX_MODEL_DS), intent(inout) :: PELAGIC_BOX_MODEL_DATA

        integer :: i

        deallocate(PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENTS)
        deallocate(PELAGIC_BOX_MODEL_DATA % AVERAGED_EXERGY_COMPONENTS)
        deallocate(PELAGIC_BOX_MODEL_DATA % EXERGY_COMPONENT_NAMES)


        do i = 1, PELAGIC_BOX_MODEL_DATA % NUM_PELAGIC_BOXES
            deallocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % EXERGY_COMPONENTS)
            deallocate(PELAGIC_BOX_MODEL_DATA % PELAGIC_BOXES(i) % AVERAGED_EXERGY_COMPONENTS)
            PELAGIC_BOX_MODEL_DATA % NUM_EXERGY_COMPONENTS = 0
        end do

    end subroutine DEALLOCATE_PELAGIC_EXERGY

end module PELAGIC_EXERGY