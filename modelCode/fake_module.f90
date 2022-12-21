!==================================================================
module para_aqua
!==================================================================

!==================================================================
contains
!==================================================================
    subroutine para_get_value(name,value)
       character*(*) name
        double precision value
        value = 0.0D0
    end subroutine
end module para_aqua

    
function getpar(name)
    use para_aqua
    implicit none
    real getpar
    character*(*) name
    double precision value
    value = 0.0D0
    getpar = value
end