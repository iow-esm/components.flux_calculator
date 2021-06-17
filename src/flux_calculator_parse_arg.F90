MODULE flux_calculator_parse_arg
    IMPLICIT NONE

!!!!!!!!!! FUNCTIONS DEFINED IN THIS MODULE
    PUBLIC get_args
    PUBLIC find_argument

!!!!!!!!!! NOW EVERYTHING ELSE

    INTEGER , SAVE :: num_args = -1
    CHARACTER(len=64), DIMENSION(:), ALLOCATABLE, SAVE :: args

    CONTAINS

    SUBROUTINE get_args()
        INTEGER :: ix
        
        IF (num_args /= -1) THEN
            RETURN
        ENDIF

        num_args = command_argument_count()

        IF (num_args == 0) THEN
            RETURN
        ENDIF

        ALLOCATE(args(num_args))

        DO ix = 1, num_args
            CALL get_command_argument(ix,args(ix))
        ENDDO

    END SUBROUTINE get_args

    FUNCTION find_argument(argument) RESULT(found)
        CHARACTER(*), INTENT(IN)     :: argument
        LOGICAL                      :: found

        INTEGER :: ix
        
        found = .FALSE.
        
        CALL get_args()

        DO ix = 1, num_args
            IF (TRIM(args(ix)) == TRIM(argument)) THEN
                found = .TRUE. 
                RETURN
            ENDIF
        ENDDO

    END FUNCTION find_argument

END MODULE flux_calculator_parse_arg