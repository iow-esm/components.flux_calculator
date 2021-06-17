MODULE flux_calculator_create_namcouple
    USE flux_calculator_basic
    USE flux_calculator_io

    IMPLICIT NONE

!!!!!!!!!! FUNCTIONS DEFINED IN THIS MODULE
    PUBLIC create_namcouple

    INTEGER, SAVE :: w_namcouple  ! a logfile to write the progress and error messages

!!!!!!!!!! NOW EVERYTHING ELSE
    CONTAINS

    SUBROUTINE write_header(num_input_fields, num_output_fields, timestep, num_timesteps)
        INTEGER, INTENT(IN)                                             :: num_input_fields, num_output_fields, timestep, num_timesteps

        WRITE(w_namcouple,*) '####################################################################'
        WRITE(w_namcouple,*) ' $NFIELDS'
        WRITE(w_namcouple,*) num_input_fields + num_output_fields
        WRITE(w_namcouple,*) ' $END'

        WRITE(w_namcouple,*) '############################################'
        WRITE(w_namcouple,*) ' $RUNTIME'
        WRITE(w_namcouple,*) timestep * num_timesteps
        WRITE(w_namcouple,*) ' $END'

        WRITE(w_namcouple,*) '############################################'
        WRITE(w_namcouple,*) ' $NLOGPRT'
        WRITE(w_namcouple,*) '1 1'
        WRITE(w_namcouple,*) ' $END'
        WRITE(w_namcouple,*) '############################################'
        WRITE(w_namcouple,*) ' $STRINGS'

    END SUBROUTINE write_header

    SUBROUTINE create_namcouple_entry(io_field, name_atmos_model, name_bottom_model, letter_bottom_model, timestep)
        TYPE(io_fields_type),    INTENT(IN)               :: io_field
        CHARACTER(*), INTENT(IN)                          :: name_atmos_model
        CHARACTER(*), DIMENSION(*), INTENT(IN)            :: name_bottom_model
        CHARACTER(*), DIMENSION(*), INTENT(IN)            :: letter_bottom_model
        INTEGER, INTENT(IN) :: timestep

        
        CHARACTER(len=32) :: my_model_name  ! name of the model form/to the filed is received/sent
        CHARACTER(len=8)  :: my_grid_name   ! grid on which the filed lives
        CHARACTER :: my_io  ! Role of flux_calculator variable (can be 'R' if the variable is received by the flux_calculator or 'S' if it is sent by flux_calculato)
        CHARACTER :: other_io ! Opposite of my_io
        CHARACTER (len=8) :: counterpart ! name of the variable that is received/sent from/by the model

        CHARACTER(len=128) :: mapping_file ! name of the mapping file used for this coupling entry
        TYPE(remapping_info_type) :: remapping_info ! info we get from the mapping file

        INTEGER :: i    ! counter

        ! find out from or to which model this variable is received or sent
        IF (io_field%name(2:2) == 'A') THEN
            my_model_name = name_atmos_model
        ELSE
            DO i = 1, MAX_BOTTOM_MODELS
                IF (io_field%name(2:2) == letter_bottom_model(i)) THEN
                    my_model_name = name_bottom_model(i)
                    EXIT
                ENDIF
            ENDDO
        ENDIF

        ! find out on which grid the field lives
        my_grid_name = grid_name(io_field%which_grid)
    
        ! if we receive, then this variable must be sent to us and vice versa
        my_io = io_field%name(1:1)
        IF (my_io == 'R') THEN
            other_io = 'S'
        ELSE
            other_io = 'R'
        ENDIF    
        
        ! contruct the name of the mapping file for this coupling
        IF (my_io == 'R') THEN
            ! if we receive from a model we have to map from model grid to exchange grid
            mapping_file = "mappings/remap_" // TRIM(my_grid_name) // "_" // TRIM(my_model_name) // "_to_exchangegrid.nc"
        ELSE
            ! if we send to a model we have to map from exchange grid to model grid
            mapping_file = "mappings/remap_" // TRIM(my_grid_name) // "_exchangegrid_to_" // TRIM(my_model_name) // ".nc"
        ENDIF
        
        ! get information on the mapping
        CALL read_remapping(mapping_file, remapping_info)

        ! construct counterpart of variable from/for the model
        counterpart = io_field%name
        counterpart(1:1) = io_field%name(2:2)
        counterpart(2:2) = other_io

        ! write entry
        IF (my_io == 'R') THEN
            WRITE(w_namcouple,*) counterpart, ' ', io_field%name, ' 1 ', timestep, ' 2 restart_flc_'//TRIM(io_field%name(3:6))//'_'//TRIM(my_model_name)//'.nc EXPOUT'
        ELSE
            WRITE(w_namcouple,*) io_field%name, ' ', counterpart, ' 1 ', timestep, ' 2 restart_flc_'//TRIM(io_field%name(3:6))//'_'//TRIM(my_model_name)//'.nc EXPOUT'
        ENDIF
        WRITE(w_namcouple,*) remapping_info%src_grid_dims, remapping_info%dst_grid_dims, '___ ___ LAG=0'   ! TODO: get rid off string literals here
        WRITE(w_namcouple,*) "R 0 R 0"
        WRITE(w_namcouple,*) "LOCTRANS MAPPING"
        WRITE(w_namcouple,*) "INSTANT"
        WRITE(w_namcouple,*) TRIM(mapping_file)
        WRITE(w_namcouple,*) "####"


    END SUBROUTINE create_namcouple_entry

    SUBROUTINE create_namcouple(input_field, num_input_fields, output_field, num_output_fields, &
                                name_atmos_model, name_bottom_model, letter_bottom_model,       &
                                timestep, num_timesteps)
        INTEGER, INTENT(IN)                                             :: num_input_fields, num_output_fields, timestep, num_timesteps
        TYPE(io_fields_type),    DIMENSION(*), INTENT(IN)               :: input_field
        TYPE(io_fields_type),    DIMENSION(*), INTENT(IN)               :: output_field
        CHARACTER(*), INTENT(IN)                                        :: name_atmos_model
        CHARACTER(*), DIMENSION(*), INTENT(IN)                          :: name_bottom_model
        CHARACTER(*), DIMENSION(*), INTENT(IN)                          :: letter_bottom_model

        INTEGER :: i 
        CHARACTER(len=16) :: namcouple_filename = "namcouple"

        INTEGER :: sys_status

        w_namcouple = 200 ! TODO find a better solution
        OPEN(w_namcouple,file=namcouple_filename,form='formatted')
        
        CALL write_header(num_input_fields, num_output_fields, timestep, num_timesteps)

        DO i = 1, num_output_fields
            CALL create_namcouple_entry(output_field(i), name_atmos_model, name_bottom_model, letter_bottom_model, timestep)
        ENDDO

        DO i = 1, num_input_fields
            CALL create_namcouple_entry(input_field(i), name_atmos_model, name_bottom_model, letter_bottom_model, timestep)
        ENDDO

        CLOSE(w_namcouple)

        CALL SYSTEM("cp "//namcouple_filename//" ../"//TRIM(name_atmos_model)//"/")
        WRITE(w_unit,*) "cp "//namcouple_filename//" ../"//TRIM(name_atmos_model)//"/"
        ! IF(sys_status /= 0) THEN
        !     WRITE(w_unit,*) "cp "//namcouple_filename//" ../"//name_atmos_model//"/"
        ! ENDIF
        DO i = 1, MAX_BOTTOM_MODELS
            IF (name_bottom_model(i) /= '') THEN
                CALL SYSTEM("cp "//namcouple_filename//" ../"//TRIM(name_bottom_model(i))//"/")
                WRITE(w_unit,*) "cp "//namcouple_filename//" ../"//TRIM(name_bottom_model(i))//"/"
        !         IF(sys_status /= 0) THEN
        !             WRITE(w_unit,*) "cp "//namcouple_filename//" ../"//name_bottom_model(i)//"/"
        !         ENDIF
            ENDIF
        ENDDO

    END SUBROUTINE create_namcouple

END MODULE flux_calculator_create_namcouple