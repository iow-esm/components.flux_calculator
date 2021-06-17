MODULE flux_calculator_io
! This module contains functions for NetCDF IO which flux_calculator performs

    use flux_calculator_basic

    IMPLICIT NONE

    public read_scrip_grid_dimensions, read_regridding_matrix, read_remapping

    ! define a type to store info about remapping
    TYPE remapping_info_type
        INTEGER, DIMENSION(2)   :: src_grid_dims
        INTEGER, DIMENSION(2)   :: dst_grid_dims
    END TYPE remapping_info_type

    contains

    SUBROUTINE read_scrip_grid_dimensions(grid_filename, mype, npes, num_bottom_models, num_tasks_per_model, &
        num_grid_cells, grid_size_global, grid_size, grid_offset, grid_lon, grid_lat, grid_area)          

    ! Read a SCRIP-type NetCDF exchange grid file and return info on which gridcells we need to calculate 
    ! in this MPI instance.
        USE netcdf
                
        CHARACTER(len=50),       INTENT(IN)    :: grid_filename
        INTEGER,                 INTENT(IN)    :: mype
        INTEGER,                 INTENT(IN)    :: npes
        INTEGER,                 INTENT(IN)    :: num_bottom_models
        INTEGER, DIMENSION(:),   INTENT(IN)    :: num_tasks_per_model
        INTEGER, DIMENSION(:,:), INTENT(INOUT) :: num_grid_cells    ! this works in two ways:
                                                                    !   (a) num_bottom_models=1 and num_tasks_per_model=1, then this is an output field and its value (all points of this grid) is assigned here
                                                                    !   (b) otherwise this is an input field (that needs to be defined in the namelist and will be passed here)
        INTEGER,                 INTENT(OUT)   :: grid_size_global  ! grid size for the entire exchange grid (all bottom models)
        INTEGER,                 INTENT(OUT)   :: grid_size         ! grid size for this task only 
        INTEGER,                 INTENT(OUT)   :: grid_offset       ! sum of grid sizes from all preceding tasks 
                                                                    !   (each grid cell of the exchange grid is used by a single task only, so when one bottom model's exchange grid is split to different tasks, 
                                                                    !    some cells may need to be duplicated in the exchange grid file due to necessary regridding (u,v <-> t) from the "halos")
        TYPE(realarray2),        INTENT(INOUT) :: grid_lon
        TYPE(realarray2),        INTENT(INOUT) :: grid_lat
        TYPE(realarray2),        INTENT(INOUT) :: grid_area
        
        INTEGER                               :: nc              ! NetCDF file id
        INTEGER                               :: dimid_grid_size ! NetCDF dimension id
        INTEGER                               :: varid           ! NetCDF variable id
      
        CALL hdlerr(NF90_OPEN(grid_filename, NF90_NOWRITE, nc), __LINE__ )
        CALL hdlerr( NF90_INQ_DIMID(nc, 'grid_size' ,  dimid_grid_size),    __LINE__ )   ! get variable id
        CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid=nc,dimid=dimid_grid_size,len=grid_size_global), __LINE__ )

        ALLOCATE(grid_lon%field(grid_size_global,1));  grid_lon%allocated = .TRUE.
        ALLOCATE(grid_lat%field(grid_size_global,1));  grid_lat%allocated = .TRUE.
        ALLOCATE(grid_area%field(grid_size_global,1)); grid_area%allocated = .TRUE.
        
        CALL hdlerr( NF90_INQ_VARID(nc, 'grid_center_lon' ,  varid),    __LINE__ )                 ! get variable id
        CALL hdlerr( NF90_GET_VAR (nc, varid, grid_lon%field(:,1), [1], [grid_size_global]), __LINE__ ) ! get variable values
        CALL hdlerr( NF90_INQ_VARID(nc, 'grid_center_lat' ,  varid),    __LINE__ )                 
        CALL hdlerr( NF90_GET_VAR (nc, varid, grid_lat%field(:,1), [1], [grid_size_global]), __LINE__ ) 
        CALL hdlerr( NF90_INQ_VARID(nc, 'grid_area' ,  varid),    __LINE__ )                        
        CALL hdlerr( NF90_GET_VAR (nc, varid, grid_area%field(:,1), [1], [grid_size_global]), __LINE__ ) 

        CALL hdlerr(NF90_CLOSE(nc),    __LINE__ )

        ! most simple case: single model, single task -> use full grid
        IF (num_bottom_models == 1) THEN
            IF (num_tasks_per_model(1) == 1) THEN
                num_grid_cells      = 0
                num_grid_cells(1,1) = grid_size_global
                grid_size           = grid_size_global
                grid_offset         = 0
            ENDIF
        ENDIF

    END SUBROUTINE read_scrip_grid_dimensions

    SUBROUTINE read_regridding_matrix(regridding_filename,            &
                                      src_grid_size, src_grid_offset, &
                                      dst_grid_size, dst_grid_offset, &
                                      mymatrix)
        USE netcdf
        INCLUDE 'mpif.h'
        
    ! read the sparse regridding matrix (in coordinates format) for the full exchange grid and extract
    ! those elements which are for the current task
        CHARACTER(len=50),              INTENT(IN)    :: regridding_filename
        INTEGER,                        INTENT(IN)    :: src_grid_size         
        INTEGER,                        INTENT(IN)    :: src_grid_offset       
        INTEGER,                        INTENT(IN)    :: dst_grid_size         
        INTEGER,                        INTENT(IN)    :: dst_grid_offset
        TYPE(sparse_regridding_matrix), INTENT(INOUT) :: mymatrix
        
        INTEGER                               :: i, num_total, num_matching
        INTEGER                               :: nc              ! NetCDF file id
        INTEGER                               :: dimid_num_links ! NetCDF dimension id
        INTEGER                               :: varid           ! NetCDF variable id
        TYPE(integerarray)                    :: all_src_index
        TYPE(integerarray)                    :: all_dst_index
        TYPE(realarray2)                      :: all_weight
        
        ! read in the overall exchange grid remapping matrix from the file
        CALL hdlerr(NF90_OPEN(regridding_filename, NF90_NOWRITE, nc), __LINE__ )
        CALL hdlerr( NF90_INQ_DIMID(nc, 'num_links' ,  dimid_num_links),    __LINE__ )   ! get dimension id
        CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid=nc,dimid=dimid_num_links,len=num_total), __LINE__ )
        ALLOCATE(all_src_index%field(num_total));   all_src_index%allocated = .TRUE.
        ALLOCATE(all_dst_index%field(num_total));   all_dst_index%allocated = .TRUE.
        ALLOCATE(all_weight%field(1,num_total));    all_weight%allocated = .TRUE.
        CALL hdlerr( NF90_INQ_VARID(nc, 'src_address' ,  varid),    __LINE__ )                      ! get variable id
        CALL hdlerr( NF90_GET_VAR (nc, varid, all_src_index%field(:), [1], [num_total]), __LINE__ ) ! get variable values
        CALL hdlerr( NF90_INQ_VARID(nc, 'dst_address' ,  varid),    __LINE__ )                      ! get variable id
        CALL hdlerr( NF90_GET_VAR (nc, varid, all_dst_index%field(:), [1], [num_total]), __LINE__ ) ! get variable values
        CALL hdlerr( NF90_INQ_VARID(nc, 'remap_matrix' ,  varid),    __LINE__ )                     ! get variable id
        CALL hdlerr( NF90_GET_VAR (nc, varid, all_weight%field(:,:), [1,1], [1,num_total]), __LINE__ )  ! get variable values
        CALL hdlerr(NF90_CLOSE(nc),    __LINE__ )

        ! find out how many entries are matching, so we can allocate the final array to the correct size
        num_matching = 0
        DO i=1,num_total
            IF ((all_dst_index%field(i) .gt. dst_grid_offset) .AND. &
            (all_dst_index%field(i) .le. dst_grid_offset + dst_grid_size)) THEN
                num_matching = num_matching + 1
            END IF
        END DO
        mymatrix%num_elements = num_matching
        ALLOCATE(mymatrix%src_index%field(num_matching));   mymatrix%src_index%allocated = .TRUE.
        ALLOCATE(mymatrix%dst_index%field(num_matching));   mymatrix%dst_index%allocated = .TRUE.
        ALLOCATE(mymatrix%weight%field(num_matching));      mymatrix%weight%allocated = .TRUE.
        
        ! walk through the full matrix and copy only the matching entries to the matrix that will be used in this task
        num_matching = 0
        DO i=1,num_total
            IF ((all_dst_index%field(i) .gt. dst_grid_offset) .AND. &
            (all_dst_index%field(i) .le. dst_grid_offset + dst_grid_size)) THEN
                num_matching = num_matching + 1
                mymatrix%dst_index%field(num_matching)   = all_dst_index%field(i)
                mymatrix%src_index%field(num_matching)   = all_src_index%field(i)
                mymatrix%weight%field(num_matching)      = all_weight%field(1,i)
            END IF
        END DO

        ! correct for the offset and check if all source cells are present at this PE
        DO i=1,num_matching
            mymatrix%dst_index%field(i) = mymatrix%dst_index%field(i) - dst_grid_offset
            mymatrix%src_index%field(i) = mymatrix%src_index%field(i) - src_grid_offset
            IF ((mymatrix%src_index%field(i) .LT. 1) .OR. (mymatrix%src_index%field(i) .GT. src_grid_size)) THEN
                write ( * , * ) 'Regridding matrix did not match task decomposition in file ',regridding_filename
                write ( * , * ) 'Stopped '
                call MPI_Abort ( MPI_COMM_WORLD, 1, 1 )
            END IF
        END DO

        ! clean up: deallocate full arrays
        DEALLOCATE(all_src_index%field)
        DEALLOCATE(all_dst_index%field)
        DEALLOCATE(all_weight%field)
      
    END SUBROUTINE read_regridding_matrix

    SUBROUTINE read_remapping(remapping_filename, remapping_info)
        
        USE netcdf
        
        CHARACTER(len=128),              INTENT(IN)     :: remapping_filename
        TYPE(remapping_info_type),       INTENT(OUT)    :: remapping_info

        INTEGER                               :: nc              ! NetCDF file id
        INTEGER                               :: varid           ! NetCDF variable id

        INTEGER :: src_grid_rank, dst_grid_rank

        ! read in the dimensions for remapping
        CALL hdlerr(NF90_OPEN(remapping_filename, NF90_NOWRITE, nc), __LINE__ )

        CALL hdlerr( NF90_INQ_DIMID(nc, 'src_grid_rank' ,  varid),    __LINE__ )                      ! get variable id
        CALL hdlerr( NF90_INQUIRE_DIMENSION(nc, varid, len=src_grid_rank), __LINE__ ) ! get variable values
        WRITE (*,*) 'src_grid_rank ', src_grid_rank

        CALL hdlerr( NF90_INQ_VARID(nc, 'src_grid_dims' ,  varid),    __LINE__ )                      ! get variable id
        CALL hdlerr( NF90_GET_VAR (nc, varid, remapping_info%src_grid_dims(:), [1], [src_grid_rank]), __LINE__ ) ! get variable values
        IF (src_grid_rank == 1) THEN
            remapping_info%src_grid_dims(2) = 1
        ENDIF

        CALL hdlerr( NF90_INQ_DIMID(nc, 'dst_grid_rank' ,  varid),    __LINE__ )                      ! get variable id
        CALL hdlerr( NF90_INQUIRE_DIMENSION(nc, varid, len=dst_grid_rank), __LINE__ ) ! get variable values

        CALL hdlerr( NF90_INQ_VARID(nc, 'dst_grid_dims' ,  varid),    __LINE__ )                      ! get variable id 
        CALL hdlerr( NF90_GET_VAR (nc, varid, remapping_info%dst_grid_dims(:), [1], [dst_grid_rank]), __LINE__ ) ! get variable values
        IF (dst_grid_rank == 1) THEN
            remapping_info%dst_grid_dims(2) = 1
        ENDIF
        
        CALL hdlerr(NF90_CLOSE(nc),    __LINE__ )


    END SUBROUTINE read_remapping

    SUBROUTINE hdlerr(istatus, line)
    ! handle errors during NetCDF calls
        use netcdf
        implicit none
        !
        INCLUDE 'mpif.h'
        
        integer, intent(in) :: istatus, line
        !
        IF (istatus .NE. NF90_NOERR) THEN
            write ( * , * ) 'NetCDF problem at line',line
            write ( * , * ) 'Stopped '
            call MPI_Abort ( MPI_COMM_WORLD, 1, 1 )
        ENDIF
        !
        RETURN
    END SUBROUTINE hdlerr

END MODULE flux_calculator_io