!****************************************************************************************
SUBROUTINE read_dimgrid (nlon,nlat,data_filename,w_unit)
  !**************************************************************************************
  USE netcdf
  IMPLICIT NONE
  !
  INTEGER                  :: i,j,k,w_unit
  !
  INTEGER                  :: il_file_id,il_grid_id,il_lon_id, &
     il_lat_id,il_indice_id, &
     lon_dims,lat_dims,imask_dims
  !
  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: lon_dims_ids,lat_dims_ids,&
     imask_dims_ids,lon_dims_len,&
     lat_dims_len,imask_dims_len  
  !               
  INTEGER, INTENT(out)     :: nlon,nlat
  !
  CHARACTER(len=50)        :: data_filename
  !
  ! Dimensions
  !
  CALL hdlerr(NF90_OPEN(data_filename, NF90_NOWRITE, il_file_id), __LINE__ )
  !
  CALL hdlerr( NF90_INQ_VARID(il_file_id, 'lon' ,  il_lon_id),    __LINE__ )
  CALL hdlerr( NF90_INQ_VARID(il_file_id, 'lat' ,  il_lat_id),    __LINE__ )
  CALL hdlerr( NF90_INQ_VARID(il_file_id, 'imask', il_indice_id), __LINE__ )
  !
  CALL hdlerr( NF90_INQUIRE_VARIABLE(il_file_id, varid=il_lon_id, ndims=lon_dims, dimids=lon_dims_ids), __LINE__ )
  !
  ! The value lon_dims_len(i) is obtained thanks to the lon_dims_ids ID already obtained from the file
  DO i=1,lon_dims
    CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid=il_file_id,dimid=lon_dims_ids(i),len=lon_dims_len(i)), __LINE__ )
  ENDDO
  !
  nlon=lon_dims_len(1)
  nlat=lon_dims_len(2)
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  CALL hdlerr( NF90_INQUIRE_VARIABLE(ncid=il_file_id, varid=il_lat_id, ndims=lat_dims, dimids=lat_dims_ids), __LINE__ )
  !
  ! The value lat_dims_len(i) is obtained thanks to the lat_dims_ids ID already obtained from the file
  DO i=1,lat_dims
    CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid=il_file_id,dimid=lat_dims_ids(i),len=lat_dims_len(i)), __LINE__ )
  ENDDO
  !
  IF ( (lat_dims_len(1) .NE. lon_dims_len(1)).OR.(lat_dims_len(2) .NE. lon_dims_len(2)) ) THEN
      WRITE(w_unit,*) 'Problem model1 in read_dimgrid'
      WRITE(w_unit,*) 'Dimensions of the latitude are not the same as the ones of the longitude'
      CALL flush(w_unit)
      STOP
  ENDIF
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  CALL hdlerr( NF90_INQUIRE_VARIABLE(ncid=il_file_id, varid=il_indice_id, ndims=imask_dims, dimids=imask_dims_ids), __LINE__ )
  !
  ! The value imask_dims_len(i) is obtained thanks to the imask_dims_ids ID already obtained from the file
  DO i=1,imask_dims
    CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid=il_file_id,dimid=imask_dims_ids(i),len=imask_dims_len(i)), __LINE__ )
  ENDDO
  !
  CALL hdlerr(NF90_CLOSE(il_file_id), __LINE__ )
  !
  WRITE(w_unit,*) 'Reading input file ',data_filename
  WRITE(w_unit,*) 'Global dimensions nlon=',nlon,' nlat=',nlat
  CALL flush(w_unit)
  !
  !
END SUBROUTINE read_dimgrid

SUBROUTINE read_scrip_grid_dimensions(grid_filename, mype, npes, num_bottom_models, num_tasks_per_model, &
                                      num_grid_cells, grid_size_global, grid_size, grid_offset)          

  ! Read a SCRIP-type NetCDF exchange grid file and return info on which gridcells we need to calculate 
  ! in this MPI instance.
  USE netcdf
  IMPLICIT NONE
                                      
  CHARACTER(len=50),       INTENT(IN)    :: grid_filename
  INTEGER,                 INTENT(IN)    :: mype
  INTEGER,                 INTENT(IN)    :: npes
  INTEGER,                 INTENT(IN)    :: num_bottom_models
  INTEGER, DIMENSION(:),   INTENT(IN)    :: num_tasks_per_model
  INTEGER, DIMENSION(:,:), INTENT(INOUT) :: num_grid_cells
  INTEGER,                 INTENT(OUT)   :: grid_size_global
  INTEGER,                 INTENT(OUT)   :: grid_size
  INTEGER,                 INTENT(OUT)   :: grid_offset
  INTEGER,                 INTENT(IN)    :: w_unit

  INTEGER                               :: nc              ! NetCDF file id
  INTEGER                               :: dimid_grid_size ! NetCDF dimension id
  !INTEGER                               :: ndims           ! number of dimensions
  !INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimids          ! array if dimensions
  !INTEGER                               :: dimid_GRID_SIZE ! NetCDF dimension id
  
  CALL hdlerr(NF90_OPEN(grid_filename, NF90_NOWRITE, nc), __LINE__ )
  CALL hdlerr( NF90_INQ_DIMID(nc, 'grid_size' ,  dimid_grid_size),    __LINE__ )   ! get variable id
  !CALL hdlerr( NF90_INQUIRE_VARIABLE(nc, varid=varid_GRID_SIZE, ndims=ndims, dimids=dimids), __LINE__ )
  CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid=nc,dimid=dimid_grid_size,len=grid_size_global), __LINE__ )
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