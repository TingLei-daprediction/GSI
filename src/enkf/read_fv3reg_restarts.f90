    module read_fv3regional_restarts
! modified from read_fv3_restarts.f90

! ifort -I${NETCDF}/include -O2 -traceback read_fv3_restarts.f90 kinds.o
! netcdf_mod.o -L/${NETCDF}/lib -lnetcdf -lnetcdff

! read data from FV3 restart files.


    use kinds, only: i_kind,r_single,r_kind
    use netcdf, only: nf90_open,nf90_close,nf90_get_var,nf90_noerr
    use netcdf, only: nf90_inq_dimid,nf90_inq_varid
    use netcdf, only: nf90_nowrite,nf90_inquire,nf90_inquire_dimension
    use netcdf_mod, only: nc_check
    public read_fv3_restart_data1d,read_fv3_restart_data2d
    public read_fv3_restart_data3d,read_fv3_restart_data4d

    contains

    subroutine read_fv3_restart_data1d(varname,filename,file_id,data_arr)
       real(r_single), intent(inout), dimension(:) :: data_arr
       character(len=24),parameter :: myname_ = 'read_fv3_restart_data1d'
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       integer(i_kind), intent(in) :: file_id
       integer(i_kind) :: var_id
       call nc_check( nf90_inq_varid(file_id,trim(adjustl(varname)),var_id),&
       myname_,'inq_varid '//trim(adjustl(varname))//' '//trim(filename) )
       call nc_check( nf90_get_var(file_id,var_id,data_arr),&
       myname_,'get_var '//trim(adjustl(varname))//' '//trim(filename) )
       data_arr=data_arr(ubound(data_arr,1):lbound(data_arr,1):-1)
    end subroutine read_fv3_restart_data1d

    subroutine read_fv3_restart_data2d(varname,filename,file_id,data_arr)
       real(r_single), intent(inout), dimension(:,:) :: data_arr
       character(len=24),parameter :: myname_ = 'read_fv3_restart_data2d'
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       integer(i_kind), intent(in) :: file_id
       integer(i_kind) :: var_id
       call nc_check( nf90_inq_varid(file_id,trim(adjustl(varname)),var_id),&
       myname_,'inq_varid '//trim(adjustl(varname))//' '//trim(filename) )
       call nc_check( nf90_get_var(file_id,var_id,data_arr),&
       myname_,'get_var '//trim(adjustl(varname))//' '//trim(filename) )
    end subroutine read_fv3_restart_data2d

    subroutine read_fv3_restart_data3d(varname,filename,file_id,data_arr)
       real(r_single), intent(inout), dimension(:,:,:) :: data_arr
       character(len=24),parameter :: myname_ = 'read_fv3_restart_data3d'
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       integer(i_kind), intent(in) :: file_id
       real(r_single) :: rtmp
       integer(i_kind) :: var_id
       integer(i_kind)::i,j,k,it
       integer(i_kind):: ilow,iup, jlow,jup,klow,kup,klev
       ilow=lbound(data_arr,1)
       iup=ubound(data_arr,1)
       jlow=lbound(data_arr,2)
       jup=ubound(data_arr,2)
       klow=lbound(data_arr,3)
       kup=ubound(data_arr,3)
       call nc_check( nf90_inq_varid(file_id,trim(adjustl(varname)),var_id),&
       myname_,'inq_varid '//trim(adjustl(varname))//' '//trim(filename) )
       call nc_check( nf90_get_var(file_id,var_id,data_tmp),&
       myname_,'get_var '//trim(adjustl(varname))//' '//trim(filename) )
!$omp parallel do schedule(dynamic,1) private(k,j,i,klev,rtmp)
       do k=klow, klow+ ceiling((kup-klow+1)/2.0)
          klev=kup-k+klow
          do j=jlow, jup
            do i=ilow, iup
             rtmp=data_arr(i,j,k)
             data_arr(i,j,k)=data_arr(i,j,klev)
             data_arr(i,j,klev)=rtmp
            enddo
         enddo
       enddo 

!       data_arr=data_arr(:,:, &
!                          ubound(data_arr,3):lbound(data_arr,3):-1)
    end subroutine read_fv3_restart_data3d

    subroutine read_fv3_restart_data4d(varname,filename,file_id,data_arr)
       real(r_single), intent(inout), dimension(:,:,:,:) :: data_arr
       character(len=24),parameter :: myname_ = 'read_fv3_restart_data4d'
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       integer(i_kind), intent(in) :: file_id
       integer(i_kind) :: var_id
       real(r_single)  :: rtmp
       integer(i_kind)::i,j,k,it
       integer(i_kind):: ilow,iup, jlow,jup,klow,kup,tlow,tup,klev
       ilow=lbound(data_arr,1)
       iup=ubound(data_arr,1)
       jlow=lbound(data_arr,2)
       jup=ubound(data_arr,2)
       klow=lbound(data_arr,3)
       kup=ubound(data_arr,3)
       tlow=lbound(data_arr,4)
       tup=ubound(data_arr,4)
       call nc_check( nf90_inq_varid(file_id,trim(adjustl(varname)),var_id),&
       myname_,'inq_varid '//trim(adjustl(varname))//' '//trim(filename) )
       call nc_check( nf90_get_var(file_id,var_id,data_tmp),&
       myname_,'get_var '//trim(adjustl(varname))//' '//trim(filename) )
       do it=tlow, tup  !normally, this dimnsions is very small or even 1, 
                        ! hence, the openmp parallelization is started from next
                        ! loop
!$omp parallel do schedule(dynamic,1) private(k,j,i,klev,rtmp)
        do k=klow, klow+ ceiling((kup-klow+1)/2.0)
          klev=kup-k+klow
          do j=jlow, jup
            do i=ilow, iup
             rtmp=data_arr(i,j,k,it)
             data_arr(i,j,k,it)=data_arr(i,j,klev,it)
             data_arr(i,j,klev,it)=rtmp
            enddo
         enddo
        enddo 
       enddo 
!       data_arr=data_arr(:,:, &
!                          ubound(data_arr,3):lbound(data_arr,3):-1,:)
    end subroutine read_fv3_restart_data4d

    end module read_fv3regional_restarts
