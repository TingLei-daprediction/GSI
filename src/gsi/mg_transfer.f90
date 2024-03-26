!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_transfer 
!***********************************************************************
!                                                                      !
!  Transfer data between analysis and filter grid                      !
!                                                                      !
! Modules: kinds, mg_parameter, mg_intstate, mg_bocos, mg_interpolate, !
!          mg_timers, mg_mppstuff                                      !
!                                                     M. Rancic (2021) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_parameter
!STA use mg_intstate, only: WORKA
use mg_intstate, only: VALL
use mg_timers
use mg_mppstuff, only: nx,my,mpi_comm_comp,ninc2,minc2
use mpimod, only: mype
use mg_mppstuff, only: finishMPI
use pre_multigrid, only: WORKA
use gridmod, only: lat1,lon1,lat2,lon2,nlat,nlon

implicit none
integer(i_kind):: n,m,l,k,i,j

public anal_to_filt
public filt_to_anal

public stack_to_composite
public composite_to_stack

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_adjoint_offset
use mg_bocos, only:  bocoT_2d
implicit none


real(r_kind),allocatable,dimension(:,:):: WORKAT
real(r_kind),allocatable,dimension(:,:,:):: WORKS

real(r_kind),allocatable,dimension(:,:,:):: VLOC  

!----------------------------------------------------------------------
!TEST
!   if(mype==0) then
!     write(0,*)'FROM ANAL_TO_FILT: lat2,lon2,nm,mm=',lat2,lon2,nm,mm
!   endif
!   call finishMPI
!TEST


    allocate(WORKAT(lon2,lat2))
    allocate(WORKS(km,1:nm,1:mm))

  do k=1,km
    do m=1,lat2
    do n=1,lon2
      WORKAT(n,m)=WORKA(k,m,n)
    enddo
    enddo
    do m=1,mm
    do n=1,nm
      WORKS(k,n,m)=WORKAT(n+1,m+1)
    enddo
    enddo
  enddo

    deallocate(WORKAT)

    allocate(VLOC(km,1-ib:im+ib,1-jb:jm+jb))                      


!T                                                 call btim(  aintp_tim)

      VLOC=0.
         call lsqr_adjoint_offset(WORKS,VLOC,km)


!T                                                 call etim(  aintp_tim)


!***
!***  Apply adjoint of lateral bc 
!***
    

         call bocoT_2d(VLOC,km,im,jm,ib,jb)
 
       VALL=0.
       VALL(1:km,1:im,1:jm)=VLOC(1:km,1:im,1:jm)
      

    deallocate(VLOC)

!                                            call etim(   btrns1_tim)

!----------------------------------------------------------------------
                        endsubroutine anal_to_filt

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_direct_offset
use mg_bocos, only:  boco_2d
implicit none


real(r_kind),allocatable,dimension(:,:,:):: WORKS
real(r_kind),allocatable,dimension(:,:,:):: WORKAT

real(r_kind),allocatable,dimension(:,:,:):: VLOC   


!----------------------------------------------------------------------

!T                                            call btim(   btrns2_tim)

!***
!***  Define VLOC
!***

    allocate(VLOC(1:km,1-ib:im+ib,1-jb:jm+jb))                     

      VLOC=0.
      VLOC(1:km,1:im,1:jm)=VALL(1:km,1:im,1:jm)
        

!***
!***  Supply boundary conditions for VLOC
!***
         call boco_2d(VLOC,km,im,jm,ib,jb)


!***
!*** Interpolate to analysis grid composite variables
!***

       allocate(WORKS(km,1:nm,1:mm))

!T                                                 call btim(   intp_tim)

         call lsqr_direct_offset(VLOC,WORKS,km)

!T                                                 call etim(   intp_tim)
    deallocate(VLOC)

!***
!*** Return WORKS to WORKA
!***

    allocate(WORKAT(km,0:nm+1,0:mm+1))


     WORKAT(:,:,:)=0.
  do k=1,km
    do m=1,mm
    do n=1,nm
      WORKAT(k,n,m) = WORKS(k,n,m)
    enddo
    enddo
  enddo


                     call boco_2d(WORKAT,km,nm,mm,1,1)

  do k=1,km
    do m=1,lat2
    do n=1,lon2
      WORKA(k,m,n)=WORKAT(k,n-1,m-1)
    enddo
    enddo
  enddo


    deallocate(WORKS)
    deallocate(WORKAT)




!T                                                 call etim(   btrns2_tim)

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine stack_to_composite                   &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(ARR_ALL,A2D,A3D)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km ,1-hx:im+hx,1-hy:jm+hy),   intent(in):: ARR_ALL
real(r_kind),dimension(km3,1-hx:im+hx,1-hy:jm+hy,lm),intent(out):: A3D
real(r_kind),dimension(km2,1-hx:im+hx,1-hy:jm+hy)   ,intent(out):: A2D
!----------------------------------------------------------------------
    do k=1,km3
      do L=1,lm
        A3D(k,:,:,L)=ARR_ALL((k-1)*lm+L,:,:)
      enddo
    enddo

    do k=1,km2
      A2D(k,:,:)=ARR_ALL(km3*lm+k,:,:)
    enddo

!----------------------------------------------------------------------
                        endsubroutine stack_to_composite

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine composite_to_stack                   &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(A2D,A3D,ARR_ALL)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km2,1-hx:im+hx,1-hy:jm+hy),   intent(in):: A2D
real(r_kind),dimension(km3,1-hx:im+hx,1-hy:jm+hy,lm),intent(in):: A3D
real(r_kind),dimension(km ,1-hx:im+hx,1-hy:jm+hy),   intent(out):: ARR_ALL
!----------------------------------------------------------------------
  do k=1,km3
    do L=1,lm
      ARR_ALL((k-1)*lm+L,:,:)=A3D(k,:,:,L)
    enddo
  enddo

  do k=1,km2
    ARR_ALL(km3*lm+k,:,:)= A2D(k,:,:)
  enddo

!----------------------------------------------------------------------
                        endsubroutine composite_to_stack 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_transfer
