!----------------------------------------------------------------------------------
module pre_multigrid
!$$$ module documentation block
!           .      .    .                                       .
! module:   pre_multigrid
!   prgmmr: pondeca          org: np23                date: 2020-09-28
!
! abstract: convert gradient vector in gsi_bundle format to 3D and 2D
!           fields that can be more conveniently accessed by the multigrid-beta
!           function code for bckg error covariances
!
! program history log:
!   2020-11-04  pondeca
!
!
! subroutines included:
!   sub init_transf_flds
!   sub mk_transf_fields
!   sub destroy_transf_flds
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
! subroutines included:
!   sub init_transf_flds
!   sub mk_transf_fields
!   sub destroy_transf_flds
!
!$$$ end documentation block

  use kinds, only: r_kind,i_kind
!!  use mg_intstate, only: WORKA

  implicit none

  integer(i_kind), parameter :: MAXSTR=256

  integer(i_kind),allocatable:: nvdim(:)
  integer(i_kind),allocatable:: latglobal_index(:)
  integer(i_kind),allocatable:: longlobal_index(:)
  integer(i_kind),allocatable:: latdecomp12(:,:)
  integer(i_kind),allocatable:: londecomp12(:,:)

  real(r_kind),allocatable:: worka(:,:,:)

  character(len=MAXSTR),allocatable:: vname_in_worka(:)
  character(len=MAXSTR),allocatable:: vorder_in_worka(:)

  real(r_kind),allocatable:: zfld(:,:)
 
! set default to private
  private
! set subroutines to public
  public :: init_transf_flds
  public :: latdecomp12 
  public :: londecomp12 
  public :: mk_transf_fields
  public :: destroy_transf_flds
  public :: nvdim
  public :: latglobal_index
  public :: longlobal_index
  public :: worka
  public :: vname_in_worka
  public :: vorder_in_worka
  public :: init_mg_get_background
  public :: mg_get_background
  public :: destroy_mg_get_background
  public :: mkgrads4ctrvec
  public :: zfld

contains

subroutine init_transf_flds(cvec)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    create_transf_flds 
!   prgmmr: pondeca          org: np22                date: 2020-09-28
!
! abstract: 
!
! program history log:
!   2020-09-28  pondeca
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use gsi_bundlemod, only : gsi_bundle
  use gridmod, only: nsig,lat2,lon2,lat1,lon1,nlat,nlon,istart,jstart
  use mpimod, only: mype,npe,ierror,mpi_itype,mpi_sum,mpi_comm_world
  implicit none

! Declare passed variables
  type(gsi_bundle),intent(in) :: cvec

! Declare local variables
  integer(i_kind) i,j,n,k,kk,mm1
  integer(i_kind):: num_fields
  integer(i_kind),allocatable:: latdecomp12_aux(:,:)
  integer(i_kind),allocatable:: londecomp12_aux(:,:)

  num_fields=cvec%n3d*nsig+cvec%n2d
  if (mype==0) print*,'in create_transf_flds: num_fields=',num_fields

!cltorg   allocate(worka(lat2,lon2,num_fields))
  allocate(worka(num_fields,lat1,lon1))
  allocate(vname_in_worka(cvec%n3d*nsig+cvec%n2d))
  allocate(vorder_in_worka(cvec%n3d+cvec%n2d))
  allocate(nvdim(cvec%n3d+cvec%n2d))
  allocate(latglobal_index(lat1))
  allocate(longlobal_index(lon1))
  allocate(latdecomp12(npe,2))
  allocate(londecomp12(npe,2))

  mm1=mype+1

  kk=0
  do n=1,cvec%n3d
     do k=1,nsig
        kk=kk+1
        vname_in_worka(kk)=cvec%r3(n)%shortname
     enddo
  end do

  do n=1,cvec%n2d
     kk=kk+1
     vname_in_worka(kk)=cvec%r2(n)%shortname
  end do

  do n=1,cvec%n3d+cvec%n2d
     if (n <= cvec%n3d) then 
        vorder_in_worka(n)=cvec%r3(n)%shortname
        nvdim(n)=3
      else
        vorder_in_worka(n)=cvec%r2(n)%shortname
        nvdim(n)=2
     endif
  end do
       
  do i=2,lon2-1
     longlobal_index(i-1)=i+jstart(mm1)-2
  enddo

  do j=2,lat2-1
     latglobal_index(j-1)=j+istart(mm1)-2
  end do

  allocate(latdecomp12_aux(npe,2))
  allocate(londecomp12_aux(npe,2))

  latdecomp12_aux=0
  latdecomp12_aux(mm1,1)=istart(mm1)
  latdecomp12_aux(mm1,2)=istart(mm1)+lat1-1
  latdecomp12=0
  call mpi_allreduce(latdecomp12_aux,latdecomp12,2*npe,mpi_itype,mpi_sum,mpi_comm_world,ierror)

  londecomp12_aux=0
  londecomp12_aux(mm1,1)=jstart(mm1)
  londecomp12_aux(mm1,2)=jstart(mm1)+lon1-1
  londecomp12=0
  call mpi_allreduce(londecomp12_aux,londecomp12,2*npe,mpi_itype,mpi_sum,mpi_comm_world,ierror)

  deallocate(latdecomp12_aux)
  deallocate(londecomp12_aux)

end subroutine init_transf_flds

subroutine mk_transf_fields(cvec,iflg)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    mk_transf_fields
!   prgmmr: pondeca          org: np22                date: 2020-09-28
!
! abstract: 
!
! program history log:
!   2020-09-28  pondeca
!   input argument list:
!     t        - t grid values
!     p        - p surface grid values
!     q        - q grid values
!     oz       - ozone grid values
!     skint    - skin temperature grid values
!     cwmr     - cloud water mixing ratio grid values
!     st       - streamfunction grid values
!     vp       - velocity potential grid values
!     sst      - sst grid values
!     slndt    - land surface temperature grid values
!     sicet    - snow/ice covered surface temperature grid values
!
!   output argument list:
!     t        - t grid values
!     p        - p surface grid values
!     q        - q grid values
!     oz       - ozone grid values
!     skint    - skin temperature grid values
!     cwmr     - cloud water mixing ratio grid values
!     st       - streamfunction grid values
!     vp       - velocity potential grid values
!     sst      - sst grid values
!     slndt    - land surface temperature grid values
!     sicet    - snow/ice covered surface temperature grid values
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use constants, only: zero
  use gsi_bundlemod, only : gsi_bundle
  use gridmod, only: nsig,lat2,lon2,nlat,nlon
  use gridmod, only: strip,istart,jstart
  use mpimod, only: mype,ierror,mpi_rtype,mpi_sum,mpi_comm_world
  use gsi_bundlemod, only: gsi_bundlegetpointer
  implicit none

! Declare passed variables
  integer(i_kind),intent(in)  :: iflg
  type(gsi_bundle),intent(in) :: cvec

! Declare local variables
  integer(i_kind) i,j,k,n,kk,mm1,istatus,iglob,jglob
  real(r_kind),pointer,dimension(:,:,:)::ptr3d=>NULL()
  real(r_kind),pointer,dimension(:,:)  ::ptr2d=>NULL()

  mm1=mype+1

  if (iflg == 0 ) then    !adjoint of operator "worka to bundle"
     worka=zero
     kk=0
     do n=1,cvec%n3d
        call gsi_bundlegetpointer ( cvec,cvec%r3(n)%shortname,ptr3d,istatus )
        do k=1,nsig
           kk=kk+1
           do i=2,lon2-1
              do j=2,lat2-1
                 worka(kk,j,i)=worka(kk,j,i)+ptr3d(j,i,k) !clt can we directly use assignment as below? 
              end do
           end do
        end do
     end do

     do n=1,cvec%n2d
        call gsi_bundlegetpointer(cvec,cvec%r2(n)%shortname,ptr2d,istatus)
        kk=kk+1
        do i=2,lon2-1
           do j=2,lat2-1
              worka(kk,j,i)=ptr2d(j,i)
           end do
        end do
     end do
  else                !forward operator "worka to bundle"
!cltthink  how to define halo points 
     kk=0
     do n=1,cvec%n3d
        call gsi_bundlegetpointer ( cvec,cvec%r3(n)%shortname,ptr3d,istatus )
        do k=1,nsig
           kk=kk+1
           do i=2,lon2-1
              do j=2,lat2-1
                 ptr3d(kk,j,i)=worka(kk,j,i)
              end do
           end do
        end do
     end do

     do n=1,cvec%n2d
        call gsi_bundlegetpointer(cvec,cvec%r2(n)%shortname,ptr2d,istatus)
        kk=kk+1
        do i=2,lon2-1
           do j=2,lat2-1
              ptr2d(j,i)=worka(j,i,kk)
           end do
        end do
     end do

  end if

end subroutine mk_transf_fields

subroutine destroy_transf_flds
  deallocate(worka)
  deallocate(vname_in_worka)
  deallocate(vorder_in_worka)
  deallocate(nvdim)
  deallocate(latglobal_index)
  deallocate(longlobal_index)
  deallocate(latdecomp12)
  deallocate(londecomp12)
end subroutine destroy_transf_flds

subroutine init_mg_get_background
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    create_transf_flds 
!   prgmmr: pondeca          org: np22                date: 2021-01-26
!
! abstract: 
!
! program history log:
!   2021-01-25  pondeca
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use gridmod, only: lat2,lon2
  implicit none

! Declare passed variables

! Declare local variables


  allocate(zfld(lat2,lon2))

end subroutine init_mg_get_background

subroutine mg_get_background
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    mg_get_background
!   prgmmr: pondeca          org: np22                date: 2021-01-25
!
! abstract: 
!
! program history log:
!   2020-09-28  pondeca
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use guess_grids, only: ntguessig
  use gsi_bundlemod, only : gsi_bundle
  use gridmod, only: strip
  use gsi_metguess_mod, only: gsi_metguess_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  implicit none

! Declare passed variables

! Declare local variables
  integer(i_kind) istatus,it
  real(r_kind),pointer,dimension(:,:)  ::ptr2d=>NULL()


  it=ntguessig
  call gsi_bundlegetpointer (gsi_metguess_bundle(it),'z' ,ptr2d,   istatus)
  zfld=ptr2d

end subroutine mg_get_background

subroutine destroy_mg_get_background
  deallocate(zfld)
end subroutine destroy_mg_get_background
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine mkgrads4ctrvec(cvec)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    mk_transf_fields
!   prgmmr: pondeca          org: np22                date: 2020-09-28
!
! abstract: 
!
! program history log:
!   2021-04-20  pondeca
!   input argument list:

!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_single,i_kind
  use constants, only: zero_single,max_varname_length
  use gsi_bundlemod, only : gsi_bundle
  use gridmod, only: nsig,lat2,lon2,nlat,nlon
  use gridmod, only: istart,jstart
  use mpimod, only: mype,ierror,mpi_real4,mpi_sum,mpi_comm_world
  use gsi_bundlemod, only: gsi_bundlegetpointer
  implicit none

! Declare passed variables
  type(gsi_bundle),intent(in) :: cvec

! Declare local variables
  integer(i_kind) i,j,k,n,mm1,istatus,iglob,jglob
  real(r_kind),pointer,dimension(:,:,:)::ptr3d=>NULL()
  real(r_kind),pointer,dimension(:,:)  ::ptr2d=>NULL()
  real(r_single),allocatable,dimension(:,:):: outwork,outwork0
  character(max_varname_length) label

! Declare local parmeters
  integer(i_kind),parameter:: lungrds=16

  allocate(outwork(nlon,nlat))
  allocate(outwork0(nlon,nlat))

  mm1=mype+1

  do n=1,cvec%n3d
     call gsi_bundlegetpointer ( cvec,cvec%r3(n)%shortname,ptr3d,istatus )
     label='outgrads_'//trim(cvec%r3(n)%shortname)
     if (mype==0) open (lungrds,file=trim(label)//'.dat',form='unformatted')
     do k=1,nsig
        outwork=zero_single
!        outwork0=zero_single    !T misha
        do j=2,lon2-1
           jglob=jstart(mm1)-2+j
           do i=2,lat2-1
              iglob=istart(mm1)-2+i
              outwork(jglob,iglob)=real(ptr3d(i,j,k),kind=r_single)
           enddo
        enddo
        call mpi_reduce(outwork,outwork0,nlon*nlat,mpi_real4,mpi_sum,0,mpi_comm_world,ierror)
        if(mype==0) write(lungrds) outwork0
     enddo
     if (mype==0) then
         close(lungrds)
         call mkgrdsctrl(nlon,nlat,nsig,label)
     endif
  enddo
  call mpi_barrier(mpi_comm_world,ierror)

  do n=1,cvec%n2d
     call gsi_bundlegetpointer(cvec,cvec%r2(n)%shortname,ptr2d,istatus)
     label='outgrads_'//trim(cvec%r2(n)%shortname)
     if (mype==0) open (lungrds,file=trim(label)//'.dat',form='unformatted')
     outwork=zero_single
!     outwork0=zero_single     !T Misha
     do j=2,lon2-1
        jglob=jstart(mm1)-2+j
        do i=2,lat2-1
           iglob=istart(mm1)-2+i
           outwork(jglob,iglob)=real(ptr2d(i,j),kind=r_single)
        enddo
     enddo
     call mpi_reduce(outwork,outwork0,nlon*nlat,mpi_real4,mpi_sum,0,mpi_comm_world,ierror)
     if(mype==0) then
        write(lungrds) outwork0
        close(lungrds)
        call mkgrdsctrl(nlon,nlat,1,label)
     end if
  enddo
  call mpi_barrier(mpi_comm_world,ierror)


  deallocate(outwork)
  deallocate(outwork0)

end subroutine  mkgrads4ctrvec
!----------------------------------------------------------------------------------
subroutine mkgrdsctrl(nx,ny,np,label)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    mkgrdsctrl2d
!   prgmmr:
!
! abstract:
!
! program history log:
!   2021-02-22  pondeca - adapted from parrish outgrads1
!
!   input argument list:
!    label
!    nx,ny
!    f
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  use kinds, only: i_kind,r_single
  implicit none

  character(*)   ,intent(in   ) :: label
  integer(i_kind),intent(in   ) :: nx,ny,np

  integer(i_kind) i,l,next,last,ntime,ioutcor,koutmax
  real(r_single) rlonmap0,undef,dlonmap,pinc,startp,rlatmap0,dlatmap
  character(80) dsdes,dsdat
  character(80) datdes(1000)
  character(1) blank
  data blank/' '/
  data undef/-9.99e33_r_single/

  ioutcor=10

  write(dsdes,'(a,".des")')trim(label)
  write(dsdat,'(a,".dat")')trim(label)
  open(unit=ioutcor,file=dsdes,form='formatted')
  ntime=1
  rlonmap0=1._r_single
  dlonmap=1._r_single
  rlatmap0=1._r_single
  dlatmap=1._r_single
  startp=1._r_single
  pinc=1._r_single
  koutmax=1
  do i=1,1000
     write(datdes(i),'(80a1)')(blank,l=1,80)
  end do
  write(datdes(1),'("DSET ",a)')trim(dsdat)
  write(datdes(2),'("options big_endian sequential")')
  write(datdes(3),'("TITLE ",a)')trim(label)
  write(datdes(4),'("UNDEF ",e11.2)')undef
  write(datdes(5),'("XDEF ",i5," LINEAR ",f7.2,f7.2)')nx,rlonmap0,dlonmap
  write(datdes(6),'("YDEF ",i5," LINEAR ",f7.2,f7.2)')ny,rlatmap0,dlatmap
  next=7
  write(datdes(next),'("ZDEF ",i5," LINEAR ",f7.2,f7.2)')np,startp,pinc
  next=next+1
  write(datdes(next),'("TDEF ",i5," LINEAR 0Z23may1992 24hr")')koutmax
  next=next+1
  write(datdes(next),'("VARS 1")')
  next=next+1
  write(datdes(next),'("f   ",i5," 99 f   ")')np
  next=next+1
  write(datdes(next),'("ENDVARS")')
  last=next
  write(ioutcor,'(a80)')(datdes(i),i=1,last)
  close(ioutcor)

end subroutine mkgrdsctrl
end module pre_multigrid
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
