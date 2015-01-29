subroutine read_amsr2(mype,val_amsr2,ithin,isfcalc,rmesh,gstime,&
     infile,lunout,obstype,nread,ndata,nodata,twind,sis,&
     mype_root,mype_sub,npe_sub,mpi_comm_sub)

! subprogram:    read_amsr2                  read bufr format amsr2 data
! prgmmr :   e.jones         org: copied from read_amsre.f90   date: 2014-10-17
!
! abstract:  This routine reads BUFR format AMSR2 radiance (brightness
!            temperature) files.  Optionally, the data are thinned to 
!            a specified resolution using simple quality control checks.
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! Special Notes:
!     10/14/04  okamoto  looks like AMSRE overlap problem is not as bad as SSM/I
!
! program history log:
!   2014-04-25  ejones - copy read_amsre.f90, and modify to ingest amsr2 data.
!                  obstype=amsr2
!   2014-06-20  ejones - modify so amsr2 channels are read in all at once
!                  (single bufr file),getting rid of the need for classifications 
!                  for LR and HR channels (from separate bufr files).
!                  obstype=amsr2
!
! input argument list:
!     mype     - mpi task id
!     val_amsr2- weighting factor applied to super obs
!     ithin    - flag to thin data
!     isfcalc  - flag to specify method to calculate sfc fields within FOV
!                when set to one, account for size/shape of FOV.  otherwise
!                use bilinear interpolation.
!     rmesh    - thinning mesh size (km)
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!     mype_root - "root" task for sub-communicator
!     mype_sub - mpi task id within sub-communicator
!     npe_sub  - number of data read tasks
!     mpi_comm_sub - sub-communicator for data read
!
! output argument list:
!     nread    - number of BUFR AMSR2 observations read
!     ndata    - number of BUFR AMSR2 profiles retained for further processing
!     nodata   - number of BUFR AMSR2 observations retained for further processing
!
! attributes:
!     language: f90
!
!$$$
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,map2tgrid,destroygrids, &
      checkob,finalcheck,score_crit
  use radinfo, only: iuse_rad,nusis,jpch_rad,nst_gsi,nstinfo 
  use gridmod, only: diagnostic_reg,regional,nlat,nlon,rlats,rlons,&
      tll2xy
  use constants, only: deg2rad,rad2deg,zero,one,three,r60inv,two
  use gsi_4dvar, only: l4dvar, iwinbgn, winlen
  use calc_fov_conical, only: instrument_init
  use deter_sfc_mod, only: deter_sfc_fov,deter_sfc,deter_sfc_amsre_low
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth, gsi_nstcoupler_deter

  implicit none

! Input variables
  character(len=*) ,intent(in   ) :: infile
  character(len=*) ,intent(in   ) :: obstype
  integer(i_kind)  ,intent(in   ) :: mype
  integer(i_kind)  ,intent(inout) :: isfcalc
  integer(i_kind)  ,intent(in   ) :: ithin
  integer(i_kind)  ,intent(in   ) :: lunout
  real(r_kind)     ,intent(inout) :: val_amsr2
  real(r_kind)     ,intent(in   ) :: gstime,twind
  real(r_kind)     ,intent(in   ) :: rmesh
  character(len=20),intent(in   ) :: sis
  integer(i_kind)  ,intent(in   ) :: mype_root
  integer(i_kind)  ,intent(in   ) :: mype_sub
  integer(i_kind)  ,intent(in   ) :: npe_sub
  integer(i_kind)  ,intent(in   ) :: mpi_comm_sub

! Output variables
  integer(i_kind)  ,intent(inout) :: nread
  integer(i_kind)  ,intent(inout) :: ndata,nodata

! Number of channels for sensors in BUFR
  integer(i_kind),parameter :: N_AMSR2CH  =  14
  integer(i_kind) :: said, GCOMW1_SAID  = 122   !WMO satellite identifier 
  integer(i_kind) :: siid, AMSR2_SIID = 478   !WMO instrument identifier 
  integer(i_kind),parameter :: maxinfo    =  33

! BUFR file sequencial number
  character(len=8)  :: subset
  character(len=4)  :: senname
  integer(i_kind)   :: lnbufr = 10
  integer(i_kind)   :: nchanl
  integer(i_kind)   :: iret,isflg,idomsfc(1)

! Work variables for time
  integer(i_kind)   :: idate
  integer(i_kind)   :: idate5(5)
  integer(i_kind)   :: nmind
  real(r_kind)      :: sstime, tdiff, t4dv

! Other work variables
  logical           :: outside,iuse,assim
  integer(i_kind)   :: nreal, kidsat
  integer(i_kind)   :: itx, k, nele, itt
  integer(i_kind)   :: ifov, ilat, ilon
  integer(i_kind)   :: i, l, n
  integer(i_kind),dimension(n_amsr2ch) :: kchamsr2
  real(r_kind)     :: sfcr
  real(r_kind)     :: dlon, dlat
  real(r_kind)     :: dlon_earth,dlat_earth
  real(r_kind)     :: timedif, pred, crit1, dist1
  real(r_kind),allocatable,dimension(:,:):: data_all
  integer(i_kind),allocatable,dimension(:)::nrec
  integer(i_kind):: irec,isub,next
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind),dimension(0:3):: ts
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr

  character(len=7),parameter:: fov_flag="conical"
  character(len=3) :: fov_satid
  
  integer(i_kind) :: ichan, instr, idum

  logical :: valid

  real(r_kind) :: clath_sun_glint_calc , clonh_sun_glint_calc 
  real(r_kind) :: date5_4_sun_glint_calc
  real(r_kind) :: expansion, dlat_earth_deg, dlon_earth_deg

! Set standard parameters
!  logical       :: amsr2_low
!  logical       :: amsr2_mid
!  logical       :: amsr2_hig
  integer(i_kind) ntest
  integer(i_kind) :: nscan,iskip,kskip,kch,kchanl
! real(r_kind),parameter :: TEN      =  10._r_kind
! real(r_kind),parameter :: R45      =  45._r_kind
  real(r_kind),parameter :: R90      =  90._r_kind
! real(r_kind),parameter :: R180     = 180._r_kind
  real(r_kind),parameter :: R360     = 360._r_kind
  real(r_kind),parameter :: tbmin    = 3._r_kind           
  real(r_kind),parameter :: tbmax    = 340._r_kind         
  real(r_kind) disterrmax
  real(r_kind),dimension(N_AMSR2CH) :: tbob_org
  real(r_kind) :: clath, clonh, fovn, iang, aang, soel, solazi         !????
  real(r_kind) :: flgch  !used for thinning priority  range:1-36


! BUFR format for AMSR2SPOT
  integer(i_kind),parameter :: N_AMSR2SPOT_LIST = 12

! BUFR format for AMSR2CHAN
  integer(i_kind),parameter :: N_AMSR2CHAN_LIST = 3

! Variables for BUFR IO
  real(r_double),dimension(4):: gcomspot_d
  real(r_double),dimension(12):: amsr2spot_d                !????
  real(r_double),dimension(3,14):: amsr2chan_d              !????

! ---- sun glint ----
  integer(i_kind) doy,mlen(12),mday(12),mon,m
  real(r_kind) sun_elevation,sun_azimuth,bearaz,sun_zenith
  data  mlen/31,28,31,30,31,30, &
             31,31,30,31,30,31/ 

! Orbit
  integer(i_kind) :: orbit, old_orbit, iorbit, ireadsb, ireadmg
  real(r_kind) :: saz


  real(r_kind),parameter:: one_minute=0.01666667_r_kind
  real(r_kind),parameter:: minus_one_minute=-0.01666667_r_kind


! ---- erin's debug shtuff ----
integer(i_kind) :: land_ob, sea_ob, seaice_ob, snow_ob, mixed_ob
integer(i_kind) :: iuse1_yes, iuse1_no, iuse2_yes, iuse2_no, iuse3_yes, iuse3_no
integer(i_kind) :: weird_soel, weird_time, weird_tdiff, weird_latlon,kskipped, kskipped2


! data selection
! Initialize variables
! ---------- erin debug stuff -------
sea_ob=0
land_ob=0
seaice_ob=0
snow_ob=0
mixed_ob=0

iuse1_yes=0
iuse1_no=0
iuse2_yes=0
iuse2_no=0
iuse3_yes=0
iuse3_no=0

weird_soel=0
weird_time=0
weird_tdiff=0
weird_latlon=0
kskipped=0
kskipped2=0
!----------- / erin debug stuff ---------

  ilon = 3
  ilat = 4

  if (nst_gsi > 0 ) then
     call gsi_nstcoupler_skindepth(obstype, zob)         ! get penetration depth (zob) for the obstype
  endif

  m = 0
  do mon=1,12 
     mday(mon) = m 
     m = m + mlen(mon) 
  end do 
  disterrmax=zero
  ntest = 0
  nreal = maxinfo+nstinfo
  ndata = 0
  nodata = 0
  orbit = -1
  old_orbit=-1
  iorbit = 0
  sstime = zero
  kchanl=14
  kchamsr2(1:14)=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14/)

  senname = 'AMSR'
  nchanl  = N_AMSR2CH
  nscan  = 243  !for the lower frequency channels
!     nscan  = 486  !for amsr2 89.0GHz ch         
  kidsat = 549
  rlndsea(0) = zero
  rlndsea(1) = 15._r_kind
  rlndsea(2) = 10._r_kind
  rlndsea(3) = 15._r_kind
  rlndsea(4) = 100._r_kind

! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  assim=.false.
  search: do i=1,jpch_rad
     if ((nusis(i)==sis) .and. (iuse_rad(i)>0)) then
        assim=.true.
        exit search
     endif
  end do search
  if (.not.assim) val_amsr2=zero

! note, fov-based surface code does not have equations for amsr-e or amsr2, but
! the fov size/shape is similar to ssmi/s.  so call
! fov code using f16 specs.

  if(isfcalc==1)then
     instr=26
     fov_satid='f16'
     idum = -999  ! dummy variable for fov number. not used for conical instr.
     ichan=14
     expansion=2.9_r_kind
     call instrument_init(instr, fov_satid, expansion, valid)
     if (.not. valid) then
       if (assim) then
         write(6,*)'READ_AMSR2:  ***ERROR*** IN SETUP OF FOV-SFC CODE. STOP'
         call stop2(71)
       else
         isfcalc = 0
         write(6,*)'READ_AMSR2:  ***ERROR*** IN SETUP OF FOV-SFC CODE'
       endif
    endif
  endif

! Make thinning grids
  call makegrids(rmesh,ithin)


! Open BUFR file
  open(lnbufr,file=infile,form='unformatted')
  call openbf(lnbufr,'IN',lnbufr)
  call datelen(10)

!write(6,*)'amsr2 check3'

! Allocate local array to contain observation information
  nele=nreal+nchanl
write(6,*)'AMSR2 NREAL=',nreal
  allocate(data_all(nele,itxmax),nrec(itxmax))

! Big loop to read data file
  next=0
  nrec=999999
  irec=0
  do while(ireadmg(lnbufr,subset,idate)>=0)
     irec=irec+1
     next=next+1
     if(next == npe_sub)next=0
     if(next /= mype_sub)cycle
     read_loop: do while (ireadsb(lnbufr)==0)

!    Retrieve bufr 1/4 :get gcomspot (said,orbn,soza)
        call ufbint(lnbufr,gcomspot_d,4,1,iret,'SAID ORBN SOLAZI SOEL')

        said = nint(gcomspot_d(1))

!       Retrieve bufr 2/4 :get amsr2spot (siid,ymdhs,lat,lon)
        call ufbrep(lnbufr,amsr2spot_d,N_AMSR2SPOT_LIST,1,iret, &
           'SIID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH AANG IANG FOVN')

        siid = nint(amsr2spot_d(1)) 
        if(siid /= AMSR2_SIID)   cycle read_loop


!       Check obs time
        idate5(1) = amsr2spot_d(02)! year
        idate5(2) = amsr2spot_d(03)! month
        idate5(3) = amsr2spot_d(04)! day
        idate5(4) = amsr2spot_d(05)! hour
        idate5(5) = amsr2spot_d(06)! min
        if( idate5(1) < 1900 .or. idate5(1) > 3000 .or. &
            idate5(2) < 1    .or. idate5(2) >   12 .or. &
            idate5(3) < 1    .or. idate5(3) >   31 .or. &
            idate5(4) < 0    .or. idate5(4) >   24 .or. &
            idate5(5) < 0    .or. idate5(5) >   60 )then
            write(6,*)'READ_AMSR2:  ### ERROR IN READING BUFR DATA:', &
              ' STRANGE OBS TIME (YMDHM):', idate5(1:5)
! ---- erin debug stuff ----
! weird_time = weird_time+1
! ---- /erin debug stuff ----
            cycle read_loop     ! remove filter for debugging
        endif


        call w3fs21(idate5,nmind)
        t4dv = (real((nmind-iwinbgn),r_kind) + amsr2spot_d(7)*r60inv)*r60inv ! add in seconds

        if (l4dvar) then
           if (t4dv<zero .OR. t4dv>winlen) cycle read_loop
        else
           sstime = real(nmind,r_kind) + amsr2spot_d(7)*r60inv ! add in seconds
           tdiff  = (sstime - gstime)*r60inv
           if (abs(tdiff)>twind) then !cycle read_loop !then 
! ---- erin debug stuff ----
! weird_tdiff = weird_tdiff+1
! ---- /erin debug stuff ----       
 cycle read_loop              ! remove filter for debugging
endif
         endif
        if (l4dvar) then
           timedif = zero
        else
           timedif = two*abs(tdiff)        ! range:  0 to 6     - used in atms          
                                           !6.0_r_kind*abs(tdiff) ! range:  0 to 18     -original
        endif


!     --- Check observing position -----
        clath= amsr2spot_d(08)
        clonh= amsr2spot_d(09)

        if( abs(clath) > R90  .or. abs(clonh) > R360 .or. &
          ( abs(clath) == R90 .and. clonh /= ZERO) )  then
           write(6,*)'READ_AMSR2:  ### ERROR IN READING BUFR DATA:',&
              ' STRANGE OBS POINT (LAT,LON):', clath, clonh
              cycle read_loop
        endif

! ---- erin debug stuff ----
! weird_latlon = weird_latlon+1
! ---- /erin debug stuff ----
!              cycle read_loop        ! remove filter for debugging
!            endif 

!    Pick up every three scene  3,6,9,,,,195 (num=65)
!    because of too small scene size and too many scene numbers
!    (low-freq ch FOV are overlapped side by side)
        fovn = amsr2spot_d(12)

     
!    Set position in a given region
        if(clonh >= R360)then
           clonh = clonh - R360
        else if(clonh < ZERO)then
           clonh = clonh + R360
        endif
     
!    If regional, map obs lat,lon to rotated grid.
        dlat_earth = clath * deg2rad
        dlon_earth = clonh * deg2rad
        dlat_earth_deg = clath
        dlon_earth_deg = clonh

        if(regional)then
!       Convert to rotated coordinate.  dlon centered on 180 (pi),
!       so always positive for limited area
           call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)

!       Check to see if in domain.  outside=.true. if dlon_earth,
!       dlat_earth outside domain, =.false. if inside
           if(outside) cycle read_loop

!    Global case 
        else
           dlat=dlat_earth
           dlon=dlon_earth
           call grdcrd1(dlat,rlats,nlat,1)
           call grdcrd1(dlon,rlons,nlon,1)
        endif

!    Sum number of read obs before thinning step.  Note that this number will contain
!    some observations that may be rejected later due to bad BTs.
        nread=nread+kchanl
    
        crit1 = 0.01_r_kind+timedif

!----erin debug stuff----
!write(6,*)'AMSR2_CRIT1 #1-',crit1
!-----/erin debug stuff----

        call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis)
! ---- erin debug stuff ----
!if(iuse)then
!  iuse1_yes = iuse1_yes+1
!else if (.not. iuse) then
!  iuse1_no = iuse1_no+1
!  cycle read_loop                   ! remove filtering for debugging
!endif
! ---- /erin debug stuff ----
        if (.not.iuse) cycle read_loop


!    QC:  "Score" observation.  We use this information to identify "best" obs

!       Locate the observation on the analysis grid.  Get sst and land/sea/ice
!       mask   

!     isflg    - surface flag
!                0 sea
!                1 land
!                2 sea ice
!                3 snow
!                4 mixed                       

!       When isfcalc is one, calculate surface fields based on size/shape of fov.
!       Otherwise, use bilinear interpolation.

        if (isfcalc==1)then
           call deter_sfc_fov(fov_flag,idum,instr,ichan,real(amsr2spot_d(11),r_kind),dlat_earth_deg,&
              dlon_earth_deg,expansion,t4dv,isflg,idomsfc(1), &
              sfcpct,vfr,sty,vty,stp,sm,ff10,sfcr,zz,sn,ts,tsavg)
        else
           call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc(1),sfcpct, &
              ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)

! ---- erin debug stuff ----
!if(isflg==0)then
!  sea_ob = sea_ob+1
!else if (isflg==1)then
!  land_ob = land_ob+1
!else if (isflg==2)then
!  seaice_ob = seaice_ob+1
!else if (isflg==3)then
!  snow_ob = snow_ob+1
!else if (isflg==4)then
!  mixed_ob = mixed_ob+1
!endif
! ---- /erin debug stuff ----

         ! move this to QC - 
           if(isflg/=0) cycle read_loop                            ! use data over water only
        endif


!---- erin debug stuff ----
!write(6,*)'AMSR2_CRIT1 #2-',crit1
!-----/erin debug stuff----
        crit1 = crit1 +rlndsea(isflg)
!----erin debub stuff----
!write(6,*)'AMSR2_CRIT1 #3-',crit1
!----/erin debug stuff----
!write(6,*)'AMSR2_DIST1',dist1

        call checkob(dist1,crit1,itx,iuse)
        if(.not. iuse)cycle read_loop

! ---- erin debug stuff ----
!if(iuse)then
!  iuse2_yes = iuse2_yes+1
!else if (.not. iuse) then
!  iuse2_no = iuse2_no+1
!  cycle read_loop          ! remove filtering for debugging
!endif
! ---- /erin debug stuff ----

!       Retrieve bufr 3/4 : get amsr2chan (chnm,tbb)
        call ufbrep(lnbufr,amsr2chan_d,3,14,iret,'SCCF ACQF TMBR')   

        tbob_org(1)=amsr2chan_d(3,1)         
        tbob_org(2)=amsr2chan_d(3,2)
        tbob_org(3)=amsr2chan_d(3,3)
        tbob_org(4)=amsr2chan_d(3,4)
        tbob_org(5)=amsr2chan_d(3,5)
        tbob_org(6)=amsr2chan_d(3,6)
        tbob_org(7)=amsr2chan_d(3,7)
        tbob_org(8)=amsr2chan_d(3,8)
        tbob_org(9)=amsr2chan_d(3,9)
        tbob_org(10)=amsr2chan_d(3,10)
        tbob_org(11)=amsr2chan_d(3,11)
        tbob_org(12)=amsr2chan_d(3,12)
        tbob_org(13)=amsr2chan_d(3,13)
        tbob_org(14)=amsr2chan_d(3,14)


!       Set obs information

        iskip = 0
        do l=1,nchanl
           if(tbob_org(l)<tbmin .or. tbob_org(l)>tbmax)then
              iskip = iskip + 1
           end if
        end do
        kskip = 0
        do l=1,kchanl
           kch=kchamsr2(l)
           if(tbob_org(kch)<tbmin .or. tbob_org(kch)>tbmax)then
              kskip = kskip + 1
! ---- erin debug stuff ----
! kskipped = kskipped +1
! ----/erin debug stuff ----
           endif
        end do
        if(kskip == kchanl .or. iskip == nchanl) then
!---- erin debug stuff-----
!kskipped2 = kskipped2+1
!-----/erin debug stuff---- 
           cycle read_loop                ! remove filtering for debugging
endif

! ** INVESIGATE THIS flgch=iskip*three
        flgch=iskip*three  !used for thin, range 0 to 36
        crit1 = crit1 + flgch

!    Set data quality predictor ***NEED TO COME UP WITH THIS***
        pred = zero


!    Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"

!    Map obs to grids
        crit1 = crit1+pred
        call finalcheck(dist1,crit1,itx,iuse)
        if(.not. iuse)cycle read_loop

! ---- erin debug stuff ----
!if(iuse)then
!  iuse3_yes = iuse3_yes+1
!else if (.not. iuse) then
!  iuse3_no = iuse3_no+1
!  cycle read_loop          ! remove filtering for debugging
!endif
! ---- /erin debug stuff ----
     
        soel = gcomspot_d(4)

!    Check observational info 

        if( soel < -180._r_kind .or. soel > 180._r_kind )then
           write(6,*)'READ_AMSR2:  ### ERROR IN READING BUFR DATA:', &
              ' STRANGE OBS INFO(FOV,SOLAZI,SOEL):', fovn, solazi, soel
! ---- erin debug stuff ----
!            weird_soel=weird_soel+1
! ---- /erin debug stuff ----
           cycle read_loop                  ! remove filtering for debugging
        endif

!  -------- Retreive Sun glint angle -----------
!----- erin stuff - modify/get rid of this for QC?
        doy = mday( int(idate5(2)) ) + int(idate5(3))
        if ((mod( int(idate5(1)),4)==0).and.( int(idate5(2)) > 2))  then 
           doy = doy + 1
        end if 

        ifov = nint(fovn)


        clath_sun_glint_calc = clath
        clonh_sun_glint_calc = clonh
        if(clonh>180_r_kind) clonh_sun_glint_calc = clonh -360.0_r_kind
!       date5_4_sun_glint_calc = idate5(4)  
        date5_4_sun_glint_calc =  &                                                                                                
        real(idate5(4),r_kind)+real(idate5(5),r_kind)*r60inv+real(amsr2spot_d(7),r_kind)*r60inv*r60inv   


!write(6,*)'amsr2 check6'
!----- erin stuff - is call zensun needed??      
        call zensun(doy,date5_4_sun_glint_calc,clath_sun_glint_calc,clonh_sun_glint_calc,sun_zenith,sun_azimuth)

        saz = amsr2spot_d(11)*deg2rad    ! satellite zenith/incidence angle (rad) ???

!          saz = amsr2spot(10,1)*deg2rad   ! satellite zenith angle (rad) 
!             ==> not use this value but fixed values(55.0 deg)   10/12/04
!             because BUFR saza value looks strange (raging -3 to 25),
!
!       interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
!
        if ( nst_gsi > 0 ) then
           tref  = ts(0)
           dtw   = zero
           dtc   = zero
           tz_tr = one
           if ( sfcpct(0) > zero ) then
             call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr) 
           endif
        endif

        data_all(1,itx) = 122.00                     ! satellite ID
        data_all(2,itx) = t4dv                       ! time diff (obs - anal) (hours)
        data_all(3,itx) = dlon                       ! grid relative longitude
        data_all(4,itx) = dlat                       ! grid relative latitude
        data_all(5,itx) = saz                        ! satellite zenith angle (rad)
        data_all(6,itx) = amsr2spot_d(10)             ! satellite azimuth angle
        data_all(7,itx) = zero                       ! look angle (rad)
        data_all(8,itx) = ifov                       ! fov number    1-196
!        data_all(8,itx) = ifov/3 + 1                 ! fov number/3  1-65 !kozo
        data_all(9,itx) = sun_zenith                 ! solar zenith angle (deg)
        data_all(10,itx) = sun_azimuth                ! solar azimuth angle (deg)
        data_all(11,itx) = sfcpct(0)                 ! sea percentage of
        data_all(12,itx) = sfcpct(1)                 ! land percentage
        data_all(13,itx) = sfcpct(2)                 ! sea ice percentage
        data_all(14,itx) = sfcpct(3)                 ! snow percentage
        data_all(15,itx) = ts(0)                      ! ocean skin temperature
        data_all(16,itx) = ts(1)                      ! land skin temperature
        data_all(17,itx) = ts(2)                      ! ice skin temperature
        data_all(18,itx) = ts(3)                      ! snow skin temperature
        data_all(19,itx) = tsavg                      ! average skin temperature
        data_all(20,itx) = vty                        ! vegetation type
        data_all(21,itx) = vfr                        ! vegetation fraction
        data_all(22,itx) = sty                        ! soil type
        data_all(23,itx) = stp                        ! soil temperature
        data_all(24,itx) = sm                         ! soil moisture
        data_all(25,itx) = sn                         ! snow depth
        data_all(26,itx) = zz                         ! surface height
        data_all(27,itx) = idomsfc(1) + 0.001_r_kind  ! dominate surface type
        data_all(28,itx) = sfcr                       ! surface roughness
        data_all(29,itx) = ff10                       ! ten meter wind factor
        data_all(30,itx) = dlon_earth*rad2deg         ! earth relative longitude (degrees)
        data_all(31,itx) = dlat_earth*rad2deg         ! earth relative latitude (degrees)

        data_all(32,itx)= val_amsr2
        data_all(33,itx)= itt

        if ( nst_gsi > 0 ) then
           data_all(maxinfo+1,itx) = tref                ! foundation temperature
           data_all(maxinfo+2,itx) = dtw                 ! dt_warm at zob
           data_all(maxinfo+3,itx) = dtc                 ! dt_cool at zob
           data_all(maxinfo+4,itx) = tz_tr               ! d(Tz)/d(Tr)
        endif

        do l=1,nchanl
           data_all(l+nreal,itx) = tbob_org(l)
        end do
        nrec(itx)=irec

     enddo read_loop
  enddo
  call closbf(lnbufr)

! If multiple tasks read input bufr file, allow each tasks to write out
! information it retained and then let single task merge files together


! ------------------- erin debug stuff ------------
!write(6,*)'**SFC TYPE STUFF** land_ob, sea_ob, seaice_ob, snow_ob, mixed_ob='&
!          ,land_ob, sea_ob, seaice_ob, snow_ob, mixed_ob
! write(6,*)'**WEIRD STUFF** weird_tdiff', weird_tdiff
!write(6,*)'**WEIRD STUFF** weird_latlon', weird_latlon
!write(6,*)'**IUSE STUFF** iuse1_yes, iuse1_no, iuse2_yes, iuse2_no, iuse3_yes, iuse3_no=',iuse1_yes, iuse1_no, iuse2_yes, iuse2_no, iuse3_yes, iuse3_no

!write(6,*)'**WEIRD STUFF** weird_soel, weird_time, weird_tdiff, weird_latlon, kskipped=',weird_soel,weird_time, weird_tdiff, weird_latlon, kskipped

!write(6,*)'**PRE COMBINERADOBS** mype, nread, ndata=', mype, nread, ndata
!----------------- /erin debug stuff -------------

  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
     nele,itxmax,nread,ndata,data_all,score_crit,nrec)

!write(6,*) '**POST COMBINERADOBS** nread,ndata=', nread, ndata


! Allow single task to check for bad obs, update superobs sum,
! and write out data to scratch file for further processing.
  if (mype_sub==mype_root.and.ndata>0) then

!    Identify "bad" observation (unreasonable brightness temperatures).
!    Update superobs sum according to observation location
     do n=1,ndata
        do i=1,nchanl
!           if(data_all(i+nreal,n) > tbmin .and. &         !remove gross check for debugging
!              data_all(i+nreal,n) < tbmax)nodata=nodata+1
            nodata=nodata+1                                ! each ob goes through     
        end do
        itt=nint(data_all(maxinfo,n))
        super_val(itt)=super_val(itt)+val_amsr2

     end do

!    Write final set of "best" observations to output file
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
  endif

  deallocate(data_all,nrec) ! Deallocate data arrays
  call destroygrids    ! Deallocate satthin arrays
  if(diagnostic_reg.and.ntest>0 .and. mype_sub==mype_root) &
     write(6,*)'READ_AMSR2:  ',&
        'mype,ntest,disterrmax=',mype,ntest,disterrmax

  return


end subroutine read_amsr2

