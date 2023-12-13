!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module  mg_timers
!***********************************************************************
!                                                                      !
!  Measure cpu and wallclock timing                                    !
!                                                     D. Jovic  (2017) !
!                                        Adjusted:    M. Rancic (2020) !
!***********************************************************************
  use kinds, only: r_kind
  use berror, only: multigrid_betafct
  implicit none

  private
  public :: btim, etim, print_mg_timers

  type timer
    logical :: running = .false.
    real(r_kind) :: start_clock = 0.0
    real(r_kind) :: start_cpu = 0.0
    real(r_kind) :: time_clock = 0.0
    real(r_kind) :: time_cpu = 0.0
  end type timer

  type(timer),save,public ::      total_tim
  type(timer),save,public ::      input_tim
  type(timer),save,public ::     output_tim
  type(timer),save,public ::       init_tim
  type(timer),save,public ::      part1_tim
  type(timer),save,public ::      part2_tim
  type(timer),save,public ::      gsub1_tim
  type(timer),save,public ::      gsub2_tim
  type(timer),save,public ::      gsub3_tim
  type(timer),save,public ::      gsub4_tim
  type(timer),save,public ::      obsr1_tim
  type(timer),save,public ::      obsr2_tim
  type(timer),save,public ::      cldan_tim
  type(timer),save,public ::      guess_tim
  type(timer),save,public ::   dynamics_tim
  type(timer),save,public ::    upsend_tim
  type(timer),save,public ::    upsend1_tim
  type(timer),save,public ::    upsend2_tim
  type(timer),save,public ::    upsend3_tim
  type(timer),save,public ::    an2filt_tim
  type(timer),save,public ::    filt2an_tim
  type(timer),save,public ::    weight_tim
  type(timer),save,public ::    bfiltT_tim
  type(timer),save,public ::      vadv1_tim
  type(timer),save,public ::      bfilt_tim
  type(timer),save,public ::      mgbft_tim
  type(timer),save,public ::      rfilt_tim
  type(timer),save,public ::      binit_tim
  type(timer),save,public ::      bfint_tim
  type(timer),save,public ::     btrns1_tim
  type(timer),save,public ::     btrns2_tim
  type(timer),save,public ::       adv2_tim
  type(timer),save,public ::       vtoa_tim
  type(timer),save,public ::    dnsend_tim
  type(timer),save,public ::    dnsend1_tim
  type(timer),save,public ::    dnsend2_tim
  type(timer),save,public ::    dnsend3_tim
  type(timer),save,public ::     update_tim
  type(timer),save,public ::    physics_tim
  type(timer),save,public ::  radiation_tim
  type(timer),save,public :: convection_tim
  type(timer),save,public :: turbulence_tim
  type(timer),save,public ::  microphys_tim
  type(timer),save,public ::       pack_tim
  type(timer),save,public ::       arrn_tim
  type(timer),save,public ::      aintp_tim
  type(timer),save,public ::       intp_tim
  type(timer),save,public ::       boco_tim
  type(timer),save,public ::      graph_tim

  integer, parameter, public :: print_clock = 1,                        &
                                print_cpu   = 2,                        &
                                print_clock_pct = 3,                    &
                                print_cpu_pct   = 4

contains

!-----------------------------------------------------------------------
  subroutine btim(t)
    implicit none
    type(timer), intent(inout) :: t

    if (t%running) then
      write(0,*)'btim: timer is already running'
      STOP
    end if
    t%running = .true.

    t%start_clock = wtime()
    t%start_cpu = ctime()

  endsubroutine btim
!-----------------------------------------------------------------------
  subroutine etim(t)
    implicit none
    type(timer), intent(inout) :: t
    real(r_kind) :: wt, ct

    wt = wtime()
    ct = ctime()

    if (.not.t%running) then
      write(0,*)'etim: timer is not running'
      STOP
    end if
    t%running = .false.

    t%time_clock = t%time_clock + (wt - t%start_clock)
    t%time_cpu = t%time_cpu + (ct - t%start_cpu)
    t%start_clock = 0.0
    t%start_cpu = 0.0

  endsubroutine etim
!-----------------------------------------------------------------------
  subroutine print_mg_timers(filename, print_type)
    use mpi
    use mpimod, only: mype
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: print_type

    integer :: fh
    integer :: ierr
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(MPI_STATUS_SIZE) :: stat
    character(len=1024) :: buffer, header
    integer :: bufsize

    call MPI_File_open(MPI_COMM_WORLD, filename, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, fh, ierr)

    buffer = ' '
    if ( print_type == print_clock ) then
!    write(buffer,"(I6,12(',',F10.4))") mype,                            &
    write(buffer,"(I6,2(',',F10.4))") mype,                            &
!                                       init_tim%time_clock,             &
!                                       upsend_tim%time_clock,           &
!                                       dnsend_tim%time_clock,           &
!                                       weight_tim%time_clock,           &
!                                       bfiltT_tim%time_clock,           &
                                       bfilt_tim%time_clock,            &
!                                       filt2an_tim%time_clock,          &
!                                       aintp_tim%time_clock,            &
!                                       intp_tim%time_clock,             &
!                                       an2filt_tim%time_clock,          &
!                                       output_tim%time_clock,           &
                                       total_tim%time_clock
    else if ( print_type == print_cpu ) then
      if(multigrid_betafct) then
    write(buffer,"(I6,9(',',F10.4))") mype,                            &
!                                       init_tim%time_cpu,               &
                                       binit_tim%time_cpu,               &
!                                       upsend_tim%time_cpu,             &
!                                       dnsend_tim%time_cpu,             &
                                       input_tim%time_cpu,             &
                                      output_tim%time_cpu,             &
!                                       dnsend1_tim%time_cpu,             &
!                                       dnsend2_tim%time_cpu,             &
!                                       dnsend3_tim%time_cpu,             &
                                       weight_tim%time_cpu,             &
                                       bfiltT_tim%time_cpu,             &
                                       bfilt_tim%time_cpu,              &
                                       mgbft_tim%time_cpu,               &
!                                       btrns1_tim%time_cpu,              &
!                                       btrns2_tim%time_cpu,              &
!                                       pack_tim%time_cpu,              &
!                                       filt2an_tim%time_cpu,            &
!                                       aintp_tim%time_cpu,              &
!                                       intp_tim%time_cpu,               &
!                                       an2filt_tim%time_cpu,            &
!                                       output_tim%time_cpu,             &
                                       bfint_tim%time_cpu,               &
!                                       graph_tim%time_cpu,               &
                                       total_tim%time_cpu
      else
    write(buffer,"(I6,4(',',F10.4))") mype,                            &
!                                        init_tim%time_cpu,               &
!                                        part1_tim%time_cpu,               &
!                                        part2_tim%time_cpu,               &
!                                        gsub1_tim%time_cpu,               &
!                                        obsr1_tim%time_cpu,               &
!                                        obsr2_tim%time_cpu,               &
!                                        cldan_tim%time_cpu,               &
!                                        guess_tim%time_cpu,               &
!                                        gsub2_tim%time_cpu,               &
!                                        gsub3_tim%time_cpu,               &
!                                        gsub4_tim%time_cpu,               &
!                                       upsend_tim%time_cpu,             &
!                                       dnsend_tim%time_cpu,             &
!                                       weight_tim%time_cpu,             &
!                                       bfiltT_tim%time_cpu,             &
                                       input_tim%time_cpu,             &
                                      output_tim%time_cpu,             &
                                       rfilt_tim%time_cpu,              &
!                                       filt2an_tim%time_cpu,            &
!                                       aintp_tim%time_cpu,              &
!                                       intp_tim%time_cpu,               &
!                                       an2filt_tim%time_cpu,            &
!                                       output_tim%time_cpu,             &
!                                       graph_tim%time_cpu,               &
                                       total_tim%time_cpu
      endif
    end if

    bufsize = LEN(TRIM(buffer)) + 1
    buffer(bufsize:bufsize) = NEW_LINE(' ')

      if(multigrid_betafct) then
    write(header,"(A6,9(',',A10))") "mype",                            &
!                                     "init",                            &
                                     "binit",                           &
!!                                     "upsend",                          &
!!                                     "dnsend",                          &
                                     "input",                          &
                                     "output",                          &
!                                     "dnsen1",                          &
!                                     "dnsen2",                          &
!                                     "dnsen3",                          &
                                     "weight",                          &
                                     "bfiltT",                          &
                                     "bfilt",                           &
                                     "mgbft",                           &
!                                     "brns1",                           &
!                                     "brns2",                           &
!                                     "pack",                            &
!                                     "filt2an",                         &
!                                     "aintp",                           &
!                                     "intp",                            &
!                                     "an2filt",                         &
!                                     "output",                          &
                                     "bfint",                          &
!                                     "graph",                          &
                                     "total"
      else
    write(header,"(A6,4(',',A10))") "mype",                            &
!                                     "init",                            &
!                                     "part1",                            &
!                                     "part2",                            &
!                                     "gsub1",                            &
!                                     "obsr1",                            &
!                                     "obsr2",                            &
!                                     "cldan",                            &
!                                     "guess",                            &
!                                     "gsub2",                            &
!                                     "gsub3",                            &
!                                     "gsub4",                            &
!                                     "upsend",                          &
!                                     "dnsend",                          &
!                                     "weight",                          &
!                                     "bfiltT",                          &
                                     "input",                          &
                                     "outpyt",                          &
                                     "rfilt",                           &
!                                     "filt2an",                         &
!                                     "aintp",                           &
!                                     "intp",                            &
!                                     "an2filt",                         &
!                                     "output",                          &
!                                     "graph",                          &
                                     "total"
      endif

    header(bufsize:bufsize) = NEW_LINE(' ')
    disp = 0
    call MPI_File_write_at(fh, disp, header, bufsize, MPI_BYTE, stat, ierr)

    disp = (mype+1)*bufsize
    call MPI_File_write_at(fh, disp, buffer, bufsize, MPI_BYTE, stat, ierr)

    call MPI_File_close(fh, ierr)

  endsubroutine print_mg_timers
!-----------------------------------------------------------------------
  function wtime()
    use mpi
    real(r_kind) :: wtime
    wtime = MPI_Wtime()
  endfunction wtime
!-----------------------------------------------------------------------
  function ctime()
    real(r_kind) :: ctime
    call CPU_TIME(ctime)
  endfunction ctime
!-----------------------------------------------------------------------
                        endmodule mg_timers
