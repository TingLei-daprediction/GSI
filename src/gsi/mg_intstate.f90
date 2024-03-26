!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_intstate
!***********************************************************************
!                                                                      !
! Contains declarations and allocations of internal state variables    !
! use for filtering                                                    !
!                       - offset version -                             !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use jp_pkind2, only: fpi
use mpimod, only: mype
!STA use mg_mppstuff, only: mype
use mg_parameter, only: im,jm,nh,hx,hy,pasp01,pasp02,pasp03
use mg_parameter, only: lm,hz,p,km,km2,km3,km,nm,mm,ib,jb,nb,mb
use mg_parameter, only: lmf,lmh,km,km
use berror, only: mg_weig1,mg_weig2,mg_weig3,mg_weig4
!STA use mg_parameter, only: mg_weig1,mg_weig2,mg_weig3,mg_weig4
use mg_mppstuff, only: my_hgen,finishMPI,barrierMPI
use jp_pbfil,only: cholaspect
use jp_pbfil,only: getlinesum
use jp_pbfil3, only: inimomtab,t22_to_3,tritform,t33_to_6,hextform
!TEST
!use gridmod, only: lat1,lon1
!TEST
implicit none

real(r_kind), allocatable,dimension(:,:,:):: V
!
! Composite control variable on first generation o filter grid
!
real(r_kind), allocatable,dimension(:,:,:):: VALL
real(r_kind), allocatable,dimension(:,:,:):: HALL
!
! Composite control variable on high generations of filter grid
!
!
!FOR ADJOINT TEST
!
!real(r_kind), allocatable,dimension(:,:):: A
!real(r_kind), allocatable,dimension(:,:):: B
!real(r_kind), allocatable,dimension(:,:):: A0
!real(r_kind), allocatable,dimension(:,:):: B0
!
real(r_kind), allocatable,dimension(:,:,:):: a_diff_f
real(r_kind), allocatable,dimension(:,:,:):: a_diff_h
real(r_kind), allocatable,dimension(:,:,:):: b_diff_f
real(r_kind), allocatable,dimension(:,:,:):: b_diff_h

real(r_kind), allocatable,dimension(:,:):: p_eps
real(r_kind), allocatable,dimension(:,:):: p_del
real(r_kind), allocatable,dimension(:,:):: p_sig
real(r_kind), allocatable,dimension(:,:):: p_rho

real(r_kind), allocatable,dimension(:,:,:):: paspx
real(r_kind), allocatable,dimension(:,:,:):: paspy
real(r_kind), allocatable,dimension(:,:,:):: pasp1
real(r_kind), allocatable,dimension(:,:,:,:):: pasp2
real(r_kind), allocatable,dimension(:,:,:,:,:):: pasp3

real(r_kind), allocatable,dimension(:,:,:):: vpasp2
real(r_kind), allocatable,dimension(:,:,:):: hss2
real(r_kind), allocatable,dimension(:,:,:,:):: vpasp3
real(r_kind), allocatable,dimension(:,:,:,:):: hss3

real(r_kind), allocatable,dimension(:):: ssx
real(r_kind), allocatable,dimension(:):: ssy
real(r_kind), allocatable,dimension(:):: ss1
real(r_kind), allocatable,dimension(:,:):: ss2
real(r_kind), allocatable,dimension(:,:,:):: ss3

integer(fpi), allocatable,dimension(:,:,:):: dixs
integer(fpi), allocatable,dimension(:,:,:):: diys
integer(fpi), allocatable,dimension(:,:,:):: dizs

integer(fpi), allocatable,dimension(:,:,:,:):: dixs3
integer(fpi), allocatable,dimension(:,:,:,:):: diys3
integer(fpi), allocatable,dimension(:,:,:,:):: dizs3

integer(fpi), allocatable,dimension(:,:,:,:):: qcols

!
!
! Composite stacked variable
!

!STA real(r_kind), allocatable,dimension(:,:,:):: WORKA


integer(i_kind),allocatable,dimension(:):: iref,jref
integer(i_kind),allocatable,dimension(:):: Lref,Lref_h
real(r_kind),allocatable,dimension(:):: cvf1,cvf2,cvf3,cvf4
real(r_kind),allocatable,dimension(:):: cvh1,cvh2,cvh3,cvh4

real(r_kind),allocatable,dimension(:):: cx0,cx1,cx2,cx3
real(r_kind),allocatable,dimension(:):: cy0,cy1,cy2,cy3

real(r_kind),allocatable,dimension(:):: p_coef,q_coef
real(r_kind),allocatable,dimension(:):: a_coef,b_coef

real(r_kind),allocatable,dimension(:,:):: cf00,cf01,cf02,cf03           &
                                         ,cf10,cf11,cf12,cf13           &
                                         ,cf20,cf21,cf22,cf23           &
                                         ,cf30,cf31,cf32,cf33

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine allocate_mg_intstate
!***********************************************************************
!                                                                      !
! Allocate internal state variables                                    !
!                                                                      !
!***********************************************************************

allocate(V(1-hx:im+hx,1-hy:jm+hy,lm))                     ; V=0.
allocate(VALL(km,1-hx:im+hx,1-hy:jm+hy))                  ; VALL=0.
allocate(HALL(km,1-hx:im+hx,1-hy:jm+hy))                  ; HALL=0.


allocate(a_diff_f(km,1-hx:im+hx,1-hy:jm+hy))                ; a_diff_f=0. 
allocate(a_diff_h(km,1-hx:im+hx,1-hy:jm+hy))                ; a_diff_h=0. 
allocate(b_diff_f(km,1-hx:im+hx,1-hy:jm+hy))                ; b_diff_f=0. 
allocate(b_diff_h(km,1-hx:im+hx,1-hy:jm+hy))                ; b_diff_h=0. 

allocate(p_eps(1-hx:im+hx,1-hy:jm+hy))                            ; p_eps=0.
allocate(p_del(1-hx:im+hx,1-hy:jm+hy))                            ; p_del=0.
allocate(p_sig(1-hx:im+hx,1-hy:jm+hy))                            ; p_sig=0.
allocate(p_rho(1-hx:im+hx,1-hy:jm+hy))                            ; p_rho=0.

allocate(paspx(1,1,1:im))                                          ; paspx=0.
allocate(paspy(1,1,1:jm))                                          ; paspy=0.

allocate(pasp1(1,1,1:lm))                                           ; pasp1=0.
allocate(pasp2(2,2,1:im,1:jm))                                    ; pasp2=0.
allocate(pasp3(3,3,1:im,1:jm,1:lm))                               ; pasp3=0.

allocate(vpasp2(0:2,1:im,1:jm))                                   ; vpasp2=0.
allocate(hss2(1:im,1:jm,1:3))                                     ; hss2= 0.

allocate(vpasp3(1:6,1:im,1:jm,1:lm))                              ; vpasp3= 0.
allocate(hss3(1:im,1:jm,1:lm,1:6))                                ; hss3= 0.

allocate(ssx(1:im))                                             ; ssx=0.
allocate(ssy(1:jm))                                             ; ssy=0.
allocate(ss1(1:lm))                                             ; ss1=0.
allocate(ss2(1:im,1:jm))                                        ; ss2=0.
allocate(ss3(1:im,1:jm,1:lm))                                   ; ss3=0.

allocate(dixs(1:im,1:jm,3))                                     ; dixs=0
allocate(diys(1:im,1:jm,3))                                     ; diys=0

allocate(dixs3(1:im,1:jm,1:lm,6))                               ; dixs3=0
allocate(diys3(1:im,1:jm,1:lm,6))                               ; diys3=0
allocate(dizs3(1:im,1:jm,1:lm,6))                               ; dizs3=0

allocate(qcols(0:7,1:im,1:jm,1:lm))                             ; qcols=0

!
! In standalone version
!

!
!STA allocate(WORKA(km,1:nm,1:mm))                               ; WORKA=0.

!
! for re-decomposition
!

allocate(iref(1:nm))                                     ; iref=0
allocate(jref(1:mm))                                     ; jref=0

allocate(cx0(1:nm))                                      ; cx0=0.
allocate(cx1(1:nm))                                      ; cx1=0.
allocate(cx2(1:nm))                                      ; cx2=0.
allocate(cx3(1:nm))                                      ; cx3=0.

allocate(cy0(1:mm))                                      ; cy0=0.
allocate(cy1(1:mm))                                      ; cy1=0.
allocate(cy2(1:mm))                                      ; cy2=0.
allocate(cy3(1:mm))                                      ; cy3=0.

allocate(p_coef(4))                                      ; p_coef=0.
allocate(q_coef(4))                                      ; q_coef=0.

allocate(a_coef(3))                                      ; a_coef=0.
allocate(b_coef(3))                                      ; b_coef=0.


allocate(cf00(1:nm,1:mm))                            ; cf00=0.
allocate(cf01(1:nm,1:mm))                            ; cf01=0.
allocate(cf02(1:nm,1:mm))                            ; cf02=0.
allocate(cf03(1:nm,1:mm))                            ; cf03=0.
allocate(cf10(1:nm,1:mm))                            ; cf10=0.
allocate(cf11(1:nm,1:mm))                            ; cf11=0.
allocate(cf12(1:nm,1:mm))                            ; cf12=0.
allocate(cf13(1:nm,1:mm))                            ; cf13=0.
allocate(cf20(1:nm,1:mm))                            ; cf20=0.
allocate(cf21(1:nm,1:mm))                            ; cf21=0.
allocate(cf22(1:nm,1:mm))                            ; cf22=0.
allocate(cf23(1:nm,1:mm))                            ; cf23=0.
allocate(cf30(1:nm,1:mm))                            ; cf30=0.
allocate(cf31(1:nm,1:mm))                            ; cf31=0.
allocate(cf32(1:nm,1:mm))                            ; cf32=0.
allocate(cf33(1:nm,1:mm))                            ; cf33=0.

allocate(Lref(1:lm))                                 ; Lref=0
allocate(Lref_h(1:lmf))                              ; Lref_h=0

allocate(cvf1(1:lm))                                 ; cvf1=0.
allocate(cvf2(1:lm))                                 ; cvf2=0.
allocate(cvf3(1:lm))                                 ; cvf3=0.
allocate(cvf4(1:lm))                                 ; cvf4=0.

allocate(cvh1(1:lmf))                                ; cvh1=0.
allocate(cvh2(1:lmf))                                ; cvh2=0.
allocate(cvh3(1:lmf))                                ; cvh3=0.
allocate(cvh4(1:lmf))                                ; cvh4=0.


!-----------------------------------------------------------------------
                        endsubroutine allocate_mg_intstate

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine def_mg_weights
!***********************************************************************
!                                                                      !
! Define weights and scales                                            !
!                                                                      !
!***********************************************************************
integer(i_kind):: i,j,L,k,n
!-----------------------------------------------------------------------

      p_eps(:,:)=0.0
      p_del(:,:)=0.0
      p_sig(:,:)=0.0
      p_rho(:,:)=0.0

!--------------------------------------------------------
    do n=1,km3
      do k=(n-1)*lm+1,n*lm
        a_diff_f(k,:,:)=mg_weig1(n)
        b_diff_f(k,:,:)=0.
      enddo
   enddo
   do n=km3+1,km3+km2
      k=km3*(lm-1)+n
        a_diff_f(k,:,:)=mg_weig1(n)
        b_diff_f(k,:,:)=0.
   enddo
   

      select case(my_hgen)
        case(2) 
          do n=1,km3
            do k=(n-1)*lm+1,n*lm
              a_diff_h(k,:,:)=mg_weig2(n)
              b_diff_h(k,:,:)=0.
            enddo
          enddo
          do n=km3+1,km3+km2
             k=km3*(lm-1)+n
               a_diff_h(k,:,:)=mg_weig2(n)
               b_diff_h(k,:,:)=0.
          enddo
        case(3) 
          do n=1,km3
            do k=(n-1)*lm+1,n*lm
              a_diff_h(k,:,:)=mg_weig3(n)
              b_diff_h(k,:,:)=0.
            enddo
          enddo
          do n=km3+1,km3+km2
             k=km3*(lm-1)+n
               a_diff_h(k,:,:)=mg_weig3(n)
               b_diff_h(k,:,:)=0.
          enddo
        case default 
          do n=1,km3
            do k=(n-1)*lm+1,n*lm
              a_diff_h(k,:,:)=mg_weig4(n)
              b_diff_h(k,:,:)=0.
            enddo
          enddo
          do n=km3+1,km3+km2
             k=km3*(lm-1)+n
               a_diff_h(k,:,:)=mg_weig4(n)
               b_diff_h(k,:,:)=0.
          enddo
      end select

          do L=1,lm
           pasp1(1,1,L)=pasp01
          enddo

          do i=1,im
            paspx(1,1,i)=pasp02
          enddo  
          do j=1,jm
            paspy(1,1,j)=pasp02
          enddo  

          do j=1,jm
          do i=1,im
            pasp2(1,1,i,j)=pasp02*(1.+p_del(i,j))
            pasp2(2,2,i,j)=pasp02*(1.-p_del(i,j))
            pasp2(1,2,i,j)=pasp02*p_eps(i,j)     
            pasp2(2,1,i,j)=pasp02*p_eps(i,j)     
          end do
          end do

        do L=1,lm
          do j=1,jm
          do i=1,im
            pasp3(1,1,i,j,l)=pasp03*(1+p_del(i,j))
            pasp3(2,2,i,j,l)=pasp03
            pasp3(3,3,i,j,l)=pasp03*(1-p_del(i,j))
            pasp3(1,2,i,j,l)=pasp03*p_eps(i,j)
            pasp3(2,1,i,j,l)=pasp03*p_eps(i,j)
            pasp3(2,3,i,j,l)=pasp03*p_sig(i,j)
            pasp3(3,2,i,j,l)=pasp03*p_sig(i,j)
            pasp3(1,3,i,j,l)=pasp03*p_rho(i,j)
            pasp3(3,1,i,j,l)=pasp03*p_rho(i,j)
          end do
          end do
        end do


        call cholaspect(1,lm,pasp1)
        call cholaspect(1,im,1,jm,pasp2)
        call cholaspect(1,im,1,jm,1,lm,pasp3)


        call getlinesum(hx,1,im,paspx,ssx)
        call getlinesum(hy,1,jm,paspy,ssy)
        call getlinesum(hz,1,lm,pasp1,ss1)
        call getlinesum(hx,1,im,hy,1,jm,pasp2,ss2)
        call getlinesum(hx,1,im,hy,1,jm,hz,1,lm,pasp3,ss3)
!-----------------------------------------------------------------------
                        endsubroutine def_mg_weights

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine init_mg_line
!***********************************************************************
!                                                                      !
! Inititate line filters                                               !
!                                                                      !
!***********************************************************************
integer(i_kind):: i,j,L,icol
logical:: ff
!-----------------------------------------------------------------------

  do j=1,jm
  do i=1,im
    call t22_to_3(pasp2(:,:,i,j),vpasp2(:,i,j))
  enddo
  enddo

  do l=1,lm
  do j=1,jm
  do i=1,im
    call t33_to_6(pasp3(:,:,i,j,l),vpasp3(:,i,j,l))
  enddo
  enddo
  enddo



  call inimomtab(p,nh,ff)

  call tritform(1,im,1,jm,vpasp2, dixs,diys, ff)

  do icol=1,3
    hss2(:,:,icol)=vpasp2(icol-1,:,:)
  enddo  


  call hextform(1,im,1,jm,1,lm,vpasp3,qcols,dixs3,diys3,dizs3, ff)


  do icol=1,6
    hss3(:,:,:,icol)=vpasp3(icol,:,:,:)
  enddo
 

!-----------------------------------------------------------------------
                        endsubroutine init_mg_line

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine deallocate_mg_intstate
!***********************************************************************
!                                                                      !
! Deallocate internal state variables                                  !
!                                                                      !
!***********************************************************************

deallocate(V)

deallocate(HALL,VALL)

deallocate(a_diff_f,b_diff_f)
deallocate(a_diff_h,b_diff_h)
deallocate(p_eps,p_del,p_sig,p_rho,pasp1,pasp2,pasp3,ss1,ss2,ss3)
deallocate(dixs,diys)
deallocate(dixs3,diys3,dizs3)
deallocate(qcols)
!
! for testing
!
!STAT deallocate(WORKA)

!
! for re-decomposition
!
deallocate(iref,jref)

deallocate(cf00,cf01,cf02,cf03,cf10,cf11,cf12,cf13)
deallocate(cf20,cf21,cf22,cf23,cf30,cf31,cf32,cf33)

deallocate(Lref,Lref_h)

deallocate(cvf1,cvf2,cvf3,cvf4)

deallocate(cvh1,cvh2,cvh3,cvh4)

deallocate(cx0,cx1,cx2,cx3)
deallocate(cy0,cy1,cy2,cy3)

deallocate(p_coef,q_coef)
deallocate(a_coef,b_coef)



                        endsubroutine deallocate_mg_intstate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_intstate
