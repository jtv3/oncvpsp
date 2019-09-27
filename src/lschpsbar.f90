!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine lschpsbar(nn,ll,ierr,ee,emin,emax,rr,vv,uu,up,mmax,mbar)

! Finds bound states of a semi-local pseudopotential

!nn  principal quantum number
!ll  angular-momentum quantum number
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value (in/our)
!emin  lower bound, potential minimum if ==0.0
!emax  upper bound
!rr  log radial mesh
!vv  local psp
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mbar  mesh point for infinite barrier

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr(mmax),vv(mmax)
 real(dp) :: emin,emax
 integer :: nn,ll,mmax,mbar

!Output variables
 real(dp) :: uu(mmax),up(mmax)
 real(dp) :: ee
 integer :: ierr

!Local Variables

 real(dp) :: aei,aeo,aii,aio !functions in aeo.f90
 real(dp) :: dlte,de
 real(dp) :: eps,exp,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap,xkap2
 real(dp) :: amesh,al,als
 integer :: ii, it,mch,mchb,nint,node,nin
 logical :: maxset

 real(dp), allocatable :: upp(:),cf(:)
 allocate(upp(mmax),cf(mmax))

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
 eps=1.0d-8
 ierr = 60

 sls=ll*(ll+1)

 dlte=emax-emin
 if(emin==0.0d0) then
   do ii=1,mbar
     emin=dmin1(emin,vv(ii)+0.5d0*sls/rr(ii)**2)
   end do
   emax=emin+dlte
 end if

 ee=0.5d0*(emin+emax)
 node=0
 de=0.0d0

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0
 maxset=.false.

 als=al**2

! return point for bound state convergence
 do nint=1,60
  
! coefficient array for u in differential eq.
   do ii=1,mmax
     cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
   end do
  
! find classical turning point for matching
   if(cf(mbar-1)<0.0d0) then
     mch=mbar-1
   else
     do ii=mbar-1,2,-1
       if(cf(ii-1)<=0.d0 .and. cf(ii)>0.d0) then
         mch=ii
         exit
       end if
     end do
   end if

   if(mch==0) then
    write(6,'(/a)') 'lschpbar: ERROR no classical turning point'
    stop
   end if
  
   de=0.0d0
! start wavefunction with series
  
   do ii=1,4
     uu(ii)=rr(ii)**(ll+1)
     up(ii)=al*(ll+1)*rr(ii)**(ll+1)
     upp(ii)=al*up(ii)+cf(ii)*uu(ii)
   end do
  
! outward integration using predictor once, corrector
! twice
   node=0
  
   do ii=4,mch-1
     uu(ii+1)=uu(ii)+aeo(up,ii)
     up(ii+1)=up(ii)+aeo(upp,ii)
     do it=1,2
       upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1)
       up(ii+1)=up(ii)+aio(upp,ii)
       uu(ii+1)=uu(ii)+aio(up,ii)
     end do
     if(uu(ii+1)*uu(ii) .le. 0.0d0) node=node+1
   end do
  
   mchb=min(mch,mbar-1)
   uout=uu(mchb)
   upout=up(mchb)
  
  
   if(node-nn+ll+1==0) then
    if(maxset) then

  
! start inward integration at barrier
! point with sinh or sin
  
     xkap2=sls/rr(mbar)**2 + 2.0d0*(vv(mbar)-ee)
  
     do ii=mbar,mbar+4
       if(xkap2>0) then
         xkap=sqrt(xkap2)
         uu(ii)=sinh(-xkap*(rr(ii)-rr(mbar)))
         up(ii)=-rr(ii)*al*xkap*cosh(-xkap*(rr(ii)-rr(mbar)))
       else
         xkap=sqrt(-xkap2)
         uu(ii)=sin(-xkap*(rr(ii)-rr(mbar)))
         up(ii)=-rr(ii)*al*xkap*cos(-xkap*(rr(ii)-rr(mbar)))
       end if
       upp(ii)=al*up(ii)+cf(ii)*uu(ii)
     end do
  
! integrate inward
  
     do ii=mbar,mchb+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2
         upp(ii-1)=al*up(ii-1)+cf(ii-1)*uu(ii-1)
         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do

     do ii=mbar+1,mbar+4
       uu(ii)=0.0d0
       up(ii)=0.0d0
       upp(ii)=0.0d0
     end do
  
! scale outside wf for continuity
  
     sc=uout/uu(mchb)
  
     do ii=mchb,mbar
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do
  
     upin=up(mchb)
  
! perform normalization sum
     nin=mbar

     ro=rr(1)/dsqrt(amesh)
     sn=ro**(2*ll+3)/dfloat(2*ll+3)
  
     do ii=1,nin-3
       sn=sn+al*rr(ii)*uu(ii)**2
     end do
  
     sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2)**2 &
&              + 28.0d0*rr(nin-1)*uu(nin-1)**2 &
&              +  9.0d0*rr(nin  )*uu(nin  )**2)/24.0d0
  
! normalize u
  
     cn=1.0d0/dsqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin
  
     do ii=1,nin
       up(ii)=cn*up(ii)
       uu(ii)=cn*uu(ii)
     end do
     do ii=nin+1,mmax
       uu(ii)=0.0d0
     end do
  
! perturbation theory for energy shift
  
     de=0.5d0*uout*(upout-upin)/(al*rr(mchb))

! convergence test and possible exit
  
     if(dabs(de)<eps) then
       ierr = 0
       exit
     end if
  
     if(de>0.0d0) then
       emin=ee
     else
       emax=ee
     end if
     ee=ee+de
     if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)
    else !.not. maxset
     emax=emax+dlte
     ee=0.5d0*(emin+emax)
    end if !maxset
   else if(node-nn+ll+1<0) then

! too few nodes
     emin=ee
     if(.not. maxset) then
      emax=emax+dlte
      ee=emax
     end if
     ee=0.5d0*(emin+emax)
   else

! too many nodes
     maxset=.true.
     emax=ee
     ee=0.5d0*(emin+emax)
   end if
  
 end do !nint

 deallocate(upp,cf)
 return

 end subroutine lschpsbar
