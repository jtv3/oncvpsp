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
 subroutine wellstate_r(nn,ll,kap,irc,ep,rr,vfull,uu,up,zz,mmax,mch)

!creates quantum well to confine positive-energy state, and calculates
!the resulting all-electron wave function
!Fully relativistic case using Dirac equation

!nn  principal quantum number of well state
!ll  angular momentum
!kap =l, -(l+1) for j=l -/+ 1/2
!irc  index of core radius
!ep  target energy for well state (>0)
!rr  log radial mesh
!vfull  all-electron potential
!uu  all-electron well-bound wave function
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log radial mesh 
!mch matching mesh point for inward-outward integrations
 
 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr(mmax),vfull(mmax)
 real(dp) :: ep,zz
 integer :: nn,ll,kap,irc,mmax

!Output variables
 real(dp) :: uu(mmax,2),up(mmax,2)
 integer :: mch
 
!Local variables
 real(dp), allocatable :: vwell(:)
 real(dp) :: al,cwell,et,xx,rwell
 real(dp) :: cwmax,cwmin

 real(dp),parameter :: eps=1.0d-8
 integer :: ii,iter,ierr,l1,itrwell
 logical :: convg
 
 l1=ll+1
 al = 0.01d0 * dlog(rr(101) / rr(1))

 allocate(vwell(mmax))

 rwell=8.0d0*rr(irc)

 convg=.false.

 do itrwell=1,20

  cwmax=ep+1.0d0
  cwmin=ep
  
  do iter=1,100

!create well potential
   do ii=1,mmax
    vwell(ii)=vfull(ii)
   end do
   cwell=0.5d0*(cwmax+cwmin)

! start outside rc to keep numerical derivatives at rc accurate
   do ii=irc+6,mmax
    xx=(rr(ii)-rr(irc+5))/rwell
    vwell(ii)=vwell(ii)+cwell*xx**3/(1.0d0+xx**3)
   end do

!find bound state in well
   et=ep
   call ldiracfb(nn,ll,kap,ierr,et,rr,zz,vwell,uu,up,mmax,mch)

   if(ierr .ne. 0) then
    write(6,'(a,i3,a,i3)') 'wellstate: lschfb convergence ERROR, n=',nn,' l=',ll
    stop
   end if

   if(abs(et-ep)<eps) then
    ep=et
    convg=.true.
    exit
   end if

   if(et<0.0d0 .and. cwell==0.0d0) then
    ep=et
    write(6,'(/a,i3,a/a,f12.6,a)') 'wellstate: Negative-energy state for l=', &
&         ll,' found for maximum barrier.',  'Normal bound state at ep=',ep, &
&        ' will be used.'
    convg=.true.
    exit
   end if

!if a bound state is found for the highest repulsive barrier to be used,
!set the barrier to zero and use the actual bound state
!otherwise, interval-halve

   if(et<0.0d0 .and. iter==1 .and. itrwell==20) then
    cwmax=0.0d0
    cwmin=0.0d0
   else
    if(et>ep) then
     cwmax=cwell
    else
     cwmin=cwell
    end if
   end if

  end do
  if(convg) then
   exit
  end if

!scale rwell if well potential wasn't converged
  rwell=0.75d0*rwell
 end do

 if(.not. convg) then
  write(6,'(/a)') 'wellstate: ERROR iteration failed to  converge'
   write(6,'(a,a,i3,a,i3)') 'wellstate: ERROR well potential iteration failed', &
&   ' to converge, n=',nn,' l=',ll
  stop
 end if

 if(cwell>0.0d0) then
  write(6,'(/a,i3,a,i3,a,i3/a,f8.4,/a,f8.4/a,f8.4)') &
&          '  Wellstate for l =',ll,'  n =',nn,'  kappa=',kap,  &
&          '    eigenvalue = ',et, &
&          '    asymptotic potential = ',cwell, &
&          '    half-point radius =',rr(irc)+rwell
 end if

 deallocate(vwell)
 return
 end subroutine wellstate_r
