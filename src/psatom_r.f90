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
! Self-consistent pseudoatom calculation for fully-relativistic
! pseudopotentials. Orbital occupancies are always weighted with the
! scalar-relativistic averaging factors, even for partially-occupied levels.

 subroutine psatom_r(na,la,ea,fa,nv,it,rhoc,rho, &
&           rr,rcmax,mmax,iexc,etot,nproj,vpuns,lloc,vkb,evkb,ierr,okb, uua)

!na  principal quantum number array, dimension nv
!la  angular-momenta
!ea  eigenvalues (input starting guess, output)
!fa  occupancies
!rpk  radius of outermost peak of wave function
!nv  number of valence states
!it  number of iterations (output)
!rr  log radial mesh
!rcmax  maximum core radius for psp
!mmax  size of log grid
!iexc  exchange-correlation function to be used
!etot  pseudoatom total energy (output)
!nproj  number of VKB projectora to use for each l
!vpuns  unscreened semi-local pseudopotentials (plus differenv vloc if lloc==4)
!lloc  index-1 of local potential
!vkb   Vanderbilt-Kleinman-Bylander projectors
!evkb VKB projector coefficients
!okb 0,use semi-local, 1 use VKB projectors
!uua Array of pseudo-atomic orbitals (output)

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 
 integer :: mmax,iexc,nv,lloc,okb
 integer :: na(30),la(30),nproj(5)
 real(dp) :: rcmax
 real(dp) :: fa(30),rr(mmax)
 real(dp) :: vpuns(mmax,5),vkb(mmax,2,4,2),evkb(2,4,2)

!Output variables
 integer :: it
 real(dp) :: etot
 real(dp) :: ea(30,2)
 real(dp) :: rho(mmax),rhoc(mmax),vi(mmax)
 real(dp) :: uua(mmax,2,nv)

!Local variables
 integer :: nin,mch
 real(dp) :: amesh,al
 real(dp) :: dr,eeel,eexc,et,emin,emax,rl,rl1,sd,sn,sls,eeig,fj
 real(dp) :: thl,vn,zion
 integer :: ii,jj,l1,ll,ikap,kap,mkap,ierr
 logical :: convg

 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: vo(:),vi1(:),vo1(:),vxc(:),vp(:)
 real(dp), allocatable :: vtot(:)

! blend parameter for Anderson iterative potential mixing
 real(dp), parameter ::  bl=0.5d0

 allocate(uu(mmax),up(mmax))
 allocate(vo(mmax),vi1(mmax),vo1(mmax),vxc(mmax),vp(mmax))
 allocate(vtot(mmax))

! this seems necessary in sratom, so I might as well do it  here too
! don't ask why!
 uu(:)=0.d0; up(:)=0.d0
 vo(:)=0.d0; vi1(:)=0.d0; vo1(:)=0.d0; vxc(:)=0.d0; vp(:)=0.d0
 vtot(:)=0.d0
 dr=0.d0; eeel=0.d0; eexc=0.d0; et=0.d0; emin=0.d0; emax=0.d0;
 rl=0.d0; rl1=0.d0; sd=0.d0; sn=0.d0; sls=0.d0; eeig=0.d0; 
 thl=0.d0; vn=0.d0; zion=0.d0; 
 nin=0; mch=0 

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)

! total valence charge
 zion=0.0d0
 do ii=1,nv
   zion=zion+fa(ii)
 end do

! screening potential for pseudocharge
! input rho is assumed to be valence rho from all-electron calculation
! which shouldn''t be a bad approximation for a starting screened
! potential

 call vout(1,rho,rhoc,vi,vxc,zion,eeel,eexc,rr,mmax,iexc)

! big self  self-consietency loop

 do it=1,100

! add well do deal with initially unbound states which may arise from
! discrepancy between all-electron valence charge used to initiate screening
! and self-consistent pseudocharge
! start well at 0.5 Ha at infinity and scale down with iterations

   convg=.true.

! solve for bound states in turn

   eeig=0.0d0

   vtot(:)=vpuns(:,lloc+1)+vi(:)

   rho(:)=0.0d0

   do ii=1,nv
     ierr = 0
     l1=la(ii)+1
     ll = la(ii)
     if(ll==0) then
      mkap=1
     else
      mkap=2
     end if
     vp(:)=vpuns(:,l1)+vi(:)
! loop on J = ll +/- 1/2
     do ikap=1,mkap
      if(ikap==1) then
        kap=-(ll+1)
        fj=dble(ll+1)*fa(ii)/dble(2*ll+1)
      else
        kap=  ll
        fj=dble(ll)*fa(ii)/dble(2*ll+1)
      end if
      et=ea(ii,ikap)

! estimite bounds on eigenvalue based on semi-local pseudopotential
      sls=(l1-1)*l1
      emax=vtot(mmax)+0.5d0*sls/rr(mmax)**2
      emin=emax
      do jj=1,mmax
        emin=dmin1(emin,vpuns(jj,l1)+vi(jj)+0.5d0*sls/rr(jj)**2)
      end do

! option for using semi-local potentials
      if(okb==0) then
        call lschpb(na(ii),la(ii),ierr,et, &
&                   rr,vp,uu,up,mmax,mch)
      else if(okb==1) then
        call lschvkbb(na(ii),la(ii),nproj(l1),ierr,et,emin,emax, &
&                   rr,vtot,vkb(1,1,l1,ikap),evkb(1,l1,ikap),uu,up,mmax,mch)
      else
        write(6,'(a,i4)') 'psatom: argument ERROR - okb=',okb
        stop
      end if
      uua(:,ikap,ii) = uu

      if(ierr .ne. 0) then
       write(6,'(/a,4i4)') 'psatom: WARNING lschvkbb convergence error n,l,kap,iter=', &
&       na(ii),la(ii),kap,it
       exit
      end if

! overall convergence criterion based on eps within lschfb
      if(ea(ii,ikap)/= et) convg=.false.
      ea(ii,ikap)=et

! accumulate charge and eigenvalues
      eeig = eeig + fj*ea(ii,ikap)
      rho(:)=rho(:) + fj*(uu(:)/rr(:))**2

    end do !ikap
    if(ierr/=0) exit
   end do !ii

   if(ierr/=0) exit

! output potential
   call vout(1,rho,rhoc,vo,vxc,zion,eeel,eexc,rr,mmax,iexc)

! generate next iteration using d. g. anderson''s
! method
   thl=0.0d0
   if(it>1) then
     sn=0.0d0
     sd=0.0d0
     do ii=1,mmax
       rl=vo(ii)-vi(ii)
       rl1=vo1(ii)-vi1(ii)
       dr=rl-rl1
       sn=sn + rl*dr*rr(ii)**2
       sd=sd + dr*dr*rr(ii)**2
     end do
     thl=sn/sd
   end if

   do ii=1,mmax
     vn=(1.0d0-bl)*((1.0d0-thl)*vi(ii) + thl*vi1(ii)) &
  &   + bl*((1.0d0-thl)*vo(ii) + thl*vo1(ii))
     vi1(ii)=vi(ii)
     vo1(ii)=vo(ii)
     vi(ii)=vn
   end do

   if(convg) exit

   if(it==100 .and. .not. convg) then
     write(6,'(/a)') 'psatom: potential failed to converge'
     ierr=101
   end if

 end do !it

! total energy output

! output potential for e-e interactions

 call vout(1,rho,rhoc,vo,vxc,zion,eeel,eexc, &
&          rr,mmax,iexc)

 etot =  eeig + eexc - 0.5d0*eeel

 deallocate(uu,up)
 deallocate(vo,vi1,vo1,vxc,vp)
 deallocate(vtot)
 return

 end subroutine psatom_r
