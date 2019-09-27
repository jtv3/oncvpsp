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
! self-consistent all-electron and pseudopotential atomic calculations
! compared for reference and tests atomic configurations

 subroutine run_config(na,la,ea,fa,nc,nv,rhov,rhomod,rr,zz,rcmax,mmax, &
&                      iexc,etot,epstot,nproj,vpuns,lloc,vkb,evkb,srel, uua)

!na  principal quantum number array, dimension nv
!la  angular-momenta
!ea  eigenvalues (input starting guess, output)
!fa  occupancies
!rpk  radius of outermost peak of wave function
!nc  number of core states
!nv  number of valence states
!rhov valence pseudo-charge divided by reference configuration zval
!rhomod  model core charge
!rr  log radial mesh
!zz  atomic number
!rcmax  maximum core radius for psp
!mmax  size of log grid
!iexc  exchange-correlation function to be used
!etot  all-electron atom otal energy
!epstot  pseudoatom total energy
!nproj  number of VKB projectora to use for each l
!vpuns  unscreened semi-local pseudopotentials (plus different vloc if lloc==4)
!lloc  index-1 of local potential
!vkb   Vanderbilt-Kleinman-Bylander projectors
!evkb VKB projector coefficients
!srel .true. for scalar-relativistic, .false. for non-relativistic
!uua Array of pseudo-atomic orbitals (output)

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 
 integer :: mmax,iexc,nc,nv,lloc
 integer :: na(30),la(30),nproj(5)
 real(dp) :: etot,epstot,rcmax,zz
 real(dp) :: fa(30),rhov(mmax),rr(mmax)
 real(dp) :: vpuns(mmax,5),vkb(mmax,2,4),evkb(2,4),rhomod(mmax,5)
 logical :: srel

!Output variables  only printing
 real(dp) :: uua(mmax,nv)
 
!Local variables
 integer :: ii,it,l1,ierr,okb
 integer :: ntstp(10),ltstp(10),ninit(4)
 real(dp) :: eaetst,etsttot,zval
 real(dp) :: ea(30),rpk(30)
 real(dp) :: ftstp(10),etstp(10)
 real(dp),allocatable :: rho(:),rhoc(:),rhocps(:),vi(:),vfull(:)

 allocate(rho(mmax),rhoc(mmax),vi(mmax),vfull(mmax),rhocps(mmax))

! atom tests for various configurations includins the reference configuration

! full potential atom solution
 ea(:)=0.d0; rpk(:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0 
 eaetst=0.d0

 call sratom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
&            rr,vfull,zz,mmax,iexc,eaetst,ierr,srel)

 if(ierr/=0) then
  write(6,'(a/a)') 'run_config: WARNING  for AE atom,', &
&       ' no output for this configuration'
  if(ierr==-1) then
    write(6,'(a,i4)') 'run_config: WARNING no classical turning point error, iteration',it
  else
    write(6,'(a)') 'run_config WARNING self-consistency failed to converge'
  end if
  deallocate(rho,rhoc,rhocps,vi,vfull)
  return
 end if

 do l1=1,4
  ninit(l1)=l1
 end do

 zval=0.0d0
 do ii=1,nv
  ltstp(ii)=la(nc+ii)
  l1=la(nc+ii)+1
  ntstp(ii)=ninit(l1)
  ninit(l1)=ninit(l1)+1
  ftstp(ii)=fa(nc+ii)
  zval=zval+fa(nc+ii)
  etstp(ii)=ea(nc+ii)
 end do

 rhocps(:)=rhomod(:,1)

! rescaled valence charge density of reference configuration used as starting
! approximation for pseudopotential atom calculation
 rho(:)=zval*rhov(:)

 okb=1
 call psatom(ntstp,ltstp,etstp,ftstp,nv,it,rhocps,rho, &
&           rr,rcmax,mmax,iexc,etsttot,nproj,vpuns,lloc, &
&           vkb,evkb,ierr,okb, uua)

 if(ierr/=0) then
  write(6,'(a/a)') 'run_config: ERROR for fully non-local  PS atom,', &
&       ' no output for this configuration'
  deallocate(rho,rhoc,rhocps,vi,vfull)
  return
 end if

 write(6,'(/a)') '   n   l     f        eae           eps        diff'
 do ii=1,nc
  write(6,'(2i4,f8.4,f14.8)') na(ii),la(ii),fa(ii),ea(ii)
 end do
 do ii=1,nv
  write(6,'(2i4,f8.4,2f14.8,1p,d12.2)') na(ii+nc),la(ii+nc),fa(ii+nc), &
&  ea(ii+nc),etstp(ii),etstp(ii)-ea(ii+nc)
 end do

 write(6,'(/a)') '    Total energies and differences'
 write(6,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      AE_ref=',etot, &
& '  AE_tst=',eaetst,'  dif=',eaetst-etot
 write(6,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      PS_ref=',epstot, &
& '  PS_tst=',etsttot,'  dif=',etsttot-epstot
 write(6,'(a,1p,d10.2)'), '      PSP excitation error=', &
& eaetst-etot-etsttot+epstot

 deallocate(rho,rhoc,rhocps,vi,vfull)
 return
 end subroutine run_config
