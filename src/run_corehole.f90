! John Vinson, NIST
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
 subroutine run_corehole( rr,zz,mmax,srel,nc,nv,na,la,ea,fa,iexc, keep_valence )

! computes the core-hole potential, screened by the relaxation of the other core electrons
! TODO: mimic Shirley's atomic code options for 

!rr  log radial grid
!zz  atomic number
!mmax  size of radial grid
!srel .true. for scalar-relativistic, .false. for non-relativistic
!nc number of core levels
!nv number of valence levels
!na principle quantum numbers of core levels
!la angular quantum numbers of core levels
!ea energies
!fa occupations
!iexc exchange-correlation potential
!keep_valence flag sets if valence density is included in the calculation of the core relaxation

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: zz
 real(dp),intent(in) :: rr(mmax)
 integer,intent(in) :: mmax
 logical,intent(in) :: srel
 integer,intent(in) :: nc
 integer,intent(in) :: nv
 integer,intent(in) :: na(nc+nv),la(nc+nv)
 real(dp),intent(in) :: ea(nc+nv)
 real(dp),intent(in) :: fa(nc+nv)
 integer, intent(in) :: iexc
 logical,intent(in) :: keep_valence


 integer :: ii, ierr, it, jj
 real(dp) :: tmp_fa(30), rpk(30)
 real(dp), allocatable :: vfull(:), rhoc(:), rho(:), vo(:), vo2(:), rhoval(:)
 real(dp) :: etot, etot2, sume, sf_in

 character(len=17) :: filename

 allocate( vfull( mmax ), rhoc( mmax ), rho(mmax ), vo(mmax), vo2(mmax), rhoval(mmax) )

 tmp_fa(:) = 0.0_DP

 if( keep_valence) then
   call sratom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
  &            rr,vfull,zz,mmax,iexc,etot,ierr,srel)
   rhoval(:) = rho(:) - rhoc(:)
!#if DEBUG
!   open(unit=99,file='derp.txt')
!   do ii = 1, mmax
!     write(99,'(4E24.12)') rr(ii), rho(ii), rhoc(ii), rhoval(ii)
!   enddo
!   close(99)
!#endif
   if( nv .eq. 0 ) then
    sf_in = 0.0_DP
   else
     sf_in = sum( fa(nc+1:nc+nv) )
   endif
 else
   rhoval(:) = 0.0_DP
   sf_in = 0.0_DP
 endif

 sume = sum(fa(:))
 do ii = 1, nc
   write(6,*) 'OCEAN ---', na(ii), la(ii), fa(ii)  
   tmp_fa(1:nc) = fa(1:nc)
   call sratom_den(na,la,ea,tmp_fa,rpk,nc,nc+nv,it,rhoc,rho, rhoval, sf_in, &
&              rr,vfull,zz,mmax,iexc,etot,ierr,srel)
   call voutHartree(rho,vo,rr,mmax,sume)

   tmp_fa(ii) = tmp_fa(ii)-1.0_DP
   call sratom_den(na,la,ea,tmp_fa,rpk,nc,nc+nv,it,rhoc,rho, rhoval, sf_in, &
&              rr,vfull,zz,mmax,iexc,etot2,ierr,srel)

   call voutHartree(rho,vo2,rr,mmax,sume-1.0_DP)
   
   write(6,*) 'OCEAN ---', etot, etot2, it

!#if DEBUG
!   write(filename, '(A,I2.2)' ) 'test', ii
!   open(file=filename, unit=99, form='formatted' )
!   do jj = 1, mmax
!    write(99,'(4E24.12)') rr(jj), vo(jj), vo2(jj), rho(jj)
!   enddo
!   close(99)
!#endif

   write(filename, '(A8,I3.3,A1,I2.2,A1,I2.2)' ) 'vc_barez', int(zz), 'n', na(ii), 'l', la(ii )
    open(file=filename, unit=99, form='formatted' )
   do jj = 1, mmax
    write(99,'(4E24.15)') rr(jj), vo2(jj) - vo(jj)
   enddo
   close(99)
    

 enddo


 deallocate( vfull, rhoc, rho, vo, vo2, rhoval )
end subroutine run_corehole
