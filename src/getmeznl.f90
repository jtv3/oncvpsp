! John Vinson, NIST
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
 subroutine getmeznl( zz, nn, ll, irc, nopf, maxopf, rr, coreuu, aepr, mels )
! Writes out the matrix elements between the projectors and the core orbitals for 
! various powers of r. Currently 0-3 hardwired

!zz  atomic number
!nn principle quantum number
!ll angular momentum l
!irc  size of the projectors
!nopf number of optimal projectors
!maxopf max number of possible projectors (likely around 16 )
!rr  log radial grid
!coreuu core level orbital
!aepr all-electron opfs

 implicit none
 integer, parameter :: dp=kind(1.0d0)


!Input variables
 integer, intent(in) :: nn,ll,irc,nopf,maxopf
 real(dp),intent(in) :: zz
 real(dp),intent(in) :: rr(irc),coreuu(irc)
 real(dp),intent(in) :: aepr(irc,nopf)

!Output variables 
 real(dp) :: mels(maxopf,0:3)

!Local variables

 real(dp) :: al,ff(0:4)
 integer :: ip, io, i1, ii, jj

 mels(:,:) = 0.0_dp
 al = 0.01d0 * dlog(rr(101)/rr(1))

 do ip = 0, 3
   do io = 1, nopf
     do i1 = 1, irc - 4, 4
       do jj = 0, 4
         ii = i1 + jj
         ff( jj ) = coreuu( ii ) * aepr( ii, io ) * rr( ii ) ** ip
       enddo

       mels( io, ip ) = mels( io, ip ) + 14.0d0 / 45.0d0 * al * rr( i1 ) * ( ff( 0 ) + ff( 4 ) ) &
                      + 64.0d0 / 45.0d0 * al * rr( i1 ) * ( ff( 1 ) + ff( 3 ) ) &
                      + 24.0d0 / 45.0d0 * al * rr( i1 ) * ff( 2 )
     end do
   end do
 end do
 !



end subroutine getmeznl
