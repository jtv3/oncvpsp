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
 subroutine core2core( nc, ic, zz, nn, ll, mmax, irc, rr, coreuu )
! Writes out the matrix elements between core orbitals for 
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
 integer, intent(in) :: nc,ic,nn(nc),ll(nc),irc,mmax
 real(dp),intent(in) :: zz
 real(dp),intent(in) :: rr(irc),coreuu(mmax,nc)


!Local variables

 real(dp) :: al,ff(0:4)
 real(dp),allocatable :: rmel(:,:)
 integer :: ip, io, i1, ii, jj, npowr

 character(len=14) :: filnam14

 ! Would need to fix the write statement below as well
 npowr = 3

 allocate(rmel(0:npowr,ic))
 rmel(:,:) = 0.0_dp
 al = 0.01d0 * dlog(rr(101)/rr(1))

 do io = 1, ic
   do ip = 0, npowr
     do i1 = 1, irc - 4, 4
       do jj = 0, 4
         ii = i1 + jj
         ff( jj ) = coreuu( ii, ic ) * coreuu( ii, io ) * rr( ii ) ** ip
       enddo

       rmel( ip, io ) = rmel( ip, io ) + 14.0d0 / 45.0d0 * al * rr( i1 ) * ( ff( 0 ) + ff( 4 ) ) &
                      + 64.0d0 / 45.0d0 * al * rr( i1 ) * ( ff( 1 ) + ff( 3 ) ) &
                      + 24.0d0 / 45.0d0 * al * rr( i1 ) * ff( 2 )
     end do
   end do
 end do
 !

 ! output result
  write ( filnam14, '(1a3,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'c2c', 'z', nint( zz ), 'n', nn(ic), 'l', ll(ic)
  open( unit=99, file=filnam14, form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i5)' ) npowr
  do ii = 1, ic
      write ( 99, '(4(1x,1e15.8),2(1x,I3))' ) rmel( :, ii), nn(ii), ll(ii)
  end do
  close( unit=99 )

  deallocate(rmel)


end subroutine core2core
