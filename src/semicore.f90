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
 subroutine semicore( zz, nc, nn, ll, la, lmax, irc, nopf, maxopf, rr, coreuu, aepr, mels )
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
 real(dp),intent(in) :: zz
 integer, intent(in) :: nc,nn(nc),ll(nc),la,lmax,irc,nopf, maxopf
 real(dp),intent(in) :: rr(irc),coreuu(irc,nc), aepr(irc,nopf), mels(maxopf,0:3,0:lmax,0:lmax)


!Local variables

 real(dp) :: al,ff(0:4),rtmp
 real(dp),allocatable :: semi(:)
 real(dp),allocatable :: rmel(:,:)
 integer :: ls, nmax
 integer :: ip, io, i1, ii, jj, npowr, ic

 character(len=14) :: filnam14 
 character(len=17) :: filnam
 character(len=1) :: hashtag

 logical :: ex

 ! Would need to fix the write statement below as well
 npowr = 3

 allocate(rmel(0:npowr,nc),semi(irc))
 rmel(:,:) = 0.0_dp
 semi(:) = 0.0_DP
 al = 0.01d0 * dlog(rr(101)/rr(1))

 ! try a semi-core for each l
 do ls = 0, lmax
   nmax = ls
   do ic = 1, nc
     if( ll(ic) .eq. ls ) then
       nmax = max( nmax, nn(ic) )
     endif
   enddo
   nmax = nmax + 1

   write( filnam, '(A8,1i3.3,2(1a1,1i2.2))' ) 'semiorbz', nint( zz ), 'n', nmax, 'l', ls
   write(6,*) filnam

   inquire(file=filnam, exist=ex )
   if( .not. ex ) cycle 

   open(unit=99, file=filnam, form='formatted', status='old' )
   read(99,*) hashtag, jj
   semi(:) = 0.0_DP
   do ii = 1, min(irc,jj)
     read(99,*) rtmp, semi(ii)
   enddo
   close(99)

   do io = 1, nc
     do ip = 0, npowr
       do i1 = 1, irc - 4, 4
         do jj = 0, 4
           ii = i1 + jj
           ff( jj ) = semi( ii ) * coreuu( ii, io ) * rr( ii ) ** ip
         enddo

         rmel( ip, io ) = rmel( ip, io ) + 14.0d0 / 45.0d0 * al * rr( i1 ) * ( ff( 0 ) + ff( 4 ) ) &
                        + 64.0d0 / 45.0d0 * al * rr( i1 ) * ( ff( 1 ) + ff( 3 ) ) &
                        + 24.0d0 / 45.0d0 * al * rr( i1 ) * ff( 2 )
       end do
     end do
   end do
   !

   ! output result
   write ( filnam14, '(1a3,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 's2c', 'z', nint( zz ), 'n', nmax, 'l', ls
   open( unit=99, file=filnam14, form='formatted', status='unknown' )
   rewind 99
   write ( 99, '(1i5)' ) npowr
   do ii = 1, nc
       write ( 99, '(4(1x,1e15.8),2(1x,I3))' ) rmel( :, ii), nn(ii), ll(ii)
   end do
   close( unit=99 )

   call getmeznl( zz, nmax, ls, irc, nopf, maxopf, rr, semi(:), aepr, mels(:,:,la,ls) )
 enddo
 deallocate(rmel,semi)

 return

end subroutine semicore
