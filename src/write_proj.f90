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
 subroutine write_proj( irc, mmax, nopf, zz, ll, rr, pspr, aepr )
! Writes out the projectors for a particular angular momentum

!mmax  size of radial grid
!irc  size of the projectors
!nopf number of optimal projectors
!ll angular momentum l
!zz  atomic number
!rr  log radial grid
!pspr pseduo opfs
!aepr all-electron opfs


 implicit none
 integer, parameter :: dp=kind(1.0d0)


!Input variables
 integer, intent(in) :: irc,mmax,nopf,ll
 real(dp),intent(in) :: zz
 real(dp),intent(in) :: rr(mmax)
 real(dp),intent(in) :: pspr(irc,nopf),aepr(irc,nopf)

!Output variables - printing only

!Local variables
 integer :: nq, ii, jj, ir, i1
 real(dp) :: qman, dq, qq, arg, wr, jl, su(nopf), dl, qmax
 character(len=16) :: fmtString
 character(len=8) :: nam
 character(len=11) :: nam2

 dl = 0.01d0 * dlog(rr(101)/rr(1))
 qmax = 20.0_dp
 dq = 0.05_dp
 nq = 1 + ceiling( qmax / dq )
 if( nopf .gt. 99 ) then
   write(6,*) 'NOT PROGRAMMED FOR 100+ OPFs'
   return
 endif

 write(nam2,'(a,I3.3)') 'radfilez', nint( zz )
 open(unit=99, file=nam2, form='formatted', status='unknown' )
 write(99,*) rr(irc), mmax, irc
 close( 99 )

 write( fmtString, '(A,I0,A)' ) '(', nopf+1, '(1x,1e22.15))'
 write ( unit=nam, fmt='(1a2,1i1,1a1,1i3.3)' ) 'ae', ll, 'z', nint( zz )
 open( unit=99, file=nam, form='formatted', status='unknown' )
 rewind 99
 do ii = 1, irc
!    write ( 99, '(10(1x,1e22.15))' ) r( i ), proj( i, : ) / r( i )
    write ( 99, fmtString ) rr( ii ), aepr( ii, : ) / rr( ii )
 end do
 close( unit=99 )
 !
 write ( unit=nam, fmt='(1a2,1i1,1a1,1i3.3)' ) 'ps', ll, 'z', nint( zz )
 open( unit=99, file=nam, form='formatted', status='unknown' )
 rewind 99
 do ii = 1, irc
    write ( 99, fmtString ) rr( ii ), pspr( ii, : ) / rr( ii )
 end do
 close( unit=99 )
 !
 write ( unit=nam, fmt='(1a2,1i1,1a1,1i3.3)' ) 'ft', ll, 'z', nint( zz )
 open( unit=99, file=nam, form='formatted', status='unknown' )
 rewind 99
 do jj = 0, nq - 1
    qq = dq * real( jj, dp )
    su( : ) = 0
    do ir = 1, irc
       ii = ir - 1; i1 = ii / 4; i1 = ii - i1 * 4 ! use bode rule
       if ( i1 .eq. 0 ) wr = 28.0d0
       if ( i1 .eq. 1 ) wr = 64.0d0
       if ( i1 .eq. 2 ) wr = 24.0d0
       if ( i1 .eq. 3 ) wr = 64.0d0
       if ( ( ir .eq. 1 ) .or. ( ir .eq. irc ) ) wr = 14.0d0 ! with endpoints fixed
       arg = qq * rr( ir )
       if ( arg .lt. 0.000001d0 ) then
          if ( ll .eq. 0 ) jl = 1 - arg ** 2 / 6
          if ( ll .eq. 1 ) jl = arg / 3 - arg ** 3 / 30
          if ( ll .eq. 2 ) jl = arg ** 2 / 15 - arg ** 4 / 210
          if ( ll .eq. 3 ) jl = arg ** 3 / 105 - arg ** 5 / 1890
       else
          if ( ll .eq. 0 ) jl = sin( arg ) / arg
          if ( ll .eq. 1 ) jl = sin( arg ) / arg ** 2 - cos( arg ) / arg
          if ( ll .eq. 2 ) jl = 3 * ( sin( arg ) / arg ** 3 - cos( arg ) / arg ** 2 ) - sin( arg ) / arg
          if ( ll .eq. 3 ) jl = 15 * ( sin ( arg ) / arg ** 4 - cos( arg ) / arg ** 3 ) &
                             - 6 * sin( arg ) / arg ** 2 + cos( arg ) / arg
       end if
       su( : ) = su( : ) + wr * rr( ir ) ** 2 * pspr( ir, : ) * jl
    end do
    write ( 99, '(1x,10(1x,1e15.8))' ) qq, su( : ) * dl / 45.0d0
 end do
 close( unit=99 )
 !
 return
end subroutine write_proj
