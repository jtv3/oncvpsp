! John Vinson
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

subroutine corezeta( rr,zz,mmax,nc,nv,na,la,fa,iexc )


 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: Ryd2eV = 13.605693122994_DP

 real(dp),intent(in) :: zz, rr(mmax)
 integer,intent(in) :: mmax, nc, nv, iexc
 integer,intent(in) ::  na(nv+nc), la(nv+nc)
 real(dp),intent(in) :: fa(nc+nv)


 real(dp) :: ea(30,2), rpk(30,2), etot
 integer :: it, ierr
 real(dp), allocatable :: rho(:), rhoc(:), vi(:)

 real(dp) :: e1, e2
 integer :: ii
 character(len=1) :: symb(0:3) = (/ 's', 'p', 'd', 'f' /)


 allocate( rho(mmax), rhoc(mmax), vi(mmax) )
 
 call relatom(na,la,ea,fa,rpk,nc, nc+nv, it, rhoc, rho, &
              rr, vi, zz, mmax, iexc, etot, ierr )

 open(unit=99,file='xifile',form='formatted',status='unknown')
 do ii = 1, nc+nv
   if( la(ii) .eq. 0 ) then
     e1 = 0.0_DP
   else
     e1 = ea(ii,1)-ea(ii,2)
   endif
   e1 = e1 * 2.0_dp / dble( 2 * la( ii ) + 1 )
   e2 = e1 * Ryd2eV * 2.0_dp !* 4.0d0 / dble( 2 * la( ii ) + 1 )
   write(99,'(F20.11,X,F20.11,XI1,A1)') e2, e1, na(ii), symb(la(ii))
   write(6,'(A,X,F20.11,X,F20.11,XI1,A1)') 'OCEAN SO', e2, e1, na(ii), symb(la(ii))
 enddo
 close(99)

 deallocate( rho, rhoc, vi )

end subroutine corezeta
