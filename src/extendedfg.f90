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
 subroutine extendedfg( zz, nc, na, la, lmax_opf, maxopf, nopfs, mmax, irc, rr, coreuu, aepr )

!zz  atomic number
!nn principle quantum number
!lc angular momentum of core
!ll angular momentum l of opfs
!irc  size of the projectors
!nopf number of optimal projectors
!rr  log radial grid
!coreuu core level orbital
!aepr all-electron opfs

 implicit none
 integer, parameter :: dp=kind(1.0d0)


!Input variables
  real(dp),intent(in) :: zz
 integer, intent(in) :: nc, na(nc),la(nc)
 integer, intent(in) :: lmax_opf, maxopf, nopfs(0:lmax_opf),mmax,irc
 real(dp),intent(in) :: rr(mmax),coreuu(mmax,nc)
 real(dp),intent(in) :: aepr(irc,maxopf,0:lmax_opf)

!Output variables 

!Local variables
 integer :: kk,ii,i1,i2,k1,ic,l1,l2,lc,lmaxf,lmaxg
 real(dp) :: tmp,s11,s23,s12,s13,dl,rp
 real(dp),allocatable :: ff(:,:),gg(:,:)
 real(dp),allocatable :: v11(:),v12(:),v23(:),v13(:)
 character(len=18) :: filnam18

 logical :: dogg, doff

! Remove later and write these in Ha. to leave all eV conversion in ocean.x
 real(dp),parameter :: ehart = 27.21138506_DP

! allocate(ff(nopf,nopf),gg(nopf,nopf),phv(irc,nopf) )
 allocate(v11(irc),v12(irc),v13(irc),v23(irc))

 dl = 0.01d0 * dlog(rr(101)/rr(1))

! do i2 = 1, nopf
!   phv(:,i2) = aepr(:,i2) !* rr( : )
! enddo

! do kk = 0, 2 * max( lc, ll )

! Probably could be re-arranged to be better about re-using things
  do ic = 1, nc
    lc = la(ic)
    do l1 = 0, lmax_opf
      do l2 = 0, l1
        ! For exchange, the max is the smaller of (lc+l1),(lc+l2)
        ! For direct, the max is the smaller of 2*lc,l1+l2
        lmaxf = min( 2 * lc, l1+l2)
        lmaxg = min( lc+l1, lc+l2)
        do kk = 0, max(lmaxf,lmaxg)
!        do kk = 0, max(2*lc, 2*lmax_opf)
          ! Enforce 3j selection rules
          ! 1. |lc-l| <= k
          ! 2. k <= (lc+l)
          ! 3. lc+l+k = even integer
          dogg = .true.
          doff = .true.

          if( ( abs(lc-l1) .gt. kk ) .or. ( abs(lc-l2) .gt. kk ) .or. (lc+l1 .lt. kk ) & 
             .or. (lc+l2 .lt. kk ) .or. ( mod(lc+l1+kk,2) .ne. 0 ) .or. ( mod(lc+l2+kk,2) .ne. 0 ) ) then
            dogg = .false.
          endif
          if( ( abs(l1-l2) .gt. kk ) .or. ( l1+l2 .lt. kk ) .or. (lc+lc .lt. kk )  &
            .or. (mod(kk,2) .ne. 0 ) .or. mod(l1+l2,2) .ne. 0 ) then
            doff = .false.
          endif
          write(6,*) 'ZZZZ', lc, l1, l2, kk, dogg, doff
          if( (.not. doff ) .and. ( .not. dogg ) ) cycle

          allocate( ff(nopfs(l1),nopfs(l2)), gg(nopfs(l1),nopfs(l2)) )

          do i1 = 1, nopfs(l1)
            do i2 = 1, nopfs(l2)
              !
              v11( : ) = 0.0_dp; v23( : ) = 0.0_dp; v12( : ) = 0.0_dp; v13( : ) = 0.0_dp
              s11 = 0; s23 = 0; s12 = 0; s13 = 0
              do ii = irc - 1, 1, -1
                tmp = 0.5_dp * dl * rr( ii ) / rr( ii ) ** ( kk + 1 )
                s11 = s11 + tmp * coreuu( ii, ic ) * coreuu( ii, ic )
                s23 = s23 + tmp * aepr( ii, i1, l1 ) * aepr( ii, i2, l2 )
                s12 = s12 + tmp * coreuu( ii, ic ) * aepr( ii, i1, l1 )
                s13 = s13 + tmp * coreuu( ii, ic ) * aepr( ii, i2, l2 )
                tmp = 0.5_dp * dl * rr( ii + 1 ) / rr( ii + 1 ) ** ( kk + 1 )
                s11 = s11 + tmp * coreuu( ii + 1, ic ) * coreuu( ii + 1, ic )
                s23 = s23 + tmp * aepr( ii + 1, i1, l1 ) * aepr( ii + 1, i2, l2 )
                s12 = s12 + tmp * coreuu( ii + 1, ic ) * aepr( ii + 1, i1, l1 )
                s13 = s13 + tmp * coreuu( ii + 1, ic ) * aepr( ii + 1, i2, l2 )
                rp = rr( ii ) ** kk
                if ( rp .gt. 0.0_dp ) then
                  v11( ii ) = s11 * rp
                  v23( ii ) = s23 * rp
                  v12( ii ) = s12 * rp
                  v13( ii ) = s13 * rp
                end if
              end do
              s11 = 0; s23 = 0; s12 = 0; s13 = 0
              do ii = 2, irc
                tmp = 0.5_dp * dl * rr( ii - 1 ) * rr( ii - 1 ) ** kk
                s11 = s11 + tmp * coreuu( ii - 1, ic ) * coreuu( ii - 1, ic )
                s23 = s23 + tmp * aepr( ii - 1, i1, l1 ) * aepr( ii - 1, i2, l2 )
                s12 = s12 + tmp * coreuu( ii - 1, ic ) * aepr( ii - 1, i1, l1 )
                s13 = s13 + tmp * coreuu( ii - 1, ic ) * aepr( ii - 1, i2, l2 )
                tmp = 0.5_dp * dl * rr( ii ) * rr( ii ) ** kk
                s11 = s11 + tmp * coreuu( ii, ic ) * coreuu( ii, ic )
                s23 = s23 + tmp * aepr( ii, i1, l2 ) * aepr( ii, i2, l2 )
                s12 = s12 + tmp * coreuu( ii, ic ) * aepr( ii, i1, l1 )
                s13 = s13 + tmp * coreuu( ii, ic ) * aepr( ii, i2, l2 )
                rp = rr( ii ) ** ( kk + 1 )
                if ( rp .gt. 0.0_dp ) then
                  v11( ii ) = v11( ii ) + s11 / rp
                  v23( ii ) = v23( ii ) + s23 / rp
                  v12( ii ) = v12( ii ) + s12 / rp
                  v13( ii ) = v13( ii ) + s13 / rp
                end if
              end do
              !
              s11 = 0; s23 = 0; s12 = 0; s13 = 0
              do ii = 1, irc
                tmp = dl * rr( ii )
                if ( ( ii .eq. 1 ) .or. ( ii .eq. irc ) ) tmp = 0.5_dp * tmp
                s11 = s11 + tmp * v23( ii ) * coreuu( ii, ic ) * coreuu( ii, ic )
                s23 = s23 + tmp * v11( ii ) * aepr( ii, i1, l1 ) * aepr( ii, i2, l2 )
                s12 = s12 + tmp * v13( ii ) * coreuu( ii, ic ) * aepr( ii, i1, l1 )
                s13 = s13 + tmp * v12( ii ) * coreuu( ii, ic ) * aepr( ii, i2, l2 )
              end do
              !
              ff( i1, i2 ) = 0.5_dp * ( s11 + s23 ) * ehart
              gg( i1, i2 ) = 0.5_dp * ( s12 + s13 ) * ehart
              !
            end do
          end do

          if( dogg ) then
            write ( filnam18, '(1a2,4i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'gk', lc, l1, l2, kk, 'z', nint( zz ), 'n', na(ic), 'l', la(ic)
            open( unit=99, file=filnam18, form='formatted', status='unknown' )
            rewind 99
            do i2 = 1, nopfs(l2)
               write ( 99, '(9f10.4)' ) gg( :, i2 )
            end do
            write(99,*) 1.0_DP
            close( unit=99 )
          endif

          if( doff ) then
            write ( filnam18, '(1a2,4i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'fk', lc, l1, l2, kk, 'z', nint( zz ), 'n', na(ic), 'l', la(ic)
            open( unit=99, file=filnam18, form='formatted', status='unknown' )
            rewind 99
            do i2 = 1, nopfs(l2)
               write ( 99, '(9f10.4)' ) ff( :, i2 )
            end do
            write(99,*) 1.0_DP
            close( unit=99 )
          endif
      
          deallocate( ff, gg )
        enddo ! l2
      enddo !l1
    enddo !kk
  enddo ! ic

 deallocate(v11,v12,v13,v23)

end subroutine extendedfg
