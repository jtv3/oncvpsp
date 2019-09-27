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
!
! Give the various matrix elements as a function of cur-off radius
!
 subroutine diagfg( zz, nn, lc, ll, irc, nopf,  rr, coreuu, aepr  )

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
 integer, intent(in) :: nn,ll,lc,irc,nopf
 real(dp),intent(in) :: zz
 real(dp),intent(in) :: rr(irc),coreuu(irc)
 real(dp),intent(in) :: aepr(irc,nopf)

!Output variables 

!Local variables
 integer :: kk,ii,i2,i3,k1,nud,imap,jrc
 real(dp) :: tmp,s11,s23,s12,s13,dl,rp
 real(dp),allocatable :: ff(:,:),gg(:,:),phv(:,:)
 real(dp),allocatable :: afk(:,:),agk(:,:)
 real(dp),allocatable :: v11(:),v12(:),v23(:),v13(:)
 character(len=18) :: filnam18
 character(len=100) :: fmtstm

! Remove later and write these in Ha. to leave all eV conversion in ocean.x
 real(dp),parameter :: ehart = 27.21138506_DP

 allocate(ff(nopf,nopf),gg(nopf,nopf),phv(irc,nopf) )
 allocate(v11(irc),v12(irc),v13(irc),v23(irc))

 dl = 0.01d0 * dlog(rr(101)/rr(1))

 do i2 = 1, nopf
   phv(:,i2) = aepr(:,i2) !* rr( : )
 enddo

 nud = nopf * ( nopf + 1 )
 nud = nud / 2

 write(6,*) nud, irc
 allocate( afk( nud, irc ), agk( nud, irc ) )
 do kk = 0, 2 * max( lc, ll )

   do jrc = 1, irc
   imap = 0

   do i2 = 1, nopf
     do i3 = 1, nopf
       !
       v11( : ) = 0.0_dp; v23( : ) = 0.0_dp; v12( : ) = 0.0_dp; v13( : ) = 0.0_dp
       s11 = 0; s23 = 0; s12 = 0; s13 = 0
       do ii = jrc - 1, 1, -1
         tmp = 0.5_dp * dl * rr( ii ) / rr( ii ) ** ( kk + 1 )
         s11 = s11 + tmp * coreuu( ii ) * coreuu( ii )
         s23 = s23 + tmp * phv( ii, i2 ) * phv( ii, i3 )
         s12 = s12 + tmp * coreuu( ii ) * phv( ii, i2 )
         s13 = s13 + tmp * coreuu( ii ) * phv( ii, i3 )
         tmp = 0.5_dp * dl * rr( ii + 1 ) / rr( ii + 1 ) ** ( kk + 1 )
         s11 = s11 + tmp * coreuu( ii + 1 ) * coreuu( ii + 1 )
         s23 = s23 + tmp * phv( ii + 1, i2 ) * phv( ii + 1, i3 )
         s12 = s12 + tmp * coreuu( ii + 1 ) * phv( ii + 1, i2 )
         s13 = s13 + tmp * coreuu( ii + 1 ) * phv( ii + 1, i3 )
         rp = rr( ii ) ** kk
         if ( rp .gt. 0.0_dp ) then
           v11( ii ) = s11 * rp
           v23( ii ) = s23 * rp
           v12( ii ) = s12 * rp
           v13( ii ) = s13 * rp
         end if
       end do
       s11 = 0; s23 = 0; s12 = 0; s13 = 0
       do ii = 2, jrc
         tmp = 0.5_dp * dl * rr( ii - 1 ) * rr( ii - 1 ) ** kk
         s11 = s11 + tmp * coreuu( ii - 1 ) * coreuu( ii - 1 )
         s23 = s23 + tmp * phv( ii - 1, i2 ) * phv( ii - 1, i3 )
         s12 = s12 + tmp * coreuu( ii - 1 ) * phv( ii - 1, i2 )
         s13 = s13 + tmp * coreuu( ii - 1 ) * phv( ii - 1, i3 )
         tmp = 0.5_dp * dl * rr( ii ) * rr( ii ) ** kk
         s11 = s11 + tmp * coreuu( ii ) * coreuu( ii )
         s23 = s23 + tmp * phv( ii, i2 ) * phv( ii, i3 )
         s12 = s12 + tmp * coreuu( ii ) * phv( ii, i2 )
         s13 = s13 + tmp * coreuu( ii ) * phv( ii, i3 )
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
!       if( i2 .ge. i3 ) imap = imap + 1
       do ii = 1, jrc
         tmp = dl * rr( ii )
         if ( ( ii .eq. 1 ) .or. ( ii .eq. jrc ) ) tmp = 0.5_dp * tmp
         s11 = s11 + tmp * v23( ii ) * coreuu( ii ) * coreuu( ii )
         s23 = s23 + tmp * v11( ii ) * phv( ii, i2 ) * phv( ii, i3 )
         s12 = s12 + tmp * v13( ii ) * coreuu( ii ) * phv( ii, i2 )
         s13 = s13 + tmp * v12( ii ) * coreuu( ii ) * phv( ii, i3 )
!         if( i2 .ge. i3 ) then
!           afk( imap, ii ) = 0.5_dp * ( s11 + s23 ) * ehart
!           agk( imap, ii ) = 0.5_dp * ( s12 + s13 ) * ehart
!         endif
       end do
       !
       ff( i2, i3 ) = 0.5_dp * ( s11 + s23 ) * ehart
       gg( i2, i3 ) = 0.5_dp * ( s12 + s13 ) * ehart
       if( i2 .ge. i3 ) then
         imap = imap + 1
         afk( imap, jrc ) = 0.5_dp * ( s11 + s23 ) * ehart
         agk( imap, jrc ) = 0.5_dp * ( s12 + s13 ) * ehart
       endif

       !
     end do
   end do
   enddo

   if( mod(kk,2) .eq. 0 .and. kk .le. 2*min(ll,lc) ) then
     write ( filnam18, '(1a2,3i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'fk', lc, ll, kk, 'z', nint( zz ), 'n', nn, 'l', lc
     open( unit=99, file=filnam18, form='formatted', status='unknown' )
     rewind 99
     do i2 = 1, nopf
        write ( 99, '(9f8.2)' ) ff( :, i2 )
     end do
     write(99,*) 0.80
     close( unit=99 )

     write( filnam18, '(1a2,5i1,a)') 'fk', lc, ll, kk, nn, lc, 'diagnose' 
     open( unit=99, file=filnam18, form='formatted', status='unknown' )
     rewind 99
     ! write the stupid format statement to be adjustable
     write( fmtstm, '(A,I0,A)' ) '(f20.14,', nud, 'f12.4)' 
     do ii = 1, irc
       write( 99, fmtstm ) rr(ii), afk( :, ii )
     enddo
     close( 99 )
   end if

   do k1 = abs( ll - lc ), ll + lc, 2
     if ( k1 .eq. kk ) then
       write ( filnam18, '(1a2,3i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'gk', lc, ll, kk, 'z', nint( zz ), 'n', nn, 'l', lc
       open( unit=99, file=filnam18, form='formatted', status='unknown' )
       rewind 99
       do i2 = 1, nopf
          write ( 99, '(9f8.2)' ) gg( :, i2 )
       end do
       write(99,*) 0.80
       close( unit=99 )

       write( filnam18, '(1a2,5i1,a)') 'gk', lc, ll, kk, nn, lc, 'diagnose' 
       open( unit=99, file=filnam18, form='formatted', status='unknown' )
       rewind 99
       ! write the stupid format statement to be adjustable
       write( fmtstm, '(A,I0,A)' ) '(f20.14,', nud, 'f12.4)'
       do ii = 1, irc
         write( 99, fmtstm ) rr(ii), agk( :, ii )
       enddo
       close( 99 )
     end if
   enddo

 end do

 deallocate( agk, afk )
 deallocate(ff,gg,phv,v11,v12,v13,v23)

end subroutine diagfg
