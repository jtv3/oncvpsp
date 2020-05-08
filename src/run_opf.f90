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
 subroutine run_opf(lmax,lloc,nproj,ep,epsh1,epsh2,depsh,vkb,evkb, &
&                     rr,vfull,vp,zz,mmax,irc,srel,nc,na,la,ea)

! computes the OPFs and writes them to file

!lmax  maximum angular momentum
!lloc  l for local potential
!nproj  number ov V / KB projectors for  each l
!ep  bound-state or scattering state reference energies for vkb potentials
!epsh1  low energy limit for "phase shift" calculation
!epsh2  high energy limit for "phase shift" calculation
!depsh  energy increment
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!rr  log radial grid
!vfull  all-electron potential
!vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
!zz  atomic number
!mmax  size of radial grid
!irphs  index of rr beyond which all vp==vlocal
! irc
!srel .true. for scalar-relativistic, .false. for non-relativistic
!nc number of core levels
!na principle quantum numbers of core levels
!la angular quantum numbers of core levels

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer,intent(in) :: lmax,lloc,mmax,irc(6)
 integer,intent(in) :: nproj(6)
 real(dp) :: epsh1,epsh2,depsh,zz
 real(dp),intent(in) :: rr(mmax),vp(mmax,5),ep(6)
 real(dp),intent(in) :: vfull(mmax),vkb(mmax,2,4),evkb(2,4)
 logical,intent(in) :: srel
 integer,intent(in) :: nc
 integer,intent(in) :: na(nc),la(nc)
 real(dp),intent(in) :: ea(nc)

!Output variables - printing only

!Local variables
 integer :: ii,ll,l1,npsh,jj,imax,imin,j1,j2,kk,ierr,nopf,mch,ic,l2,iii
 real(dp) :: epsh,singleps,e1,e2,a1,a2,eebest,ee,sca,abest,coree
 logical :: qual
 character(len=20 ) :: fnam

 real(dp),allocatable :: pshf(:),pshp(:),aeuu(:),aeup(:),psuu(:),psup(:),coreuu(:,:),coreup(:)
 real(dp),allocatable :: phips(:,:), phirn(:,:),pspr(:,:),aepr(:,:)
 real(dp),allocatable :: mels(:,:,:,:)
 integer, allocatable :: nopfs(:)

 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: prec = 0.001_dp
 real(dp), parameter :: tol = 0.001_dp
 integer,parameter :: maxopf = 16

 integer :: irphs
 real(dp) :: targRad
 logical :: ex

 irphs = maxval( irc + 2 )

 inquire( file='overrideRadius', exist=ex )
 if( ex ) then
   open( unit=99, file='overrideRadius', form='formatted', status='old' )
   read( 99, * ) targRad
   close( 99 )
   do ii = 2, mmax
     if( rr( ii ) .gt. targRad ) then
       if( ( targRad - rr( ii - 1 ) ) .lt. ( rr( ii ) - targRad ) ) then
         irphs = ii - 1
       else
         irphs = ii
       endif
       exit !quit loop, best radius found
     endif
   enddo
   write( 6, * ) 'OCEAN: Using modified radius: ', rr( irphs ), rr( maxval( irc ) )
 endif

! loop for phase shift calculation -- full, then local or Kleinman-
! Bylander / Vanderbilt
 
 write(6,*) 'OCEAN: OPF SECTION', lmax
 npsh = 128
 allocate(pshf(npsh),pshp(npsh))
 allocate(phips(irphs,npsh),phirn(irphs,npsh),pspr(irphs,npsh),aepr(irphs,npsh))
 allocate(aeuu(mmax),aeup(mmax),psuu(mmax),psup(mmax),coreuu(mmax,nc),coreup(mmax))

 allocate( mels(maxopf,0:3,0:lmax,nc), nopfs(0:lmax) )

 epsh2 = 5.0_DP
 do l1 = 1, lmax+1
   epsh1 = ep( l1 ) - 0.3_dp
   depsh = (epsh2-epsh1)/real(npsh,dp)
   write(6,'(A7,I0,3(X,3E16.8))') 'OCEAN: ', l1-1, epsh1, epsh2, depsh

   ll = l1 - 1

   ! Only need the phase shifts for the all-electron to start
   call fphsft(ll,epsh2,depsh,pshf,rr,vfull,zz,mmax,irphs,npsh,srel)

   ! TESTING
   if(ll .eq. lloc) then  
     call  vkbphsft(ll,0,epsh2,depsh,ep(l1),pshf,pshp, &
&                   rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                   mmax,irphs,npsh)
   else
     call  vkbphsft(ll,nproj(l1),epsh2,depsh,ep(l1),pshf,pshp, &
&                   rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                   mmax,irphs,npsh)
   end if
   !\ TESTING

   if( abs( pshp( npsh ) - pshf( npsh ) ) > 2.0_dp ) then
     jj = nint( ( pshp( npsh ) - pshf( npsh ) ) / pi )
     pshp( : ) = pshp( : ) - pi * jj
   endif

   ! Now loop on energies. If the phase shift of the pseudo at this energy
   ! is greater than pshf(1) cycle

   ! The ordering of the energies is from highest to lowest
   imin = 1
   do ii = 1, npsh
     if( pshp( ii ) .lt. pshf( npsh ) ) imax = ii
   end do
   do ii = npsh, 1, -1
     if( pshp( ii ) .gt. pshf( 1 ) ) imin = ii
   end do
   write(6,*) 'OCEAN: ', pshf(1), pshf(npsh)
   write(6,*) 'OCEAN: ', pshp(1), pshp(npsh)
   write(6,'(A7,3(X,I0))') 'OCEAN: ', imin, imax, npsh
   
   write(fnam,  '(A,I3.3,A,I1.1,A)') 'z', nint( zz ), 'l', ll, '.scat' 
   open(unit=98,file=fnam)

   do ii = imin, imax
     epsh = epsh2 - depsh * dfloat(ii - 1)

!     if( pshp( ii ) .gt. pshf( 1 ) ) cycle
!     if( pshp( ii ) .lt. pshf( npsh ) ) cycle

     do jj = npsh, 1, -1
       if( pshf( jj ) .gt. pshp( ii ) ) j2 = jj
     end do
     e2 = epsh2 - depsh * dfloat(j2 - 1)
     a2 = pshf( j2 )
     do jj = 1, npsh
       if( pshf( jj ) .lt. pshp( ii ) ) j1 = jj
     end do
     e1 = epsh2 - depsh * dfloat(j1 - 1)
     a1 = pshf( j1 )
!     if( ii .eq. imin ) then
!       write(6,*) ii, j1, j2
!       write(6,*) pshp( ii ), a1, a2
!     endif

     ! refine search to use all-electron wave that matches pseudo version
!     write(6,'(X,7(3X,A12))') 'diff', 'a1', 'a2', 'abest', 'e1', 'e2', 'ee' 


     do kk = 1, 3
       ee = ( e1 + e2 ) * .5_dp
       ! get all-electron phase at this energy
       call fphsft(ll,ee,depsh,singleps,rr,vfull,zz,mmax,irphs,1,srel)
       !JTV do a better job aligning phase
       do while ( singleps .lt. ( a1 - pi ) )
         singleps = singleps + 2.0_dp * pi
         if( singleps .gt. (a2+pi/2.0_dp) ) write(6,*) 'OCEAN: ERROR'
       enddo

       if( ( singleps - pshp( ii ) ) * ( pshp( ii ) - a2 ) .gt. 0.0_dp ) then
         a1 = singleps
         e1 = ee
       else
         a2 = singleps
         e2 = ee
       end if
     enddo

     ! do a better check and exit here
     ! allow for more iterations but include an escape check when (pshp(ii) - singleps) < 1.0d-10 ?
     do kk = 1, 9
       ee = ( e1 * ( a2 - pshp( ii ) ) + e2 * ( pshp( ii ) - a1 ) ) / ( a2 - a1 )
       ! get all-electron phase at this energy
       call fphsft(ll,ee,depsh,singleps,rr,vfull,zz,mmax,irphs,1,srel)
       !JTV do a better job aligning phase
       do while ( singleps .lt. ( a1 - pi ) )
         singleps = singleps + 2.0_dp * pi
         if( singleps .gt. (a2+pi/2.0_dp) ) write(6,*) 'OCEAN: ERROR'
       enddo

       if( kk .eq. 1 ) then
         qual = .true.
       else
         qual = ( abs( singleps - pshp( ii ) ) .lt. abs( abest - pshp( ii ) ) )
       end if
       if( qual ) then
         eebest = ee
         abest = singleps
         if( abs( abest - pshp( ii ) ) .lt. 3.0d-14 ) then
           exit
         endif
       end if
!       if( ii .eq. imin ) write(6,'(I0,7(3X,E12.6))') kk, abs(abest - pshp( ii )), a1, a2, abest, e1, e2, ee
!       write(6,'(I1,7(3X,E12.5))') kk, abs(abest - pshp( ii )), a1, a2, abest, e1, e2, ee
       if( ( singleps - pshp( ii ) ) * ( pshp( ii ) - a2 ) .gt. 0.0_dp ) then
         a1 = singleps
         e1 = ee
       else
         a2 = singleps
         e2 = ee
       end if
     end do 
!     write(6,*) 'BBB', pshp( ii )-abest, kk
         
     ! Now go ahead and get orbitals
     call lschfs(ll,ierr,eebest,rr,vfull,aeuu,aeup,zz,mmax,irphs,srel)
!     call lschfs(ll,ierr,epsh,rr,vfull,aeuu,aeup,zz,mmax,irphs,srel)
     call lschvkbs(ll,nproj(l1),epsh,rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1),psuu,psup,mmax,irphs)

     ! Because of number of node mismatches, there may be a sign difference in the derivatives of psi
     ! between the pseudo and all-electron, but sca should be ~+1 so go ahead and take the abs( )
     sca = abs( ( pshp( ii )**2 + psup(irphs)**2 ) / ( pshp( ii ) * abest + psup(irphs) * aeup(irphs) ) )
!!     if( ii .eq. imin ) write(6,'(A,5(X,E12.6))') 'SCA', sca, abest, pshp( ii ), psup(irphs), aeup(irphs)
!      write(6,'(A,1(X,E12.5))') 'SCA', sca
!!     sca = 1.0_dp

     ! Apparently need to check to make sure the signs are in agreement. Might be missing node in search range?
     do iii = irphs, 1, -1
       if( abs( aeuu( iii ) ) .gt. tol .and. abs( psuu( iii ) ) .gt. tol ) then
         if( aeuu( iii ) * psuu( iii ) .lt. 0 ) sca = -sca
          exit
       endif
     enddo

     phirn( :, ii - imin + 1) = aeuu( : ) * sca
     phips( :, ii - imin + 1) = psuu( : )

     write(fnam, '(A,I3.3,A,I1.1,A,I4.4)') 'z', nint( zz ), 'l', ll, '.', imax - ii + 1
     open(file=fnam, unit=99, form='formatted', status='unknown' )
     write(99,*) '#', epsh
     write(98,*)
     do jj = 1, irphs
       write(99,*) rr(jj), phips(jj,ii), phirn(jj,ii)
       write(98,*) rr(jj), phips(jj,ii), phirn(jj,ii)
     end do
     write(98,*)
     close(99)

   end do

   close( 98 )


   ! Principle-component analysis of the set of scattering waves (N = imax-imin+1 )
   ! prec sets the eigenvalue threshold. Should be based on basis size?
   call orthred( irphs, mmax, imax-imin+1, nopf, ll, rr, phips, phirn, pspr, aepr, prec, maxopf, ierr )

   if( ierr .ne. 0 ) return

   ! Subroutine to write out the projectors
   call write_proj( irphs, mmax, nopf, zz, ll, rr, pspr, aepr )
   nopfs(ll) = nopf

   ! Calculate and write out the 
   do ic=1,nc
     coree = ea(ic)
     call lschfb(na(ic),la(ic),ierr,coree,rr,vfull,coreuu(:,ic),coreup,zz,mmax,mch,srel)

     write(6,'(A7,3(I2,X),F20.12,X,F20.12)') 'OCEAN: ', ic, na(ic), la(ic), ea(ic),coree

     write( fnam,  '(1a8,1i3.3,2(1a1,1i2.2))' ) 'coreorbz', nint( zz ), 'n', na( ic ), 'l', la( ic )
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(1a1,1i8)' ) '#', irphs
     do ii = 1, irphs
        write ( 99, '(2(1x,1e22.15))' ) rr( ii ), coreuu( ii, ic )
     end do
     close( unit=99 )

     call getmeznl( zz, na(ic), la(ic), irphs, nopf, maxopf, rr, coreuu(:,ic), aepr, mels(:,:,ll,ic) )

     call getfgnew( zz, na(ic), la(ic), ll, irphs, nopf, rr, coreuu(:,ic), aepr )
     call diagfg( zz, na(ic), la(ic), ll, irphs, nopf, rr, coreuu(:,ic), aepr )

     call core2core( nc, ic, zz, na, la, mmax, irphs, rr, coreuu )

   enddo
 end do

 do ic = 1, nc
   write( fnam, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', nint( zz ), 'n', na(ic), 'l', la(ic)
   open(unit=99, file=fnam, form='formatted', status='unknown' )
   rewind 99
   write ( 99, '(1i5)' ) 3
   do ll=0,lmax
     do ii = 0,3
       write(99,'(4(1x,1e15.8))' ) mels( 1:nopfs(ll),ii,ll,ic)
     end do
   end do
   close( 99 )
 end do

 write(fnam, '(A,I3.3)') 'prjfilez', nint(zz )
 open(unit=99, file=fnam, form='formatted', status='unknown' )
 ! FIX nq, dq
 write( 99, '(3i5,2x,1e15.8)' ) 0, lmax, 401, 0.05_dp
 do ll=0,lmax
   write(99, '(1i5)' ) nopfs(ll)
 end do
 close( 99 )


 deallocate(aeuu,aeup,psuu,psup)
 deallocate(pshf,pshp)
 deallocate(phips,phirn,pspr,aepr,coreuu,coreup)
 deallocate(mels)
 return
 end subroutine run_opf
