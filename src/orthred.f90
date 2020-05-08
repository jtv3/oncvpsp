! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine orthred( nr, mmax, ntot, nnew, ll, rr, pheps, pheae, pspr, aepr, prec, maxopf, ierr )
!
! creates the OPFs from a set of AE&PS partial waves
!
! nr 
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  
  ! Input variables
  integer,intent(in) :: nr, mmax, ntot, ll
  real(dp),intent(in) :: rr( mmax ), pheps( nr, ntot ), pheae( nr, ntot )
  real(dp),intent(in) :: prec
  integer,intent(in) :: maxopf

  ! Output variables
  real(dp),intent(out) :: pspr( nr, ntot ), aepr( nr, ntot )
  integer,intent(out) :: nnew
  integer,intent(out) :: ierr
  !
  integer :: i, j, lwork, liwork, navail, info
  real( dp ) :: tmp, tmp1, tmp2, err1, err2, err3, accuracy, dumf, su
  real( dp ) :: w( ntot ), fv1( ntot ), fv2( ntot ), fm1( 2 * ntot )
  real( dp ), allocatable, dimension( :, : ) :: ar, ai, zr, zi, eigenvectors
  real( dp ), allocatable, dimension( : ) :: eigenvalues, work
  integer, allocatable :: isuppz( : ), iwork( : )

  real(dp), external :: dlamch
  !
  nnew = 0
  ierr = 0
!  allocate( ar( ntot, ntot ), ai( ntot, ntot ), zr( ntot, ntot ), zi( ntot, ntot ) )
  allocate( ar( ntot, ntot ) )
  allocate( eigenvectors( ntot, ntot ), eigenvalues( ntot ), isuppz( 2 * ntot ) )
  !
!  ai = 0
  err1 = 0.0d0; err2 = 0.0d0
  do i = 1, ntot
     do j = 1, ntot
!        call radint( mmax, r, dl, pheps( :, i ), pheps( :, j ), tmp )
        call vpinteg(pheps( :, i ), pheps( :, j ), nr, ll, tmp, rr )
        if ( i .eq. j ) then 
           err1 = max( err1, abs( tmp - 1.0d0 ) )
        else
           err2 = max( err2, abs( tmp ) )
        end if
        ar( i, j ) = -tmp
     end do
  end do
  write ( 6, '(1a17,2(1x,1e15.8))' ) 'OCEAN: norm err =', err1, err2


!  call elsch( ntot, ntot, ar, ai, w, 1, zr, zi, fv1, fv2, fm1, ierr )
  ! Only need a few of the eigenvectors since they are ordered
  accuracy = dlamch( 'S' )
  ! First determine work sizes
  allocate( work( 1 ), iwork( 1 ) )
  lwork = -1
  liwork = -1
  call DSYEVR( 'V', 'I', 'U', ntot, ar, ntot, dumf, dumf, 1, maxopf, accuracy, navail, &
               eigenvalues, eigenvectors, ntot, isuppz, work, lwork, iwork, liwork, info )
  lwork = min( nint(work( 1 )), 26 * ntot )
  liwork = min( iwork( 1 ), 10 * ntot )
  deallocate( work, iwork )
  allocate( work( lwork ), iwork( liwork ) )
  call DSYEVR( 'V', 'I', 'U', ntot, ar, ntot, dumf, dumf, 1, maxopf, accuracy, navail, &
               eigenvalues, eigenvectors, ntot, isuppz, work, lwork, iwork, liwork, info )
  if( info .ne. 0 ) then
    write(6,*) 'OCEAN: WARNING!!!! DSYEVR failed! ', info
    ierr = info
    return
  endif
  !
  nnew = 0
  su = 0.0_dp
  do i = 1, navail
    if( (1.0_dp - su) .lt. prec ) exit 
!     if ( abs( eigenvalues( i ) ) .gt. prec ) nnew = nnew + 1
    nnew = nnew + 1
    su = su + abs(eigenvalues( i )/real(ntot,dp) ) 
  end do
  write ( 6, '(A7,6(1x,1e15.8))' ) 'OCEAN: ', eigenvalues( 1 : nnew )
  write ( 6, '(A7,2x,1a9,3(X,f10.4))' ) 'OCEAN: ', 'runsu = ', &
                sum( eigenvalues( 1 : navail ) ), sum( eigenvalues( 1 : nnew )), su
  !
!  err1 = 0.0d0; err2 = 0.d0; err3 = 0.0d0
!  do i = 1, ntot
!     err1 = max( err1, abs( sum( zr( :, i ) ** 2 ) - 1.0d0 ) )
!     err2 = max( err2, sum( zi( :, i ) ** 2 ) )
!     do j = 1, ntot
!        if ( j .ne. i ) err3 = max( err3, dot_product( zr( :, i ), zr( :, j ) ) )
!     end do
!  end do
!  write ( 6, '(1a7,5x,3(1x,1e15.8))' ) 'errs = ', err1, err2, err3
  !
  do i = 1, nnew
     eigenvectors( :, i ) = eigenvectors( :, i ) / sqrt( abs( eigenvalues( i ) ) )
  end do
  pspr( :, : ) = 0.d0; aepr( :, : ) = 0.0d0
  do i = 1, nnew
     do j = 1, ntot
        pspr( 1 : nr, i ) = pspr( 1 : nr, i ) + eigenvectors( j, i ) * pheps( 1 : nr, j )
        aepr( 1 : nr, i ) = aepr( 1 : nr, i ) + eigenvectors( j, i ) * pheae( 1 : nr, j )
     end do
  end do

  ! Check orthonormality of projectors
  err1 = 0.0d0
  err2 = 0.0d0
  do i = 1, nnew
     do j = 1, nnew
        call vpinteg( pspr( :, i ), pspr( :, j ), nr, ll, tmp1, rr )
        call vpinteg( aepr( :, i ), aepr( :, j ), nr, ll, tmp2, rr )
!        call radint( mmax, r, dl, pspr( :, i ), pspr( :, j ), tmp1 )
!        call radint( mmax, r, dl, aepr( :, i ), aepr( :, j ), tmp2 )
        tmp1 = abs( tmp1 ); tmp2 = abs( tmp2 ); 
        if ( j .eq. i ) then
           tmp1 = abs( tmp1 - 1.0d0 ); tmp2 = abs( tmp2 - 1.0d0 ); 
        end if
        err1 = max( err1, tmp1 ); err2 = max( err2, tmp2 )
     end do
  end do
  write ( 6, '(1a28,2(1x,1e15.8))' ) 'OCEAN:  new overlap error = ', err1, err2
  !
  deallocate( eigenvectors, eigenvalues, isuppz, work, iwork, ar )
  return
end subroutine orthred
