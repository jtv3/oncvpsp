!
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
 subroutine run_vkb(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,pswf,vfull,vp, &
&                   evkb,vkb,nlim)

! computes Vanderbilt / Kleinman-Bylander non-local potentials

!lmax  maximum angular momentum
!lloc  l for local potential (lloc==4 => use linear combination)
!lpopt  choice of polynomial for lloc==4
!dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
!irc  core radii indices
!nproj  number of projectors for each l
!rr  log radial grid
!mmax  size of radial grid
!pswf pseudo wave functions
!vfull  all-electron potential
!vp  semi-local pseudopotentials for first projectors 
!     (vp(:,5) is local potential if lloc=4)
!vkb  semi-local "potentials"*pswf (input); VKB projectors (output)
!evkb  coefficients of BKB projectors (output)
!nlim  index of maximum rc including that of vlocal (output)

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,lpopt,mmax
 integer :: irc(6),nproj(6)
 real(dp) :: rr(mmax),pswf(mmax,2,5),vfull(mmax),vp(mmax,5)
 real(dp) :: dvloc0

!Input/Output variables
 real(dp) :: vkb(mmax,2,4),evkb(2,4)
 integer :: nlim

!Local variables
 integer :: ii,ipk,jj,kk,ierr,l1,info
 real(dp) :: apk,sn,tt,qq12
 real(dp) :: bb(2,2),bbev(2),bbi(2,2)
 real(dp) :: sovl(2,2),sovlev(2),smhalf(2,2),sphalf(2,2)
 real(dp) :: bbist(2,2),bbistev(2),bbit(2,2)
 real(dp), allocatable :: vloc(:),vkbt(:,:),vkbst(:,:),work(:)


 allocate(vloc(mmax),vkbt(mmax,2),vkbst(mmax,2),work(10))


! calculate local potential

 if(lloc==4) then
   call vploc(rr,vfull,vp,dvloc0,irc(5),mmax,lpopt)
 end if

 vloc(:)=vp(:,lloc+1)

! Vanderbilt / Kleinman-Bylander projector construction
 nlim=irc(lloc+1)
 do l1=1,lmax+1
   nlim=max(nlim,irc(l1))
   if(l1==lloc+1) cycle
   do jj=1,nproj(l1)
     do ii=1,irc(l1)
       vkb(ii,jj,l1)=vkb(ii,jj,l1)-vloc(ii)*pswf(ii,jj,l1)
     end do
     do ii=irc(l1)+1,mmax
       vkb(ii,jj,l1)=0.0d0
     end do
   end do
 end do

! Vanderbilt B-matrix construction

 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   bb(:,:)=0.0d0
   call vpinteg(pswf(1,1,l1),vkb(1,1,l1),irc(l1),2*l1,bb(1,1),rr)
   evkb(1,l1)=1.0d0/bb(1,1)
! this is Kleinman-Bylander result if nproj==1

   if(nproj(l1)==2) then

     call vpinteg(pswf(1,1,l1),vkb(1,2,l1),irc(l1),2*l1,bb(1,2),rr)
     call vpinteg(pswf(1,2,l1),vkb(1,1,l1),irc(l1),2*l1,bb(2,1),rr)
     call vpinteg(pswf(1,2,l1),vkb(1,2,l1),irc(l1),2*l1,bb(2,2),rr)

     write(6,'(/a,i4/a,1p,d14.4)') &
&          'B matrix, ll=',l1-1,&
&          'Hermiticity error',bb(1,2)-bb(2,1)

! symmetrize exactly
     tt=0.5d0*(bb(1,2)+bb(2,1))
     bb(1,2)=tt
     bb(2,1)=tt

! find eigenvalues and eigenvectors of the B matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

     call dsyev( 'V', 'U', 2, bb, 2, bbev, work, 10, info )
     if(info .ne. 0) then
      write(6,'(a,i4)') 'run_vkb: ERROR B matrix eigenvalue error, info=',info
      stop
     end if

! take linear combinations to form diagonal projectors

     vkbt(:,:)=vkb(:,:,l1)

     evkb(1,l1)=1.0d0/bbev(1)
     evkb(2,l1)=1.0d0/bbev(2)
     vkb(:,1,l1)=bb(1,1)*vkbt(:,1) + bb(2,1)*vkbt(:,2)
     vkb(:,2,l1)=bb(1,2)*vkbt(:,1) + bb(2,2)*vkbt(:,2)
   
! Uncommennt the following line to tune away from thosts (see doc/ghosts.txt)
!    write(20,'(i4,a,3f16.6)') l1-1,'L',dvloc0,bbev(1),bbev(2)
   end if !nproj==2

! normalize projectors
   do jj=1,nproj(l1)
     call vpinteg(vkb(1,jj,l1),vkb(1,jj,l1),irc(l1),2*l1,sn,rr)
     if(vkb(irc(l1)-10,jj,l1)>=0.0d0) then
       vkb(:,jj,l1)=vkb(:,jj,l1)/sqrt(sn)
     else
       vkb(:,jj,l1)=-vkb(:,jj,l1)/sqrt(sn)
     end if
     evkb(jj,l1)=sn*evkb(jj,l1)
   end do


! calculate overlap
   if(nproj(l1)==2) then
    call vpinteg(vkb(1,1,l1),vkb(1,2,l1),irc(l1),2*l1,sn,rr)

    sovl(1,1)=1.0d0
    sovl(2,2)=1.0d0
    sovl(1,2)=sn
    sovl(2,1)=sn

    bbi(:,:)=0.0d0
    bbi(1,1)=evkb(1,l1)
    bbi(2,2)=evkb(2,l1)
   

!   write(6,'(a)') 'S matrix and eigenvalues'
!   do jj=1,2
!    write(6,'(1p,2e12.4)') (sovl(jj,ii),ii=1,2)
!   end do

! find eigenvalues and eigenvectors of the overlap matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

    call dsyev( 'V', 'U', 2, sovl, 2, sovlev, work, 10, info )
    if(info .ne. 0) then
     write(6,'(a,i4)') 'run_vkb: sovl matrix eigenvalue ERROR, info=',info
     stop
    end if

!   write(6,'(/1p,2e12.4)') (sovlev(ii),ii=1,2)

! construct S^(-1/2) AND s^(1/2)

     do jj=1,2
      tt=sqrt(sovlev(jj))
      do ii=1,2
       sphalf(ii,jj)=tt*sovl(ii,jj)
       smhalf(ii,jj)=sovl(ii,jj)/tt
      end do
     end do

! construct B^(-1)* = S^(1/2)^T B^(-1) S^(1/2)

     bbit(:,:)=0.0d0
     bbist(:,:)=0.0d0

     do ii=1,2
      do jj=1,2
       do kk=1,2
        bbit(ii,jj)=bbit(ii,jj)+bbi(ii,kk)*sphalf(kk,jj)
       end do
      end do
     end do

     do ii=1,2
      do jj=1,2
       do kk=1,2
        bbist(ii,jj)=bbist(ii,jj)+bbit(kk,jj)*sphalf(kk,ii)
       end do
      end do
     end do

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

     call dsyev( 'V', 'U', 2, bbist, 2, bbistev, work, 10, info )
     if(info .ne. 0) then
      write(6,'(a,i4)') 'run_vkb: B^(-1)* matrix eigenvalue ERROR, info=',info
      stop
     end if


! take linear combinations to form orthonormal basis functions

     vkbt(:,:)=vkb(:,:,l1)
     vkbst(:,:)=0.0d0

     do jj=1,2
      do ii=1,2
       vkbst(:,jj)=vkbst(:,jj) + smhalf(ii,jj)*vkbt(:,ii)
      end do
     end do


! take linear combinations to form orthonormal projectors

     vkb(:,:,L1)=0.0d0

     do ii=1,2
      evkb(ii,l1)=bbistev(ii)
      do jj=1,2
       vkb(:,ii,l1)=vkb(:,ii,l1) + bbist(jj,ii)*vkbst(:,jj)
      end do
     end do

! re-order if necessary so first projector has largest magnitude coefficient
! which is consistent with relativistic sr_so_r.f90 output

     if(abs(evkb(2,l1))>abs(evkb(1,l1))) then
       tt=evkb(1,l1)
       vkbt(:,1)=vkb(:,1,l1)
       evkb(1,l1)=evkb(2,l1)
       vkb(:,1,l1)=vkb(:,2,l1)
       evkb(2,l1)=tt
       vkb(:,2,l1)=vkbt(:,1)
     end if
  
     write(6,'(/a,1p,2e12.4)') '  Orthonormal projector coefficients',&
&        (evkb(jj,l1),jj=1,2)

! Set sign of projectors (physically irrelevant) so that they are positive
! at their peak (needed for compaisons apparently)

     do jj=1,nproj(l1)
       apk=0.0d0
       do ii=1,mmax
         if(abs(vkb(ii,jj,l1))>apk) then
           apk=abs(vkb(ii,jj,l1))
           ipk=ii
         end if
       end do
       if(vkb(ipk,jj,l1)<0.0d0) then
         vkb(:,jj,l1)=-vkb(:,jj,l1)
       end if
     end do

   end if !nproj(l1)==2

 end do  ! l1

 deallocate(vloc,vkbt,vkbst,work)

 return
 end subroutine run_vkb
