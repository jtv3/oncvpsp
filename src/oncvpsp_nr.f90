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
 program oncvpsp
!
! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
! and references therein.
!
!   D. R. Hamann
!   Mat-Sim Research LLC
!   P.O. Box 742
!   Murray Hill, NJ 07974
!   USA
!
!   Developed from original "gncpp" code of March 8,1987
!
!   Output format for ABINIT pspcod=8
!
 implicit none
 integer, parameter :: dp=kind(1.0d0)

!
 integer :: ii,ierr,iexc,iexct,ios,iprint,irps,it,icmod,lpopt
 integer :: jj,kk,ll,l1,lloc,lmax,lt,inline
 integer :: mch,mchf,mmax,n1,n2,nc,nlim,nlloc,nlmax,nrl
 integer :: nv,irct,ncnf
!
 integer :: dtime(8),na(30),la(30),np(6)
 integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 integer :: irc(6),nodes(4)
 integer :: ncut,nxtra,nproj(6),npx(6),lpx(6)
 integer :: ncon(6),nbas(6)

 real(dp) :: al,amesh,csc,csc1,deblt,depsh,depsht,drl,eeel
 real(dp) :: eeig,eexc
 real(dp) :: emax,epsh1,epsh1t,epsh2,epsh2t
 real(dp) :: et,etest,emin,sls
 real(dp) :: fcfact,rcfact,dvloc0
 real(dp) :: dvloc0t
 real(dp) :: rr1,rcion,rciont,rcmax,rct,rlmax,rpkt
 real(dp) :: sf,zz,zion,zval,etot
 real(dp) :: xdummy
!
 real(dp) :: cl(6),evkb(2,4),debl(6),ea(30),eacopy(30),ep(6),fa(30),facnf(30,5)
 real(dp) :: fnp(6),fp(6)
 real(dp) :: qcut(6),rc(6),rc0(6)
 real(dp) :: rpk(30),cvgplt(2,7,2,4)
 real(dp) :: epx(6),fpx(6)
 real(dp) :: qq(2,2)
 real(dp) :: epstot
 real(dp), parameter :: eps=1.0d-8

 real(dp), allocatable :: rr(:)
 real(dp), allocatable :: rho(:),rhoc(:),rhot(:)
 real(dp), allocatable :: uu(:),uu2(:),up(:),up2(:)
 real(dp), allocatable :: vp(:,:),vfull(:),vkb(:,:,:),pswf(:,:,:)
 real(dp), allocatable :: vkbt(:,:,:)
 real(dp), allocatable :: vwell(:)
 real(dp), allocatable :: vpuns(:,:)
 real(dp), allocatable :: vo(:),vxc(:)
 real(dp), allocatable :: rhomod(:,:),rhoae(:,:),rhops(:,:),rhotae(:)
 real(dp), allocatable :: uua(:,:), uuatmp(:,:) !pseudo-atomic orbitals array

 character*2 :: atsym
 character*4 :: psfile

 logical :: srel

 write(6,'(a/a//)') &
&      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&      'scalar-relativistic version 3.3.1 12/12/2017'

 write(6,'(a/a/a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&      'in any publication utilizing these pseudopotentials.'

 srel=.false.

 nproj(:)=0
 fcfact=0.0d0
 rcfact=0.0d0
 rc(:)=0.0d0
 ep(:)=0.0d0

! read input data
 inline=0

! atom and reference configuration
 call cmtskp(inline)
 read(5,*,iostat=ios) atsym,zz,nc,nv,iexc,psfile
 call read_error(ios,inline)

 call cmtskp(inline)
 do ii=1,nc+nv
   read(5,*,iostat=ios) na(ii),la(ii),fa(ii)
   call read_error(ios,inline)
 end do

! pseudopotential and optimization
 call cmtskp(inline)
 read(5,*,iostat=ios) lmax
 call read_error(ios,inline)

 call cmtskp(inline)
 do l1=1,lmax+1
   read(5,*,iostat=ios) lt,rc(l1),ep(l1),ncon(l1),nbas(l1),qcut(l1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
 end do

! local potential
 call cmtskp(inline)
 read(5,*, iostat=ios) lloc,lpopt,rc(5),dvloc0
 call read_error(ios,inline)

! Vanderbilt-Kleinman-Bylander projectors
 call cmtskp(inline)
 do l1=1,lmax+1
   read(5,*,iostat=ios) lt,nproj(l1),debl(l1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
 end do

! model core charge
 call cmtskp(inline)
 read(5,*, iostat=ios) icmod, fcfact
 if(ios==0 .and. icmod==3) then
  backspace(5)
  read(5,*, iostat=ios) icmod, fcfact, rcfact
 end if
 call read_error(ios,inline)

! log derivative analysis
 call cmtskp(inline)
 read(5,*,iostat=ios) epsh1, epsh2, depsh
 call read_error(ios,inline)

! output grid
 call cmtskp(inline)
 read(5,*,iostat=ios) rlmax,drl
 call read_error(ios,inline)

! test configurations
 call cmtskp(inline)
 read(5,*,iostat=ios) ncnf
 call read_error(ios,inline)

 do jj=2,ncnf+1
   call cmtskp(inline)
   read(5,*,iostat=ios) nvcnf(jj)
   call read_error(ios,inline)
   do ii=nc+1,nc+nvcnf(jj)
     call cmtskp(inline)
     read(5,*,iostat=ios) nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
     call read_error(ios,inline)
   end do
 end do

! end of reading input data

 nvcnf(1)=nv
 do ii=1,nc+nv
   nacnf(ii,1)=na(ii)
   lacnf(ii,1)=la(ii)
   facnf(ii,1)=fa(ii)
 end do

 do jj=2,ncnf+1
   do ii=1,nc
     nacnf(ii,jj)=na(ii)
     lacnf(ii,jj)=la(ii)
     facnf(ii,jj)=fa(ii)
   end do
 end do

 if(icmod==0) then
   fcfact=0.0d0
 end if

 call check_data(atsym,zz,fcfact,rcfact,epsh1,epsh2,depsh,rlmax,drl,fa,facnf, &
&                rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile)

 nrl=int((rlmax/drl)-0.5d0)+1

!PWSCF wants an even number of mesh pointe
!if(trim(psfile)=='upf') then
  if(mod(nrl,2)/=0) nrl=nrl+1
!end if

 amesh=1.012d0
 al=dlog(amesh)
 rr1=.0005d0/zz
 mmax=dlog(45.0d0 /rr1)/al

!calculate zion for output
 zion=zz
 do ii=1,nc
  zion=zion-fa(ii)
 end do

 allocate(rr(mmax))
 allocate(rho(mmax),rhoc(mmax),rhot(mmax))
 allocate(uu(mmax),uu2(mmax),up(mmax),up2(mmax))
 allocate(vp(mmax,5),vfull(mmax),vkb(mmax,2,4),pswf(mmax,2,4))
 allocate(vkbt(mmax,2,4))
 allocate(vwell(mmax))
 allocate(vpuns(mmax,5))
 allocate(vo(mmax),vxc(mmax))
 allocate(rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax))

 vp(:,:)=0.0d0
 vkb(:,:,:)=0.0d0
 do ii=1,mmax
  rr(ii)=rr1*exp(al*(ii-1))
 end do

!
! full potential atom solution
!
   call sratom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
&              rr,vfull,zz,mmax,iexc,etot,ierr,srel)
!
!


! collect pseudopotential information from all-electron configuration
! and eigenvalues
! nxtra counts valence states where semi-core of same l have
! been included as "valence"
 

 np(:)=0
 fp(:)=0.0d0

 nxtra=0
 npx(:)=0
 lpx(:)=0
 fpx(:)=0.0d0
 epx(:)=0.0d0
   
 
 if(nv/=0) then
   emax=0.0d0
   do ii=nc+1,nc+nv
     l1=la(ii)+1
     emax=dmax1(emax,ea(ii))
     if(np(l1)==0) then
       np(l1)=na(ii)
       fp(l1)=fa(ii)
       ep(l1)=ea(ii)
     else if (fa(ii) .gt. 0.0d0) then
       nxtra = nxtra + 1
       npx(nxtra)=na(ii)
       lpx(nxtra)=la(ii)
       fpx(nxtra)=fa(ii)
       epx(nxtra)=ea(ii)
       debl(l1)=ea(ii)-ep(l1)
     end if
   end do
 end if

! Drop digits beyond 5 decimals for input rcs before making any use of them
 do l1=1,max(lmax+1,lloc+1)
    jj=int(rc(l1)*10.0d5)
    rc(l1)=jj/10.0d5
 end do

 rcmax=0.0d0
 do l1=1,lmax+1
   rcmax=dmax1(rcmax,rc(l1))
 end do
 do l1=1,lmax+1
   if(rc(l1)==0.0d0) then
     rc(l1)=rcmax
   end if
 end do

! principal quantum number for all-electron well states
 do ii=nc,1,-1
   l1=la(ii)+1
   if(np(l1)==0) then
    np(l1)=na(ii)+1
   end if
 end do

 do l1=1,lmax+1
  if(np(l1)==0) then
   np(l1)=l1
  end if
 end do

 nproj(lloc+1)=0
 rc0(:)=rc(:)

! output printing (echos input data, with all-electron eigenvalues added)

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z   nc   nv     iexc    psfile'
 write(6,'(a,a,f6.2,2i5,i8,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',psfile
 write(6,'(a/a)') '#','#   n    l    f        energy (Ha)'
 do ii=1,nc+nv
   write(6,'(2i5,f8.2,1pd18.7)') na(ii),la(ii),fa(ii),ea(ii)
 end do

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 write(6,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2f10.5,2i5,f10.5)') l1-1,rc(l1),ep(l1),ncon(l1),&
&        nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
 end do

 write(6,'(a/a/a)') '#','# MODEL CORE CHARGE', &
&      '# icmod, fcfact, rcfact'
 write(6,'(i5,2f10.5)') icmod,fcfact,rcfact

 write(6,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
&      '# epsh1, epsh2, depsh'
 write(6,'(3f8.2)') epsh1,epsh2,depsh

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.2)') rlmax,drl

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do

 write(6,'(//a)') 'Reference configufation results'
 write(6,'(a,i6)') '  iterations',it
 if(it .ge. 100) then
  write(6,'(a)') 'oncvpsp: ERROR all-electron reference atom not converged'
  stop
 end if
 write(6,'(a,1p,d18.8)') '  all-electron total energy (Ha)',etot
  
!find log mesh point nearest input rc
 rcmax=0.0d0
 irc(:)=0
 do l1=1,max(lmax+1,lloc+1)
  if(rc(l1)==0.0d0) cycle
  rct=rc(l1)
  do ii=2,mmax
   if(rr(ii)>rct) then
    irc(l1)=ii
    rc(l1)=rr(ii)
    exit
   end if
  end do
  rcmax=dmax1(rcmax,rc(l1))
 end do

!
 cvgplt(:,:,:,:)=0.0d0
!
! loop to construct pseudopotentials for all angular momenta
!
 write(6,'(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
&      'and semi-local pseudopoentials for all angular momenta'

 do l1=1,lmax+1
   ll=l1-1
   uu(:)=0.0d0; uu2(:)=0.0d0; qq(:,:)=0.0d0

! all-electron wave function for first projector
   if(fp(l1) .gt. 0.0d0) then
     call lschfb(np(l1),ll,ierr,ep(l1), &
&                rr,vfull,uu,up,zz,mmax,mch,srel)
     if(ierr .ne. 0) then
       write(6,'(a,i2,a,i2)') 'oncvpsp: lschfb convergence ERROR, n=',&
&       np(l1),'  l=',ll
       stop
     end if
   else
     call wellstate(np(l1),ll,irc(l1),ep(l1),rr,vfull,uu,up,zz,mmax,mch,srel)
   end if

! find outermost peak of wavefunction
   do jj=mch-1,1,-1
     if(up(jj)*up(jj+1)<0.0d0) then
       rpkt=rr(jj)
       exit
     end if
   end do

   call fpovlp(uu,uu,irc(l1),ll,zz,qq(1,1),rr,srel)

! all-electron wave function for second projector
   if(nproj(l1)==2) then
     et=ep(l1)+debl(l1)
     if(et<0.0d0) then
       call lschfb(np(l1)+1,ll,ierr,et,rr,vfull,uu2,up2,zz,mmax,mch,srel)
       if(ierr .ne. 0) then
         write(6,'(a,i2,a,i2)') 'oncvpsp: lschfb convergence ERROR, n=',&
&         np(l1),'  l=',ll 
         stop
       end if
       if(abs(et-ep(l1)-debl(l1))>eps) then
         write(6,'(/a,f10.6,a/a,f10.6,a/a)') &
&              'WARNING: negative energy of ',ep(l1)+debl(l1), &
&              ' specified for second projector', &
&              ' Actual bound-state energy of ',et, ' will be used', &
&              ' Action: debl should probably be increased'
         debl(l1)=et-ep(l1)
       end if
     else
       call wellstate(np(l1)+1,ll,irc(l1),et,rr,vfull,uu2,up2,zz,mmax,mch,srel)
     end if
     call fpovlp(uu,uu2,irc(l1),ll,zz,qq(1,2),rr,srel)
     call fpovlp(uu2,uu2,irc(l1),ll,zz,qq(2,2),rr,srel)
   end if

   write(6,'(/a,f6.3,a,i2)') &
&        'First projector wave function outermost peak radius = ',&
&        rpkt,', l=',ll

   write(6,'(/a,i2/a,f12.8,a,f12.8,a,f12.8)') &
&        'All-electron projector functions [0,rc] norms and overlap, l=',&
&        ll,'  qq(1,1)=',qq(1,1),'  qq(1,2)=',qq(1,2),'  qq(2,2)=',qq(2,2)  

! calculate optimized pseudo wave functions and related potentials
   call run_optimize(ep(l1),ep(l1)+debl(l1),ll,mmax,rr,uu,uu2,qq, &
&                    irc(l1),qcut(l1),ncon(l1),nbas(l1),nproj(l1), &
&                    pswf(1,1,l1),vp(1,l1),vkb(1,1,l1),vfull,cvgplt(1,1,1,l1))
 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(6,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

! Uncomment the following lines and line 132 of run_vkb.f90 to write fort.20
! reporting the B-matrix eigenvalues for all L as dvloc0 is changed to
! eliminate a ghost.  A range of +/- 1 Ha around the input dvloc0 is
! covered and labeled 0L, 1L, 2L, etc. for easy grepping and plotting.
! do kk=0,40
!   dvloc0t=dvloc0+0.05d0*(kk-20)
!   vkbt(:,:,:)=vkb(:,:,:)
!   call run_vkb(lmax,lloc,lpopt,dvloc0t,irc,nproj,rr,mmax,pswf,vfull,vp, &
!&               evkb,vkbt,nlim)
! end do
! stop
 call run_vkb(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,pswf,vfull,vp, &
&             evkb,vkb,nlim)


! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening

! get all-electron and pseudopotential valence-state by valence-state
! charge densities
! null charge and eigenvalue accumulators
 eeig=0.0d0
 zval=0.0d0
 rho(:)=0.0d0
 nodes(:)=0
 rhotae(:)=0.0d0
 irps=0
 do kk=1,nv

   et=ea(nc+kk)
   ll=la(nc+kk)
   l1=ll+1
   call lschfb(na(nc+kk),ll,ierr,et, &
&                rr,vfull,uu,up,zz,mmax,mch,srel)

   rhoae(:,kk)=(uu(:)/rr(:))**2

   rhotae(:)=rhotae(:) + fa(nc+kk)*rhoae(:,kk)

   sls=(l1-1)*l1
   emax=vp(mmax,l1)+0.5d0*sls/rr(mmax)**2
   emin=emax
   do jj=1,mmax
     emin=dmin1(emin,vp(jj,l1)+0.5d0*sls/rr(jj)**2)
   end do

   call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                uu,up,mmax,mch)

   if(ierr/=0) then 
     write(6,'(a,3i4)') 'oncvpsp: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
     flush(6)
     stop
   end if
   rhops(:,kk)=(uu(:)/rr(:))**2
   rho(:)=rho(:)+fa(nc+kk)*rhops(:,kk)
   eeig=eeig+fa(nc+kk)*et

   zval=zval+fa(nc+kk)
   nodes(l1)=nodes(l1)+1
   irps=max(irps,irc(l1))
 end do !kk

 allocate(rhomod(mmax,5))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
! or Teter function fit

 if(icmod==1) then
   call modcore(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod==2) then
   call modcore2(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod>=3) then
   call modcore3(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

!   call modcore3(icmod,rho,rhoc,rhomod,fcfact,rcfact,mmax,rr, &
!&                na,la,ea,fa,nc,nv,iexc,zz,vfull,srel, &
!&                irc,lloc,nproj,vp,vkb,2,evkb)
 end if

! screening potential for pseudocharge

 call vout(1,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,iexc)

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(6,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 call run_diag(lmax,np,fp,ep,debl,lloc,irc,npx,lpx,fpx,epx,nxtra, &
&              vkb,evkb,nproj,rr,vfull,vp,zz,mmax,srel)

 call run_ghosts(lmax,la,ea,nc,nv,lloc,irc,qcut, &
&                vkb,evkb,nproj,rr,vp,mmax)

! unscreen semi-local potentials

 do l1=1,max(lmax+1,lloc+1)
   vpuns(:,l1)=vp(:,l1)-vo(:)
 end do

! fix unscreening error due to greater range of all-electron charge
 do ii=mmax,1,-1
   if(rho(ii)==0.0d0) then 
     do l1=1,max(lmax+1,lloc+1)
       vpuns(ii,l1)=-zion/rr(ii)
     end do
   else
     exit
   end if
 end do

! allocate atomic orbitals
allocate(uua(mmax,nv))

! loop over reference plus test atom configurations
 rhot(:)=rho(:)/zval
 do jj=1,ncnf+1
   do ii=nv+1,nvcnf(jj)
     ea(ii)=ea(nv)
   end do
   allocate(uuatmp(mmax,nvcnf(jj)))
 
   write(6,'(/a,i2)') 'Test configuration',jj-1
   call run_config(nacnf(1,jj),lacnf(1,jj),ea,facnf(1,jj),nc,nvcnf(jj),rhot, &
&                  rhomod,rr,zz,rcmax,mmax,iexc,etot,epstot,nproj,vpuns, &
&                  lloc,vkb,evkb,srel, uuatmp)

   if(jj==1) then !only for reference config
      uua = uuatmp !save atomic orbitals
      eacopy = ea  !save eigenvalues
   end if

   deallocate(uuatmp)
 end do

 call run_plot(lmax,np,fp,ep,lloc,rc,npx,lpx,fpx,epx,nxtra, &
&              vkb,evkb,nproj,rr,vfull,vp,vpuns,zz,mmax,drl,nrl, &
&              rho,rhoc,rhomod,pswf,srel,cvgplt)

 call run_phsft(lmax,lloc,nproj,ep,epsh1,epsh2,depsh,vkb,evkb, &
&               rr,vfull,vp,zz,mmax,irc,srel)
 
 call gnu_script(lmax,lloc,nproj,nxtra,lpx)

 if(trim(psfile)=='psp8' .or. trim(psfile)=='both') then
!  call linout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
!&             zz,zion,mmax,iexc,fcfact,nrl,drl,atsym)
  call linout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&             rhotae,rhoc,zz,zion,mmax,iexc,icmod,nrl,drl,atsym, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&             epsh1,epsh2,depsh,rlmax,psfile)
 end if

 if(trim(psfile)=='upf' .or. trim(psfile)=='both') then
  call upfout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&             zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&             epsh1,epsh2,depsh,rlmax,psfile, uua,eacopy)
 end if

 stop
 end program oncvpsp


 subroutine cmtskp(inline)
! skips lines of standard input (file 5) whose first character is #
 implicit none

!In/Out variable
 integer :: inline

!Local variable
 character*1 tst

 tst='#'
 do while (tst=='#')
  read(5,*) tst
  inline=inline+1
 end do
 backspace(5)
 inline=inline-1

 return
 end subroutine cmtskp

 subroutine read_error(ios,inline)
! report data read error and stop
 implicit none

!Input variables
 integer :: ios,inline

 inline=inline+1
 if(ios/=0) then
  write(6,'(a,i4)') 'read_error: ERROR data file line',inline
  write(6,'(a)') 'Program will stop'
  stop
 end if

 return
 end subroutine read_error
