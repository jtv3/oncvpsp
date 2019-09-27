!! Copyright (C) 2008-2013 M. Oliveira, F. Nogueira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: functionals.F90 817 2014-01-14 17:35:40Z micael $


module functionals_m

  use xc_f90_types_m !libxc
  use xc_f90_lib_m !libxc
  implicit none

  real(8), parameter :: pi=3.141592653589793238462643383279502884197d0

                    !---Interfaces---!

  interface assignment (=)
    module procedure functional_copy
  end interface


                    !---Derived Data Types---!

  type functional_t
    private
    integer :: family  ! LDA, GGA, etc.
    integer :: id      ! identifier
    integer :: nspin   ! XC_UNPOLARIZED | XC_POLARIZED

    type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
    type(xc_f90_pointer_t) :: info ! information about the functional

    integer  :: deriv_method
    integer  :: irel
    real(8) :: xalpha
  end type functional_t


                    !---Global Variables---!

  !
  integer, parameter :: XC_DERIV_NUMERICAL  = 1, &
                        XC_DERIV_ANALYTICAL = 2

  !Some functionals not available in libxc
  integer, parameter :: XC_OEP_X_EXX  = 600, &
                        XC_OEP_XC_SIC = 601, &
                        XC_MGGA_K_GE2 = 599

  !Some version of libxc break backward compatibility...
#if LIBXC_VERSION == 100
  integer, parameter :: XC_KINETIC = 3, &
                        XC_GGA_X_LB = XC_GGA_XC_LB
#endif


                    !---Public/Private Statements---!

  private
  public :: functional_t, &
            functional_null, &
            functional_init, &
            assignment(=), &
            functional_get_vxc, &
            functional_get_tau, &
            functional_name, &
            functional_kind, &
            functional_family, &
            functional_save, &
            functional_load, &
            functional_end,  &
            XC_DERIV_NUMERICAL, &
            XC_DERIV_ANALYTICAL, &
            XC_OEP_X_EXX, &
            XC_OEP_XC_SIC, &
            XC_MGGA_K_GE2

contains

  subroutine functional_null(functl)
    !-----------------------------------------------------------------------!
    ! Nullifies and sets to zero all the components of the functional.      !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(out) :: functl


    functl%family       = 0
    functl%id           = 0
    functl%nspin        = 0
    functl%deriv_method = 0
    functl%irel         = 0
    functl%xalpha       =0.0d0

  end subroutine functional_null

  subroutine functional_init(nspin, id, deriv_method, irel, functl)
    !-----------------------------------------------------------------------!
    ! Initializes a functional.                                             !
    !-----------------------------------------------------------------------!
    integer,            intent(in)    :: nspin
    integer,            intent(in)    :: id, deriv_method, irel
    type(functional_t), intent(inout) :: functl


    if (.not. (deriv_method == XC_DERIV_NUMERICAL .or. deriv_method == XC_DERIV_ANALYTICAL))&
&       stop 'functionals: ERROR deriv_method not well defined'

    functl%id = id
    functl%nspin = nspin
    functl%deriv_method = deriv_method
    functl%irel = irel

    if(functl%id /= 0) then
      ! get the family of the functional
      if (id == XC_OEP_X_EXX) then
        functl%family = XC_FAMILY_OEP
      elseif (id == XC_MGGA_K_GE2) then
        functl%family = XC_FAMILY_MGGA
      else
        functl%family = xc_f90_family_from_id(functl%id)
      end if

      if (functl%family == XC_FAMILY_UNKNOWN) then
        write(*, '(a,i3,a)') "'", functl%id,"' is not a known functional!"
        write(*, *) "Please check the manual for a list of possible values."
        stop 'functionals: ERROR unknown functional'
      end if

    end if

    !Extra variables
    if (functl%id == XC_LDA_C_XALPHA) then
stop 'runctionals: ERROR - parsing of xalpha coefficient not coded'
      !call oct_parse_float('Xalpha', 1.0d0, functl%xalpha)
    end if

    !Initialize
    if(functl%id /= 0 .and. functl%id /= XC_MGGA_K_GE2) then
      call functional_libxc_init(functl)

      if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        if ( ( (iand(xc_f90_info_flags(functl%info), XC_FLAGS_HAVE_FXC) == 0) .and. &
             (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA)   ) .or. &
             ( (iand(xc_f90_info_flags(functl%info), XC_FLAGS_HAVE_KXC) == 0) .and. &
             (functl%family == XC_FAMILY_MGGA)                                       )        ) then
          write(*, '(a,i3,a)') "It is not possible to compute the analytical derivatives for functional '", functl%id,"'"
          stop 'functionals ERROR line 161'
        end if
      end if

    end if

  end subroutine functional_init

  subroutine functional_libxc_init(functl)
    !-----------------------------------------------------------------------!
    ! Initialize the libxc objects of the functional.                       !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(inout) :: functl

    if (functl%family /= XC_FAMILY_OEP) then
      call xc_f90_func_init(functl%conf, functl%info, functl%id, functl%nspin)

      if (functl%id == XC_LDA_C_XALPHA) then
        call xc_f90_lda_c_xalpha_set_par(functl%conf, functl%xalpha)
      end if

      if (functl%id == XC_LDA_X) then
#if LIBXC_VERSION >= 200
        call xc_f90_lda_x_set_par(functl%conf, 4.0d0/3.0d0, functl%irel,0.0d0)
#else
        call xc_f90_lda_x_set_par(functl%conf, functl%irel)
#endif
      end if
    end if

  end subroutine functional_libxc_init

  subroutine functional_copy(functl_out, functl_in)
    !-----------------------------------------------------------------------!
    ! Copies the functional functl_in to functl_out.                        !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(inout) :: functl_out
    type(functional_t), intent(in)    :: functl_in


    call functional_end(functl_out)

    functl_out%family = functl_in%family
    functl_out%id     = functl_in%id
    functl_out%nspin  = functl_in%nspin

    functl_out%deriv_method = functl_in%deriv_method
    functl_out%irel = functl_in%irel
    functl_out%xalpha = functl_in%xalpha

    if(functl_out%id /= 0 .and. functl_out%id /= XC_MGGA_K_GE2) then
      call functional_libxc_init(functl_out)
    end if

  end subroutine functional_copy

  subroutine functional_end(functl)
    !-----------------------------------------------------------------------!
    ! Frees all memory associated to the functional.                        !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(inout) :: functl


    select case (functl%family)
    case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA)
      if (functl%id /= XC_MGGA_K_GE2) then
        call xc_f90_func_end(functl%conf)
      end if
    end select

    functl%family = 0
    functl%id     = 0
    functl%nspin  = 0

    functl%deriv_method = 0
    functl%irel         = 0
    functl%xalpha       =0.0d0

  end subroutine functional_end

  function functional_name(functl)
    !-----------------------------------------------------------------------!
    ! Returns the name of the functional.                                   !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(in) :: functl
    character(120) :: functional_name

    select case (functl%id)
    case (0)
      functional_name = "None"
    case (XC_OEP_X_EXX)
      functional_name = "Exact Exchange"
    case (XC_MGGA_K_GE2)
      functional_name = "Second-order gradient expansion of the kinetic energy density"
    case default
      call xc_f90_info_name(functl%info, functional_name)
    end select

  end function functional_name

  function functional_kind(functl)
    !-----------------------------------------------------------------------!
    ! Returns the kind of functional we have                                !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(in) :: functl
    integer :: functional_kind

    select case (functl%id)
    case (0)
      functional_kind = -1
    case (XC_OEP_X_EXX)
      functional_kind = XC_EXCHANGE
    case (XC_OEP_XC_SIC)
      functional_kind = XC_EXCHANGE_CORRELATION
    case (XC_MGGA_K_GE2)
      functional_kind = XC_KINETIC
    case default
      functional_kind = xc_f90_info_kind(functl%info)
    end select

  end function functional_kind

  elemental function functional_family(functl)
    !-----------------------------------------------------------------------!
    ! Returns the family of the functional                                  !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(in) :: functl
    integer :: functional_family

    functional_family = functl%family

  end function functional_family

  subroutine functional_save(unit, functl)
    !-----------------------------------------------------------------------!
    ! Writes the functional data to a file.                                 !
    !-----------------------------------------------------------------------!
    integer,            intent(in) :: unit
    type(functional_t), intent(in) :: functl


    write(unit) functl%family
    write(unit) functl%id
    write(unit) functl%nspin

    write(unit) functl%deriv_method
    write(unit) functl%irel
    write(unit) functl%xalpha

  end subroutine functional_save

  subroutine functional_load(unit, functl)
    !-----------------------------------------------------------------------!
    ! Reads the exchange-correlation model data from a file.                !
    !-----------------------------------------------------------------------!
    integer,            intent(in)    :: unit
    type(functional_t), intent(inout) :: functl


    read(unit) functl%family
    read(unit) functl%id
    read(unit) functl%nspin

    read(unit) functl%deriv_method
    read(unit) functl%irel
    read(unit) functl%xalpha

    if(functl%id /= 0 .and. functl%id /= XC_MGGA_K_GE2) then
      call functional_libxc_init(functl)
    end if

  end subroutine functional_load

  subroutine functional_get_vxc(functl, np, al, rr, rho, rho_grad, rho_lapl, tau, ip, v, e, vtau)
    !-----------------------------------------------------------------------!
    ! Given a density, computes the corresponding exchange/correlation      !
    ! potentials and energies.                                              !
    !                                                                       !
    !  functl   - functional                                                !
    !  np       - number of mesh points                                     !
    !  al       - grid parameter                                            !
    !  rr       - radial points                                             !
    !  rho      - electronic radial density                                 !
    !  rho_grad - gradient of the electronic radial density                 !
    !  rho_lapl - laplacian of the electronic radial density                !
    !  tau      - radial kinetic energy density                             !
    !  ip       - ionization potential                                      !
    !  v        - potential                                                 !
    !  e        - energy per-volume                                         !
    !  vtau     - extra term arising from MGGA potential                    !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(inout) :: functl
    integer     ,       intent(in)    :: np
    real(8),           intent(in)    :: al
    real(8),           intent(in)    :: rr(np)
    real(8),           intent(in)    :: rho(np, functl%nspin)
    real(8),           intent(in)    :: rho_grad(np, functl%nspin)
    real(8),           intent(in)    :: rho_lapl(np, functl%nspin)
    real(8),           intent(in)    :: tau(np, functl%nspin)
    real(8),           intent(in)    :: ip(functl%nspin)
    real(8),           intent(out)   :: v(np, functl%nspin), e(np)
    real(8),           intent(out)   :: vtau(np, functl%nspin)

    integer  :: i, is, nspin
    real(8) :: a, b, c
    real(8), parameter   :: alpha = -0.012d0, beta = 1.023d0

    ! Global variables
    real(8), allocatable :: dedrho(:,:), dedgrad(:,:), dedlapl(:,:), dedtau(:,:)
    real(8), allocatable :: d2edrhodgrad(:,:), d2edgrad2(:,:)

    ! Local variables
    real(8), allocatable :: n(:), s(:), l(:), t(:)
    real(8), allocatable :: dedn(:), deds(:), dedl(:), dedt(:)
    real(8), allocatable :: d2edn2(:), d2eds2(:), d2ednds(:)

    real(8), allocatable :: dpr(:), dppr(:), dlap(:)


    if (.not. (size(v, dim=2) == functl%nspin)) stop 'functionals: ERROR bad nspin definition'
    if (.not. (functional_kind(functl) /= XC_KINETIC)) stop 'functionals :ERROR no kinetic functional here'

    ! Initialize all output quantities to zero
    v =0.0d0 ; e =0.0d0 ; vtau =0.0d0

    ! If the functional is not set, there is nothing to be done
    if (functl%family == 0) then
      return
    end if

    ! Shortcut
    nspin = functl%nspin


    ! Compute c parameter of the TB09 functional
    if (functl%id == XC_MGGA_X_TB09) then
      if (maxval(ip) ==0.0d0) then
        c = 1.0d0
      else
        a =0.0d0
        do
          c = alpha + beta*sqrt(2.0d0*sqrt(2.0d0*(maxval(ip) + a)))
          b = (3.0d0*c - 2.0d0)/pi*sqrt(5.0d0/6.0d0*(maxval(ip) + a))
          if (abs(a - b) < 1.0e-8) exit
          a = b
        end do
      end if
      call xc_f90_mgga_x_tb09_set_par(functl%conf, c)
    end if


                    !---Allocate work arrays---!

    ! LDA
    allocate(n(nspin), dedn(nspin))
    allocate(dedrho(np, nspin))
    n =0.0d0; dedn =0.0d0
    dedrho =0.0d0

    if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
      if (nspin == 1) then
        allocate(d2edn2(1))
      else
        allocate(d2edn2(3))
      end if
      d2edn2 =0.0d0
    end if

    ! GGA
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
      if (nspin == 1) then
        allocate(s(1), deds(1))
      else
        allocate(s(3), deds(3))
      end if
      allocate(dedgrad(np, nspin))
      s =0.0d0; deds =0.0d0
      dedgrad =0.0d0

      if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        if (nspin == 1) then
          allocate(d2eds2(1), d2ednds(1))
          allocate(d2edgrad2(np, 1))
          allocate(d2edrhodgrad(np, 1))
        else
          allocate(d2eds2(6), d2ednds(6))
          allocate(d2edgrad2(np, 3))
          allocate(d2edrhodgrad(np, 4))
        end if
        d2eds2 =0.0d0; d2ednds =0.0d0
        d2edgrad2 =0.0d0; d2edrhodgrad =0.0d0
      end if

    end if

    ! MGGA
    if (functl%family == XC_FAMILY_MGGA) then
      allocate(l(nspin), dedl(nspin))
      allocate(t(nspin), dedt(nspin))
      allocate(dedlapl(np, nspin))
      allocate(dedtau(np, nspin))
      l =0.0d0; dedl =0.0d0
      t =0.0d0; dedt =0.0d0    
      dedlapl =0.0d0
      dedtau =0.0d0

      if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        !Not yet implemented
      end if

    end if


                    !---Space loop---!

    do i = 1, np
      ! make a local copy with the correct memory order
      n(1:nspin) = rho(i, 1:nspin)
      if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
        s(1) = rho_grad(i, 1)**2
        if(nspin == 2) then
          s(2) = rho_grad(i, 1)*rho_grad(i, 2)
          s(3) = rho_grad(i, 2)**2
        end if
      end if
      if (functl%family == XC_FAMILY_MGGA) then
#if LIBXC_VERSION >= 200
        t(1:nspin) = tau(i, 1:nspin)/2.0d0
#else
        t(1:nspin) = tau(i, 1:nspin)
#endif     
        l(1:nspin) = rho_lapl(i, 1:nspin)
      end if

      if (iand(xc_f90_info_flags(functl%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
 
        select case(functl%family)
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc_vxc(functl%conf, 1, n(1), e(i), dedn(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc_vxc(functl%conf, 1, n(1), s(1), e(i), dedn(1), deds(1))
        case(XC_FAMILY_MGGA)
          call xc_f90_mgga_exc_vxc(functl%conf, 1, n(1), s(1), l(1), t(1), e(i), &
                                   dedn(1), deds(1), dedl(1), dedt(1))
        end select

      else !Just get the potential

        select case(functl%family)
        case(XC_FAMILY_LDA)
          call xc_f90_lda_vxc(functl%conf, 1, n(1), dedn(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_vxc(functl%conf, 1, n(1), s(1), dedn(1), deds(1))
        case(XC_FAMILY_MGGA)
          call xc_f90_mgga_vxc(functl%conf, 1, n(1), s(1), l(1), t(1), &
                               dedn(1), deds(1), dedl(1), dedt(1))
        end select
        e(i) =0.0d0

      end if

      e(i) = e(i)*sum(n)
      dedrho(i, :) = dedn(:)
      if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
        if (nspin == 1) then
          dedgrad(i, 1) = 2.0d0*deds(1)*rho_grad(i, 1)
        else
          dedgrad(i, 1) = 2.0d0*deds(1)*rho_grad(i, 1) + deds(2)*rho_grad(i, 2)
          dedgrad(i, 2) = 2.0d0*deds(3)*rho_grad(i, 2) + deds(2)*rho_grad(i, 1)
        end if
      end if

      if(functl%family == XC_FAMILY_MGGA) then
        dedlapl(i, 1:nspin) = dedl(1:nspin)
#if LIBXC_VERSION >= 200
        dedtau(i, 1:nspin) = dedt(1:nspin)/2.0d0
#else
        dedtau(i, 1:nspin) = dedt(1:nspin)
#endif
      end if

      if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        !Evaluate second-order derivatives
        if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
          call xc_f90_gga_fxc(functl%conf, 1, n(1), s(1), d2edn2(1), d2ednds(1), d2eds2(1))

          if (nspin == 1) then
            d2edrhodgrad(i, 1) = 2.0d0*rho_grad(i, 1)*d2ednds(1)
            d2edgrad2(i, 1) = 2.0d0*deds(1) + 4.0d0*s(1)*d2eds2(1)
          else
            d2edrhodgrad(i, 1) = 2.0d0*rho_grad(i, 1)*d2ednds(1) + rho_grad(i, 2)*d2ednds(2)
            d2edrhodgrad(i, 2) = 2.0d0*rho_grad(i, 1)*d2ednds(4) + rho_grad(i, 2)*d2ednds(5)
            d2edrhodgrad(i, 3) = 2.0d0*rho_grad(i, 2)*d2ednds(3) + rho_grad(i, 1)*d2ednds(2)
            d2edrhodgrad(i, 4) = 2.0d0*rho_grad(i, 2)*d2ednds(6) + rho_grad(i, 1)*d2ednds(5)

            d2edgrad2(i, 1) = 2.0d0*deds(1) + 4.0d0*s(1)*d2eds2(1) + 4.0d0*s(2)*d2eds2(2) + s(3)*d2eds2(4)
            d2edgrad2(i, 2) = deds(2) + 4.0d0*s(2)*d2eds2(3) + 2.0d0*s(1)*d2eds2(2) + 2.0d0*s(3)*d2eds2(5) + s(2)*d2eds2(4)
            d2edgrad2(i, 3) = 2.0d0*deds(3) + 4.0d0*s(3)*d2eds2(6) + 4.0d0*s(2)*d2eds2(5) + s(1)*d2eds2(4)
          end if
          
        else if (functl%family == XC_FAMILY_MGGA) then
          !Not yet implemented
        end if
      end if

    end do ! loop over points i


                    !---Compute potentials---!

    ! LDA contribution
    v = dedrho

    ! GGA contribution
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
      if (functl%deriv_method == XC_DERIV_NUMERICAL) then
        allocate (dpr(np))
        allocate (dppr(np))
        allocate (dlap(np))
        do is = 1, nspin
          call derivs(np, dedgrad(:,is), al, rr, dpr, dppr, dlap)
!          v(:, is) = v(:, is) - mesh_divergence(m, dedgrad(:, is))
          v(:, is) = v(:, is) - (2.0d0/rr(:)*dedgrad(:,is) + dpr(:))
        end do
        deallocate (dpr, dppr, dlap)
      elseif (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        stop 'functionals: ERROR - no analytical derivatives coded'
      end if
    end if

    ! MGGA contribution
    if (functl%family == XC_FAMILY_MGGA) then
stop 'functrionals: ERROOR : meta GGA not coded yet '
      if (functl%deriv_method == XC_DERIV_NUMERICAL) then
        allocate (dpr(np))
        allocate (dppr(np))
        allocate (dlap(np))
        do is = 1, nspin
          call derivs(np, dedlapl(:,is), al, rr, dpr, dppr, dlap)
!          v(:, is) = v(:, is) + mesh_laplacian(m, dedlapl(:, is))
          v(:, is) = v(:, is) + dlap(:)
        end do
        deallocate (dpr, dppr, dlap)
      elseif (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        !Not yet implemented
      end if

      vtau = dedtau
    end if

    !Shift potentials that do not go to zero at infinity
    do is = 1, nspin
      select case (functl%id)
      case (XC_MGGA_X_BJ06)
        a = sqrt(5.0d0/6.0d0)/pi
      case (XC_MGGA_X_TB09)
        a = (3.0d0*c - 2.0d0)*sqrt(5.0d0/6.0d0)/pi
#if LIBXC_VERSION >= 210
      case (XC_GGA_X_AK13)
        a = sqrt(2.0d0)*(1.0d0/(54.0d0*pi) + 2.0d0/15.0d0)
#endif
      case default
        a =0.0d0
      end select
      do i = np, 1, -1
        if (v(i, is) /=0.0d0) then
          v(1:i, is) = v(1:i, is) - a*sqrt(ip(is))
          exit
        end if
      end do
    end do

                    !---Deallocate arrays---!

    ! LDA
    deallocate(n, dedn, dedrho)
    if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
      deallocate(d2edn2)
    end if

    ! GGA
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
      deallocate(s, deds, dedgrad)
      if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        deallocate(d2eds2, d2ednds, d2edgrad2, d2edrhodgrad)
      end if
    end if

    ! MGGA
    if (functl%family == XC_FAMILY_MGGA) then
      deallocate(l, dedl, dedlapl)
      deallocate(t, dedt, dedtau)
      if (functl%deriv_method == XC_DERIV_ANALYTICAL) then
        !Not yet implemented
      end if
    end if

  end subroutine functional_get_vxc

  subroutine functional_get_tau(functl, np, rho, rho_grad, rho_lapl, tau)
    !-----------------------------------------------------------------------!
    ! Computes the approximated kinetic energy density.                     !
    !                                                                       !
    !  functl   - functional                                                !
    !  m        - mesh                                                      !
    !  rho      - electronic radial density                                 !
    !  rho_grad - gradient of the electronic radial density                 !
    !  rho_lapl - laplacian of the electronic radial density                !
    !  tau      - radial kinetic energy density                             !
    !-----------------------------------------------------------------------!
    type(functional_t), intent(in)  :: functl
    integer     ,       intent(in)  :: np
    real(8),           intent(in)  :: rho(np, functl%nspin)
    real(8),           intent(in)  :: rho_grad(np, functl%nspin)
    real(8),           intent(in)  :: rho_lapl(np, functl%nspin)
    real(8),           intent(out) :: tau(np, functl%nspin)

    integer  :: i, is, nspin
    real(8), allocatable :: n(:), s(:), l(:), t(:)


    if (.not. (functional_kind(functl) == XC_KINETIC)) stop 'functionals: ERROR no kinetic part here'

    if (functl%id == XC_MGGA_K_GE2) then
      where (rho <= 1e-30)
        tau =0.0d0
      elsewhere        
        tau = 3.0d0/5.0d0*(3.0d0*pi**2)**(2.0d0/3.0d0)*rho**(5.0d0/3.0d0) + &
             rho_lapl/3.0d0 + rho_grad**2/rho/36.0d0
      end where
      return
    end if

    nspin = functl%nspin

    !Allocate work arrays
    allocate(n(nspin), t(nspin))
    n =0.0d0; t =0.0d0
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
      if (nspin == 1) then
        allocate(s(1))
      else
        allocate(s(3))
      end if
      s =0.0d0
    end if
    if (functl%family == XC_FAMILY_MGGA) then
      allocate(l(nspin))
      l =0.0d0
    end if

    !Spin loop
    do is = 1, nspin
      !Space loop
      do i = 1, np
        ! make a local copy with the correct memory order
        n(is) = rho(i, is)

        if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) then
          s(1) = rho_grad(i, 1)**2
          if(nspin == 2) then
            s(2) = rho_grad(i, 1)*rho_grad(i, 2)
            s(3) = rho_grad(i, 2)**2
          end if
        end if
        if (functl%family == XC_FAMILY_MGGA) then
          l(is) = rho_lapl(i, is)
        end if

        select case(functl%family)
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(functl%conf, 1, n(1), t(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(functl%conf, 1, n(1), s(1), t(1))
        end select

#if LIBXC_VERSION >= 200
        tau(i, is) = 2.0d0*t(1)*n(is)
#else
        tau(i, is) = t(1)*n(is)
#endif
      end do
    end do

    !Deallocate arrays
    deallocate(n, t)
    if (functl%family == XC_FAMILY_GGA .or. functl%family == XC_FAMILY_MGGA) deallocate(s)
    if (functl%family == XC_FAMILY_MGGA) deallocate(l)

  end subroutine functional_get_tau


end module functionals_m

