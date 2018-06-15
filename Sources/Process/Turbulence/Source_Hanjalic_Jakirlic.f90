!==============================================================================!
  subroutine Source_Hanjalic_Jakirlic(grid, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for Hanjalic-Jakirlic model.                                   !
!   Following paper "A new approach to modelling near-wall turbulence energy   !
!   and stress dissipation"                                                    !
!   DOI: 10.1017/S0022112002007905                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Work_Mod, only: l_sc_x => r_cell_01,  &
                      l_sc_y => r_cell_02,  &
                      l_sc_z => r_cell_03,  &
                      kin_x  => r_cell_04,  &
                      kin_y  => r_cell_05,  &
                      kin_z  => r_cell_06,  &
                      kin_xx => r_cell_07,  &
                      kin_yy => r_cell_08,  &
                      kin_zz => r_cell_09,  &
                      u_xx   => r_cell_10,  &
                      u_yy   => r_cell_11,  &
                      u_zz   => r_cell_12,  &
                      u_xy   => r_cell_13,  &
                      u_xz   => r_cell_14,  &
                      u_yz   => r_cell_15,  &
                      v_xx   => r_cell_16,  &
                      v_yy   => r_cell_17,  &
                      v_zz   => r_cell_18,  &
                      v_xy   => r_cell_19,  &
                      v_xz   => r_cell_20,  &
                      v_yz   => r_cell_21,  &
                      w_xx   => r_cell_22,  &
                      w_yy   => r_cell_23,  &
                      w_zz   => r_cell_24,  &
                      w_xy   => r_cell_25,  &
                      w_xz   => r_cell_26,  &
                      w_yz   => r_cell_27,  &
                      kin_sq => r_cell_28,  &
                      eps_lim=> r_cell_29,  &
                      kin_lim=> r_cell_30,  &
                      re_t   => r_cell_31
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2, iterator
  real    :: mag
  real    :: f_, f_w
  real    :: a11, a22, a33, a12, a13, a23
  real    :: e11, e22, e33, e12, e13, e23
  real    :: eps_h_11, eps_h_22, eps_h_33, eps_h_12, eps_h_13, eps_h_23
  real    :: f_eps, phi_km_2_n_k_n_m, eps_2_kin
  real    :: f_s,c_1,c_2
  real    :: p11, p22, p33, p12, p13, p23
  real    :: phi_ij_2_11, phi_ij_2_22, phi_ij_2_33
  real    :: phi_ij_2_12, phi_ij_2_13, phi_ij_2_23
  real    :: u_k_u_m_n_k_n_m
  real    :: prod_and_coriolis, phi_ij_2, phi_ij_1_w, phi_ij_2_w, eps_h, stress
  real    :: a_, a_2, a_3
  real    :: e_, e_2, e_3
  real    :: c_, c_1_w, c_2_w
  real    :: n1n1, n2n2, n3n3, n1n2, n1n3, n2n3
  real    :: eps_2_kin_wave, dissipation_2_term
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   Wall visc.    vis_wall [kg/(m*s)]  |                                       !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
! but dens > 1 mod. not applied here yet

  call Time_And_Length_Scale(grid) ! duplicated in Main_Pro

  !--------------------------------------------!
  !   preparations before filling source term  !
  !--------------------------------------------!

  do c = 1, grid % n_cells
    eps_lim(c) = max(eps % n(c), TINY) ! limited eps % n
    kin_lim(c) = max(kin % n(c), TINY) ! limited kin % n
  end do

  ! needed for eps_tot -> re_t
  call Grad_Mod_For_Phi(grid, kin % n, 1, kin_x,  .true.)  ! dK/dx
  call Grad_Mod_For_Phi(grid, kin % n, 2, kin_y,  .true.)  ! dK/dy
  call Grad_Mod_For_Phi(grid, kin % n, 3, kin_z,  .true.)  ! dK/dz

  call Grad_Mod_For_Phi(grid, kin_x,   1, kin_xx, .true.)  ! d^2 K / dx^2
  call Grad_Mod_For_Phi(grid, kin_y,   2, kin_yy, .true.)  ! d^2 K / dy^2
  call Grad_Mod_For_Phi(grid, kin_z,   3, kin_zz, .true.)  ! d^2 K / dz^2

  do c = 1, grid % n_cells
    eps_tot(c) = max(eps % n(c) + &
      0.5 * viscosity * (kin_xx(c) + kin_yy(c) + kin_zz(c)), TINY)
    ! page 142 Re_t
    re_t(c)  = kin % n(c)**2./(viscosity*eps_tot(c))
  end do

  if(name_phi .ne. 'EPS') then ! it is not eps eq.

    ! needed for n_i_n_j
    call Grad_Mod_For_Phi(grid, l_scale, 1, l_sc_x,.true.)
    call Grad_Mod_For_Phi(grid, l_scale, 2, l_sc_y,.true.)
    call Grad_Mod_For_Phi(grid, l_scale, 3, l_sc_z,.true.)

  else ! it is eps eq.

    ! for formula 2.19
    do c = 1, grid % n_cells
      kin_sq(c) = sqrt( kin % n(c) )
    end do

    call Grad_Mod_For_Phi(grid, kin_sq, 1, kin_x, .true.) ! dK/dx
    call Grad_Mod_For_Phi(grid, kin_sq, 2, kin_y, .true.) ! dK/dy
    call Grad_Mod_For_Phi(grid, kin_sq, 3, kin_z, .true.) ! dK/dz

    ! page 165 for term at -2*\nu(...)
    ! d^2U_k/dx_i dx_j can be sent outside this fucntion
    call Grad_Mod_For_Phi(grid, u % x,  1, u_xx, .true.)  ! d2U/dxdx
    call Grad_Mod_For_Phi(grid, u % y,  2, u_yy, .true.)  ! d2U/dydy
    call Grad_Mod_For_Phi(grid, u % z,  3, u_zz, .true.)  ! d2U/dzdz
    call Grad_Mod_For_Phi(grid, u % x,  2, u_xy, .true.)  ! d2U/dxdy
    call Grad_Mod_For_Phi(grid, u % x,  3, u_xz, .true.)  ! d2U/dxdz
    call Grad_Mod_For_Phi(grid, u % y,  3, u_yz, .true.)  ! d2U/dydz

    call Grad_Mod_For_Phi(grid, v % x,  1, v_xx, .true.)  ! d2V/dxdx
    call Grad_Mod_For_Phi(grid, v % y,  2, v_yy, .true.)  ! d2V/dydy
    call Grad_Mod_For_Phi(grid, v % z,  3, v_zz, .true.)  ! d2V/dzdz
    call Grad_Mod_For_Phi(grid, v % x,  2, v_xy, .true.)  ! d2V/dxdy
    call Grad_Mod_For_Phi(grid, v % x,  3, v_xz, .true.)  ! d2V/dxdz
    call Grad_Mod_For_Phi(grid, v % y,  3, v_yz, .true.)  ! d2V/dydz

    call Grad_Mod_For_Phi(grid, w % x,  1, w_xx, .true.)  ! d2W/dxdx
    call Grad_Mod_For_Phi(grid, w % y,  2, w_yy, .true.)  ! d2W/dydy
    call Grad_Mod_For_Phi(grid, w % z,  3, w_zz, .true.)  ! d2W/dzdz
    call Grad_Mod_For_Phi(grid, w % x,  2, w_xy, .true.)  ! d2W/dxdy
    call Grad_Mod_For_Phi(grid, w % x,  3, w_xz, .true.)  ! d2W/dxdz
    call Grad_Mod_For_Phi(grid, w % y,  3, w_yz, .true.)  ! d2W/dydz

    call Grad_Mod_For_Phi(grid, uu % n, 1, uu % x, .true.)  ! duu/dx
    call Grad_Mod_For_Phi(grid, uu % n, 2, uu % y, .true.)  ! duu/dy
    call Grad_Mod_For_Phi(grid, uu % n, 3, uu % z, .true.)  ! duu/dz

    call Grad_Mod_For_Phi(grid, vv % n, 1, vv % x, .true.)  ! dvv/dx
    call Grad_Mod_For_Phi(grid, vv % n, 2, vv % y, .true.)  ! dvv/dy
    call Grad_Mod_For_Phi(grid, vv % n, 3, vv % z, .true.)  ! dvv/dz

    call Grad_Mod_For_Phi(grid, ww % n, 1, ww % x, .true.)  ! dww/dx
    call Grad_Mod_For_Phi(grid, ww % n, 2, ww % y, .true.)  ! dww/dy
    call Grad_Mod_For_Phi(grid, ww % n, 3, ww % z, .true.)  ! dww/dz

    call Grad_Mod_For_Phi(grid, uv % n, 1, uv % x, .true.)  ! duv/dx
    call Grad_Mod_For_Phi(grid, uv % n, 2, uv % y, .true.)  ! duv/dy
    call Grad_Mod_For_Phi(grid, uv % n, 3, uv % z, .true.)  ! duv/dz

    call Grad_Mod_For_Phi(grid, uw % n, 1, uw % x, .true.)  ! duw/dx
    call Grad_Mod_For_Phi(grid, uw % n, 2, uw % y, .true.)  ! duw/dy
    call Grad_Mod_For_Phi(grid, uw % n, 3, uw % z, .true.)  ! duw/dz

    call Grad_Mod_For_Phi(grid, vw % n, 1, vw % x, .true.)  ! dvw/dx
    call Grad_Mod_For_Phi(grid, vw % n, 2, vw % y, .true.)  ! dvw/dy
    call Grad_Mod_For_Phi(grid, vw % n, 3, vw % z, .true.)  ! dvw/dz

  end if

  !------------------------!
  !   filling source term  !
  !------------------------!
  do  c = 1, grid % n_cells

    ! Epsilon over kinetic energy used almost 30 times in this loop
    eps_2_kin = eps % n(c) / kin_lim(c)

    ! P_k = 0.5 P_ii = - u_i u_k dU_i/dx_k
    p_kin(c) = -( uu % n(c) * u % x(c)  &
                + uv % n(c) * u % y(c)  &
                + uw % n(c) * u % z(c)  &
                + uv % n(c) * v % x(c)  &
                + vv % n(c) * v % y(c)  &
                + vw % n(c) * v % z(c)  &
                + uw % n(c) * w % x(c)  &
                + vw % n(c) * w % y(c)  &
                + ww % n(c) * w % z(c)  )

    if (name_phi .ne. 'EPS') then

      ! formula 2.4
      a11 = uu % n(c) / kin_lim(c) - TWO_THIRDS
      a22 = vv % n(c) / kin_lim(c) - TWO_THIRDS
      a33 = ww % n(c) / kin_lim(c) - TWO_THIRDS
      a12 = uv % n(c) / kin_lim(c)
      a13 = uw % n(c) / kin_lim(c)
      a23 = vw % n(c) / kin_lim(c)

      ! page 143 A2
      a_2   = a11**2. + a22**2. + a33**2.  &
        + 2.*(a12**2. + a13**2. + a23**2.  )

      ! page 143 A3
      a_3 = a11**3. + a22**3. + a33**3. &
        + 3*a12**2.*(a11+a22)           &
        + 3*a13**2.*(a11+a33)           &
        + 3*a23**2.*(a22+a33)           &
        + 6*a12*a13*a23

      ! page 143 A
      a_ = 1. - (9./8.)*(a_2-a_3)
      a_ = max(a_,0.)
      a_ = min(a_,1.)

      !---------------------------------------------------!
      !   iterative procedure to find e_, e_2, e_3, f_s   !
      !---------------------------------------------------!
      ! initial value for e_ = a_
      e_ = a_
      f_s = 1. - sqrt(a_) * e_**2.

      do iterator = 1, 6
        ! formula 2.14
        eps_h_11 = (1.-f_s) * TWO_THIRDS*eps % n(c) + f_s * uu % n(c)*eps_2_kin
        eps_h_22 = (1.-f_s) * TWO_THIRDS*eps % n(c) + f_s * vv % n(c)*eps_2_kin
        eps_h_33 = (1.-f_s) * TWO_THIRDS*eps % n(c) + f_s * ww % n(c)*eps_2_kin
        eps_h_12 =                                    f_s * uv % n(c)*eps_2_kin
        eps_h_13 =                                    f_s * uw % n(c)*eps_2_kin
        eps_h_23 =                                    f_s * vw % n(c)*eps_2_kin

        ! formula 2.4
        e11 = eps_h_11 / eps_lim(c) - TWO_THIRDS
        e22 = eps_h_22 / eps_lim(c) - TWO_THIRDS
        e33 = eps_h_33 / eps_lim(c) - TWO_THIRDS
        e12 = eps_h_12 / eps_lim(c)
        e13 = eps_h_13 / eps_lim(c)
        e23 = eps_h_23 / eps_lim(c)

        ! page 143 e2
        e_2   = e11**2. + e22**2. + e33**2.  &
          + 2.*(e12**2. + e13**2. + e23**2.  )

        ! page 143 e3
        e_3 = e11**3. + e22**3. + e33**3. &
          + 3*e12**2.*(e11+e22)           &
          + 3*e13**2.*(e11+e33)           &
          + 3*e23**2.*(e22+e33)           &
          + 6*e12*e13*e23

        ! page 143 E
        e_ = 1. - (9./8.) * (e_2 - e_3)
        e_ = max(e_, 0.) ! does not influence result
        e_ = min(e_, 1.) ! does not influence result
        ! page 143 f_s
        f_s = 1.-(a_**0.5*e_**2.)
      end do

      mag = max(l_sc_x(c)**2. + l_sc_y(c)**2. + l_sc_z(c)**2., TINY)

      ! normal unit (n_i never appears individually, only as n_i * n_j)
      n1n1 = l_sc_x(c)**2.       / mag
      n2n2 = l_sc_y(c)**2.       / mag
      n3n3 = l_sc_z(c)**2.       / mag
      n1n2 = l_sc_x(c)*l_sc_y(c) / mag
      n1n3 = l_sc_x(c)*l_sc_z(c) / mag
      n2n3 = l_sc_y(c)*l_sc_z(c) / mag

      ! page 165 f
      f_  = min((re_t(c)/150)**1.5, 1.)
      ! page 143 table (not used)
      ! f_d   = 1./(1. + 0.1*re_t)
      ! page 165 C
      c_  = 2.5*a_*min(0.6, a_2)**0.25*f_
      ! page 165 C_1
      c_1  = c_ + sqrt(a_)*e_**2.
      ! page 165 C_2
      c_2  = 0.8*sqrt(a_)
      ! page 165 C_1^w
      c_1_w  = max(1. - 0.7*c_, 0.3)
      ! page 165 C_2^w
      c_2_w  = min(a_,0.3)
      ! page 165 f_w
      f_w  = min(kin % n(c)**1.5/(2.5*eps % n(c)*grid % wall_dist(c)), 1.4)

      ! P_11 + G_11 (paper 2: formula C.1) ----- [ copied from Sources_Ebm ]
      p11 = -2.*(uu % n(c)*u % x(c) + uv % n(c)*u % y(c) + uw % n(c)*u % z(c)) &
            -2.*omega_y*2.*uw % n(c) + 2.*omega_z*2.*uv % n(c)

      p22 = -2.*(uv % n(c)*v % x(c) + vv % n(c)*v % y(c) + vw % n(c)*v % z(c)) &
            -2.*omega_x*2.*vw % n(c) + 2.*omega_z*2.*uw % n(c)

      p33 = -2.*(uw % n(c)*w % x(c) + vw % n(c)*w % y(c) + ww % n(c)*w % z(c)) &
            -2.*omega_x*2.*vw % n(c) + 2.*omega_y*2.*uw % n(c)

      p12 = &
        - uu % n(c)*v % x(c) - uv%n(c)*(v % y(c)+u % x(c)) - uw % n(c)*v % z(c)&
        - vv % n(c)*u % y(c) - vw % n(c)*u % z(c) &
        + 2.*omega_x*uw%n(c)-2.*omega_y*vw%n(c)+2.*omega_z*(vv%n(c)-uu%n(c))

      p13 = &
        - uu % n(c)*w % x(c) - uv%n(c)*w % y(c) - uw % n(c)*(w % z(c)+u % x(c))&
        - vw % n(c)*u % y(c) - ww % n(c)*u % z(c) &
        - 2.*omega_x*uv%n(c) - 2.*omega_y*(ww%n(c)-uu%n(c)) + 2.*omega_z*vw%n(c)

      p23 = &
        - uu % n(c)*w % x(c) - uv%n(c)*w % y(c) - uw % n(c)*(w % z(c)+u % x(c))&
        - vw % n(c)*u % y(c) - ww % n(c)*u % z(c) &
        - 2.*omega_x*(vv%n(c)-ww%n(c))+ 2.*omega_y*uv%n(c) - 2.*omega_z*uw%n(c)

      ! page 164 Phi_ij,2
      phi_ij_2_11 = - c_2 * (p11 - TWO_THIRDS*p_kin(c) )
      phi_ij_2_22 = - c_2 * (p22 - TWO_THIRDS*p_kin(c) )
      phi_ij_2_33 = - c_2 * (p33 - TWO_THIRDS*p_kin(c) )
      phi_ij_2_12 = - c_2 *  p12
      phi_ij_2_13 = - c_2 *  p13
      phi_ij_2_23 = - c_2 *  p23

      if(name_phi .eq. 'UU' .or. &
         name_phi .eq. 'VV' .or. &
         name_phi .eq. 'WW') then

        ! for formula page 164 Phi_ij,1^w
        u_k_u_m_n_k_n_m = uu % n(c)*n1n1 +   vv % n(c)*n2n2 +   ww % n(c)*n3n3 &
                      + 2*uv % n(c)*n1n2 + 2*uw % n(c)*n1n3 + 2*vw % n(c)*n2n3

        ! page 164 for Phi_ij,2^w
        phi_km_2_n_k_n_m =   phi_ij_2_11*n1n1  &
                         +   phi_ij_2_33*n3n3  &
                         +   phi_ij_2_22*n2n2  &
                         + 2*phi_ij_2_12*n1n2  &
                         + 2*phi_ij_2_13*n1n3  &
                         + 2*phi_ij_2_23*n2n3
      end if

      !---------------!
      !   uu stress   !
      !---------------!
      if(name_phi .eq. 'UU') then

        stress = uu % n(c)

        prod_and_coriolis = p11

        ! page 164 Phi_ij,2
        phi_ij_2 = phi_ij_2_11

        ! page 164 Phi_ij,1^w
        phi_ij_1_w = c_1_w * f_w * eps_2_kin * ( u_k_u_m_n_k_n_m - 3. *        &
          ( uu % n(c)*n1n1 + uv % n(c)*n1n2 + uw % n(c)*n1n3 )                 )

        ! page 164 Phi_ij,2^w
        phi_ij_2_w = c_2_w * f_w * ( phi_km_2_n_k_n_m - 3. *                   &
          ( phi_ij_2_11*n1n1 + phi_ij_2_12*n1n2 + phi_ij_2_13*n1n3 )           )

        eps_h = eps_h_11
      end if
      !---------------!
      !   vv stress   !
      !---------------!
      if(name_phi .eq. 'VV') then

        stress = vv % n(c)

        prod_and_coriolis = p22

        ! page 164 Phi_ij,2
        phi_ij_2 = phi_ij_2_22

        ! page 164 Phi_ij,1^w
        phi_ij_1_w = c_1_w * f_w * eps_2_kin * ( u_k_u_m_n_k_n_m - 3. *        &
          ( uv % n(c)*n1n2 + vv % n(c)*n2n2 + vw % n(c)*n2n3 )                 )

        ! page 164 Phi_ij,2^w
        phi_ij_2_w = c_2_w * f_w * ( phi_km_2_n_k_n_m - 3. *                   &
          ( phi_ij_2_12*n1n2 + phi_ij_2_22*n2n2 + phi_ij_2_23*n2n3 )           )

        eps_h = eps_h_22
      end if
      !---------------!
      !   ww stress   !
      !---------------!
      if(name_phi .eq. 'WW') then

        stress = ww % n(c)

        prod_and_coriolis = p33

        ! page 164 Phi_ij,2
        phi_ij_2 = phi_ij_2_33

        ! page 164 Phi_ij,1^w
        phi_ij_1_w = c_1_w * f_w * eps_2_kin * ( u_k_u_m_n_k_n_m - 3. *        &
          ( uw % n(c)*n1n3 + vw % n(c)*n2n3 + ww % n(c)*n3n3 )                 )

        ! page 164 Phi_ij,2^w
        phi_ij_2_w = c_2_w * f_w * ( phi_km_2_n_k_n_m - 3. *                   &
          ( phi_ij_2_13*n1n3 + phi_ij_2_23*n2n3 + phi_ij_2_33*n3n3 )           )

        eps_h = eps_h_33
      end if
      !---------------!
      !   uv stress   !
      !---------------!
      if(name_phi .eq. 'UV') then

        stress = uv % n(c)

        prod_and_coriolis = p12

        ! page 164 Phi_ij,2
        phi_ij_2 = phi_ij_2_12

        ! page 164 Phi_ij,1^w
        phi_ij_1_w = c_1_w * f_w * eps_2_kin * ( - 1.5 ) * (    &
          uu % n(c)*n1n2 + uv % n(c)*n2n2 + uw % n(c)*n2n3 +    &
          uv % n(c)*n1n1 + vv % n(c)*n1n2 + vw % n(c)*n1n3      )

        ! page 164 Phi_ij,2^w
        phi_ij_2_w = c_2_w * f_w * ( - 1.5 ) * (                   &
          phi_ij_2_11*n1n2 + phi_ij_2_12*n2n2 + phi_ij_2_13*n2n3 + &
          phi_ij_2_12*n1n1 + phi_ij_2_22*n1n2 + phi_ij_2_23*n1n3   )

        eps_h = eps_h_12

      end if
      !---------------!
      !   uw stress   !
      !---------------!
      if(name_phi .eq. 'UW') then

        stress = uw % n(c)

        prod_and_coriolis = p13

        ! page 164 Phi_ij,2
        phi_ij_2 = phi_ij_2_13

        ! page 164 Phi_ij,1^w
        phi_ij_1_w = c_1_w * f_w * eps_2_kin * ( - 1.5 ) * (    &
          uu % n(c)*n1n3 + uv % n(c)*n2n3 + uw % n(c)*n3n3 +    &
          uw % n(c)*n1n1 + vw % n(c)*n1n2 + ww % n(c)*n1n3      )

        ! page 164 Phi_ij,2^w
        phi_ij_2_w = c_2_w * f_w * ( - 1.5 ) * (                   &
          phi_ij_2_11*n1n3 + phi_ij_2_12*n2n3 + phi_ij_2_13*n3n3 + &
          phi_ij_2_13*n1n1 + phi_ij_2_23*n1n2 + phi_ij_2_33*n1n3   )

        eps_h = eps_h_13

      end if
      !---------------!
      !   vw stress   !
      !---------------!
      if(name_phi .eq. 'VW') then

        stress = vw % n(c)

        prod_and_coriolis = p23

        ! page 164 Phi_ij,1^w
        phi_ij_1_w = c_1_w * f_w * eps_2_kin * ( - 1.5 ) * (    &
          uw % n(c)*n1n2 + vw % n(c)*n2n2 + ww % n(c)*n2n3 +    &
          uv % n(c)*n1n3 + vv % n(c)*n2n3 + vw % n(c)*n3n3      )

        ! page 164 Phi_ij,2^w
        phi_ij_2_w = c_2_w * f_w * ( - 1.5 ) * (                   &
          phi_ij_2_13*n1n2 + phi_ij_2_23*n2n2 + phi_ij_2_33*n2n3 + &
          phi_ij_2_12*n1n3 + phi_ij_2_22*n2n3 + phi_ij_2_23*n3n3   )

        ! page 164 Phi_ij,2
        phi_ij_2 = phi_ij_2_23

        eps_h = eps_h_23

      end if

      !-------------------------------------!
      !   repeating part for all stresses   !
      !-------------------------------------!
      ! limit stress
      if (stress .lt. 0) stress = min(stress,-TINY)
      if (stress .ge. 0) stress = max(stress, TINY)

      ! term depending on Kronecker symbol 
      if(name_phi .eq. 'UU' .or. &
         name_phi .eq. 'VV' .or. &
         name_phi .eq. 'WW') then

        b(c) = b(c) + grid % vol(c) *  (   &
          c_1 * eps % n(c)*TWO_THIRDS      & ! part of Phi_ij_1
          - (1.-f_s)*TWO_THIRDS*eps % n(c) ) ! part of eps_ij^h

      end if

      b(c) = b(c) + grid % vol(c) *  (     &
        + max(prod_and_coriolis, 0.)       & ! P_ij + G_ij, if > 0
        + max(phi_ij_2, 0.)                & ! Phi_ij_2   , if > 0
        + max(phi_ij_1_w, 0.)              & ! Phi_ij_1_w , if > 0
        + max(phi_ij_2_w, 0.)              ) ! Phi_ij_2_w , if > 0

      A % val(A % dia(c)) = A % val(A % dia(c)) + grid % vol(c) * ( &
          c_1 * eps_2_kin            & ! part of Phi_ij_1
        + f_s * eps_2_kin            & ! part of eps_ij^h / u_iu_j
        + (                          &
        - min(prod_and_coriolis, 0.) & ! (P_ij + G_ij) / u_iu_j, if < 0
        - min(phi_ij_2, 0.)          & ! (Phi_ij_2)    / u_iu_j, if < 0
        - min(phi_ij_1_w, 0.)        & ! (Phi_ij_1_w)  / u_iu_j, if < 0
        - min(phi_ij_2_w, 0.)        & ! (Phi_ij_2_w)  / u_iu_j, if < 0
          ) / stress                                              )
  !----------------------!
  !   Epsilon equation   !
  !----------------------!

    else ! it is eps eq.
    ! page 165 second term
      dissipation_2_term = - 2 * viscosity * (                                 &

        uu % x(c) * u_xx(c) + uv % x(c) * v_xx(c) + uw % x(c) * w_xx(c)        &
      + uv % x(c) * u_xy(c) + vv % x(c) * v_xy(c) + vw % x(c) * w_xy(c)        &
      + uw % x(c) * u_xz(c) + vw % x(c) * v_xz(c) + ww % x(c) * w_xz(c)        &

      + uu % y(c) * u_xy(c) + uv % y(c) * v_xy(c) + uw % y(c) * w_xy(c)        &
      + uv % y(c) * u_yy(c) + vv % y(c) * v_yy(c) + vw % y(c) * w_yy(c)        &
      + uw % y(c) * u_yz(c) + vw % y(c) * v_yz(c) + ww % y(c) * w_yz(c)        &

      + uu % z(c) * u_xz(c) + uv % z(c) * v_xz(c) + uw % z(c) * w_xz(c)        &
      + uv % z(c) * u_yz(c) + vv % z(c) * v_yz(c) + vw % z(c) * w_yz(c)        &
      + uw % z(c) * u_zz(c) + vw % z(c) * v_zz(c) + ww % z(c) * w_zz(c)        &

      + c_3e * t_scale(c) * (                                                  &
        uu % x(c)*( u % x(c)*u_xx(c) + v % x(c)*v_xx(c) + w % x(c)*w_xx(c) )   &
      + uu % y(c)*( u % x(c)*u_xy(c) + v % x(c)*v_xy(c) + w % x(c)*w_xy(c) )   &
      + uu % z(c)*( u % x(c)*u_xz(c) + v % x(c)*v_xz(c) + w % x(c)*w_xz(c) )   &

      + uv % x(c)*( u % y(c)*u_xx(c) + v % y(c)*v_xx(c) + w % y(c)*w_xx(c) )   &
      + uv % y(c)*( u % y(c)*u_xy(c) + v % y(c)*v_xy(c) + w % y(c)*w_xy(c) )   &
      + uv % z(c)*( u % y(c)*u_xz(c) + v % y(c)*v_xz(c) + w % y(c)*w_xz(c) )   &

      + uw % x(c)*( u % z(c)*u_xx(c) + v % z(c)*v_xx(c) + w % z(c)*w_xx(c) )   &
      + uw % y(c)*( u % z(c)*u_xy(c) + v % z(c)*v_xy(c) + w % z(c)*w_xy(c) )   &
      + uw % z(c)*( u % z(c)*u_xz(c) + v % z(c)*v_xz(c) + w % z(c)*w_xz(c) )   &

      + uv % x(c)*( u % x(c)*u_xy(c) + v % x(c)*v_xy(c) + w % x(c)*w_xy(c) )   &
      + uv % y(c)*( u % x(c)*u_yy(c) + v % x(c)*v_yy(c) + w % x(c)*w_yy(c) )   &
      + uv % z(c)*( u % x(c)*u_yz(c) + v % x(c)*v_yz(c) + w % x(c)*w_yz(c) )   &

      + vv % x(c)*( u % y(c)*u_xy(c) + v % y(c)*v_xy(c) + w % y(c)*w_xy(c) )   &
      + vv % y(c)*( u % y(c)*u_yy(c) + v % y(c)*v_yy(c) + w % y(c)*w_yy(c) )   &
      + vv % z(c)*( u % y(c)*u_yz(c) + v % y(c)*v_yz(c) + w % y(c)*w_yz(c) )   &

      + vw % x(c)*( u % z(c)*u_xy(c) + v % z(c)*v_xy(c) + w % z(c)*w_xy(c) )   &
      + vw % y(c)*( u % z(c)*u_yy(c) + v % z(c)*v_yy(c) + w % z(c)*w_yy(c) )   &
      + vw % z(c)*( u % z(c)*u_yz(c) + v % z(c)*v_yz(c) + w % z(c)*w_yz(c) )   &

      + uw % x(c)*( u % x(c)*u_xz(c) + v % x(c)*v_xz(c) + w % x(c)*w_xz(c) )   &
      + uw % y(c)*( u % x(c)*u_yz(c) + v % x(c)*v_yz(c) + w % x(c)*w_yz(c) )   &
      + uw % z(c)*( u % x(c)*u_zz(c) + v % x(c)*v_zz(c) + w % x(c)*w_zz(c) )   &

      + vw % x(c)*( u % y(c)*u_xz(c) + v % y(c)*v_xz(c) + w % y(c)*w_xz(c) )   &
      + vw % y(c)*( u % y(c)*u_yz(c) + v % y(c)*v_yz(c) + w % y(c)*w_yz(c) )   &
      + vw % z(c)*( u % y(c)*u_zz(c) + v % y(c)*v_zz(c) + w % y(c)*w_zz(c) )   &

      + ww % x(c)*( u % z(c)*u_xz(c) + v % z(c)*v_xz(c) + w % z(c)*w_xz(c) )   &
      + ww % y(c)*( u % z(c)*u_yz(c) + v % z(c)*v_yz(c) + w % z(c)*w_yz(c) )   &
      + ww % z(c)*( u % z(c)*u_zz(c) + v % z(c)*v_zz(c) + w % z(c)*w_zz(c) )   &
                                         )                                     &
                                                                               )
      ! page 165 f_eps
      f_eps = 1. - (1. - 1.4/c_2e)*exp(-(re_t(c)/6.)**2.)

      ! formula 2.19, second term, other part is in A
      eps_2_kin_wave = viscosity * (kin_x(c)**2. + kin_y(c)**2. + kin_z(c)**2.)

      ! page 165, diss. eq., second term
      b(c) = b(c) + grid % vol(c) * (                             &
        dissipation_2_term                                        &
        + c_2e * f_eps * eps % n(c) * eps_2_kin_wave / kin_lim(c) )

      A % val(A % dia(c)) = A % val(A % dia(c)) + grid % vol(c) * ( &
        + c_1e * p_kin(c) / kin_lim(c)  & ! page 165, diss. eq., first term
        + c_2e * f_eps * eps % n(c) / kin_lim(c) ) ! diss. eq., third last
    end if
  end do

  ! formula 2.19 fixes boundary conditions
  if(name_phi .ne. 'EPS') then ! it is not eps eq.
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! Calculate a values of dissipation  on wall
      if(c2 < 0 ) then
        if (Grid_Mod_Bnd_Cond_Type(grid,c2)  .ne. BUFFER ) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
             Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

            eps % n(c2) = viscosity*(kin_x(c1)**2 + kin_y(c1)**2 + kin_z(c1)**2)
          end if ! end if of BC=wall
        end if ! end if of c2<0
      end if ! end if of c2<0
    end do
  end if

  end subroutine