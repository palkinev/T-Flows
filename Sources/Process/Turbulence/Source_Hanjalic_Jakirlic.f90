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
                      U_xx   => r_cell_10,  &
                      U_yy   => r_cell_11,  &
                      U_zz   => r_cell_12,  &
                      U_xy   => r_cell_13,  &
                      U_xz   => r_cell_14,  &
                      U_yz   => r_cell_15,  &
                      V_xx   => r_cell_16,  &
                      V_yy   => r_cell_17,  &
                      V_zz   => r_cell_18,  &
                      V_xy   => r_cell_19,  &
                      V_xz   => r_cell_20,  &
                      V_yz   => r_cell_21,  &
                      W_xx   => r_cell_22,  &
                      W_yy   => r_cell_23,  &
                      W_zz   => r_cell_24,  &
                      W_xy   => r_cell_25,  &
                      W_xz   => r_cell_26,  &
                      W_yz   => r_cell_27,  &
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
  real    :: phi_ij_1_w_11, phi_ij_1_w_22, phi_ij_1_w_33
  real    :: phi_ij_1_w_12, phi_ij_1_w_13, phi_ij_1_w_23
  real    :: phi_ij_2_w_11, phi_ij_2_w_22, phi_ij_2_w_33
  real    :: phi_ij_2_w_12, phi_ij_2_w_13, phi_ij_2_w_23
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

  !--------------------------------------------!
  !   preparations before filling source term  !
  !--------------------------------------------!

  do c = 1, grid % n_cells
    kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)), TINY)
    eps_lim(c) = eps % n(c) + TINY ! limited eps to divide by
    kin_lim(c) = kin % n(c) + TINY ! limited kin to divide by
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

    do c = 1, grid % n_cells
      l_scale(c) = kin % n(c)**1.5/eps_lim(c)
      t_scale(c) = kin % n(c)/eps_lim(c)
    end do

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
    call Grad_Mod_For_Phi(grid, U % x,  1, U_xx, .true.)  ! d2U/dxdx
    call Grad_Mod_For_Phi(grid, U % y,  2, U_yy, .true.)  ! d2U/dydy
    call Grad_Mod_For_Phi(grid, U % z,  3, U_zz, .true.)  ! d2U/dzdz
    call Grad_Mod_For_Phi(grid, U % x,  2, U_xy, .true.)  ! d2U/dxdy
    call Grad_Mod_For_Phi(grid, U % x,  3, U_xz, .true.)  ! d2U/dxdz
    call Grad_Mod_For_Phi(grid, U % y,  3, U_yz, .true.)  ! d2U/dydz

    call Grad_Mod_For_Phi(grid, V % x,  1, V_xx, .true.)  ! d2V/dxdx
    call Grad_Mod_For_Phi(grid, V % y,  2, V_yy, .true.)  ! d2V/dydy
    call Grad_Mod_For_Phi(grid, V % z,  3, V_zz, .true.)  ! d2V/dzdz
    call Grad_Mod_For_Phi(grid, V % x,  2, V_xy, .true.)  ! d2V/dxdy
    call Grad_Mod_For_Phi(grid, V % x,  3, V_xz, .true.)  ! d2V/dxdz
    call Grad_Mod_For_Phi(grid, V % y,  3, V_yz, .true.)  ! d2V/dydz

    call Grad_Mod_For_Phi(grid, W % x,  1, W_xx, .true.)  ! d2W/dxdx
    call Grad_Mod_For_Phi(grid, W % y,  2, W_yy, .true.)  ! d2W/dydy
    call Grad_Mod_For_Phi(grid, W % z,  3, W_zz, .true.)  ! d2W/dzdz
    call Grad_Mod_For_Phi(grid, W % x,  2, W_xy, .true.)  ! d2W/dxdy
    call Grad_Mod_For_Phi(grid, W % x,  3, W_xz, .true.)  ! d2W/dxdz
    call Grad_Mod_For_Phi(grid, W % y,  3, W_yz, .true.)  ! d2W/dydz

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
    eps_2_kin = eps % n(c) / (kin % n(c) + TINY)

    ! P_k = 0.5 P_ii = - u_i u_k dU_i/dx_k
    p_kin(c) =     -( uu % n(c) * u % x(c)  &
                    + uv % n(c) * u % y(c)  &
                    + uw % n(c) * u % z(c)  &
                    + uv % n(c) * v % x(c)  &
                    + vv % n(c) * v % y(c)  &
                    + vw % n(c) * v % z(c)  &
                    + uw % n(c) * w % x(c)  &
                    + vw % n(c) * w % y(c)  &
                    + ww % n(c) * w % z(c)  )

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
      e_ = max(e_, 0.)
      e_ = min(e_, 1.)
      ! page 143 f_s
      f_s = 1.-(a_**0.5*e_**2.)
    end do

    if (name_phi .ne. 'EPS') then

      mag = max(l_sc_x(c)**2. + l_sc_y(c)**2. + l_sc_z(c)**2., TINY)

      ! normal unit (n_i never appears individually, only as n_i * n_j)
      n1n1 = l_sc_x(c)**2.       / mag
      n2n2 = l_sc_y(c)**2.       / mag
      n3n3 = l_sc_z(c)**2.       / mag
      n1n2 = l_sc_x(c)*l_sc_y(c) / mag
      n1n3 = l_sc_x(c)*l_sc_z(c) / mag
      n2n3 = l_sc_y(c)*l_sc_z(c) / mag

      ! for formula page 164 Phi_ij,1^w
      u_k_u_m_n_k_n_m = uu % n(c)*n1n1 +    vv % n(c)*n2n2 +    ww % n(c)*n3n3 &
                   + 2.*uv % n(c)*n1n2 + 2.*uw % n(c)*n1n3 + 2.*vw % n(c)*n2n3

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
      phi_ij_2_11 = -c_2 * (p11 - TWO_THIRDS*p_kin(c) )
      phi_ij_2_22 = -c_2 * (p22 - TWO_THIRDS*p_kin(c) )
      phi_ij_2_33 = -c_2 * (p33 - TWO_THIRDS*p_kin(c) )
      phi_ij_2_12 = -c_2 * p12
      phi_ij_2_13 = -c_2 * p13
      phi_ij_2_23 = -c_2 * p23

      ! page 164 Phi_ij,1^w
      phi_ij_1_w_11 = c_1_w * f_w * eps_2_kin * (              &
        u_k_u_m_n_k_n_m - 3. *                                 &
           ( uu % n(c)*n1n1 + uv % n(c)*n1n2 + uw % n(c)*n1n3) )

      phi_ij_1_w_22 = c_1_w * f_w * eps_2_kin * (              &
        u_k_u_m_n_k_n_m - 3. *                                 &
           ( uv % n(c)*n1n2 + vv % n(c)*n2n2 + vw % n(c)*n2n3) )

      phi_ij_1_w_33 = c_1_w * f_w * eps_2_kin * (              &
        u_k_u_m_n_k_n_m - 3. *                                 &
           ( uw % n(c)*n1n3 + vw % n(c)*n2n3 + ww % n(c)*n3n3) )

      phi_ij_1_w_12 = c_1_w * f_w * eps_2_kin * ( - 1.5 ) * ( &
        uu % n(c)*n1n2 + uv % n(c)*n2n2 + uw % n(c)*n2n3 +    &
        uv % n(c)*n1n1 + vv % n(c)*n1n2 + vw % n(c)*n1n3      )

      phi_ij_1_w_13 = c_1_w * f_w * eps_2_kin * ( - 1.5 ) * ( &
        uu % n(c)*n1n3 + uv % n(c)*n2n3 + uw % n(c)*n3n3 +    &
        uw % n(c)*n1n1 + vw % n(c)*n1n2 + ww % n(c)*n1n3      )

      phi_ij_1_w_23 = c_1_w * f_w * eps_2_kin * ( - 1.5 ) * ( &
        uv % n(c)*n1n3 + vv % n(c)*n2n3 + vw % n(c)*n3n3 +    &
        uw % n(c)*n1n2 + vw % n(c)*n2n2 + ww % n(c)*n2n3      )

      ! page 164 for Phi_ij,2^w
      phi_km_2_n_k_n_m =   phi_ij_2_11*n1n1  &
                       +   phi_ij_2_33*n3n3  &
                       +   phi_ij_2_22*n2n2  &
                       + 2*phi_ij_2_12*n1n2  &
                       + 2*phi_ij_2_13*n1n3  &
                       + 2*phi_ij_2_23*n2n3

      ! page 164 Phi_ij,2^w
      phi_ij_2_w_11 = c_2_w * f_w *                                      &
        ( phi_km_2_n_k_n_m - 1.5  *                                      &
          2. * ( phi_ij_2_11*n1n1 + phi_ij_2_12*n1n2 + phi_ij_2_13*n1n3) )

      phi_ij_2_w_22 = c_2_w * f_w *                                      &
        ( phi_km_2_n_k_n_m - 1.5  *                                      &
          2. * ( phi_ij_2_12*n1n2 + phi_ij_2_22*n2n2 + phi_ij_2_23*n2n3) )

      phi_ij_2_w_33 = c_2_w * f_w *                                      &
        ( phi_km_2_n_k_n_m - 1.5  *                                      &
          2. * ( phi_ij_2_13*n1n3 + phi_ij_2_23*n2n3 + phi_ij_2_33*n3n3) )

      phi_ij_2_w_12 = c_2_w * f_w * ( - 1.5 ) * (                &
        phi_ij_2_11*n1n2 + phi_ij_2_12*n2n2 + phi_ij_2_13*n2n3 + &
        phi_ij_2_12*n1n1 + phi_ij_2_22*n1n2 + phi_ij_2_23*n1n3   )

      phi_ij_2_w_13 = c_2_w * f_w * ( - 1.5 ) * (                &
        phi_ij_2_11*n1n3 + phi_ij_2_12*n2n3 + phi_ij_2_13*n3n3 + &
        phi_ij_2_13*n1n1 + phi_ij_2_23*n1n2 + phi_ij_2_33*n1n3   )

      phi_ij_2_w_23 = c_2_w * f_w * ( - 1.5 ) * (                &
        phi_ij_2_13*n1n2 + phi_ij_2_23*n2n2 + phi_ij_2_33*n2n3 + &
        phi_ij_2_12*n1n3 + phi_ij_2_22*n2n3 + phi_ij_2_23*n3n3   )

      !---------------!
      !   uu stress   !
      !---------------!
      if(name_phi .eq. 'UU') then
        ! limited stress
        stress = max(uu % n(c), TINY)

        prod_and_coriolis = p11

        phi_ij_2   = phi_ij_2_11
        phi_ij_1_w = phi_ij_1_w_11
        phi_ij_2_w = phi_ij_2_w_11

        eps_h = eps_h_11
      end if
      !---------------!
      !   vv stress   !
      !---------------!
      if(name_phi .eq. 'VV') then
        ! limited stress
        stress = max(vv % n(c), TINY)

        prod_and_coriolis = p22

        phi_ij_2   = phi_ij_2_22
        phi_ij_1_w = phi_ij_1_w_22
        phi_ij_2_w = phi_ij_2_w_22

        eps_h = eps_h_22
      end if
      !---------------!
      !   ww stress   !
      !---------------!
      if(name_phi .eq. 'WW') then
        ! limited stress
        stress = max(ww % n(c), TINY)

        prod_and_coriolis = p33

        phi_ij_2   = phi_ij_2_33
        phi_ij_1_w = phi_ij_1_w_33
        phi_ij_2_w = phi_ij_2_w_33

        eps_h = eps_h_33
      end if
      !---------------!
      !   uv stress   !
      !---------------!
      if(name_phi .eq. 'UV') then
        ! limited stress
        stress = max(uv % n(c), TINY)

        prod_and_coriolis = p12

        phi_ij_2   = phi_ij_2_12
        phi_ij_1_w = phi_ij_1_w_12
        phi_ij_2_w = phi_ij_2_w_12

        eps_h = eps_h_12

        b(c) = b(c) + grid % vol(c) * ( &
          - c_1 * eps % n(c)*TWO_THIRDS ) ! part of Phi_ij_1 to discard Kroneker
      end if
      !---------------!
      !   uw stress   !
      !---------------!
      if(name_phi .eq. 'UW') then
        ! limited stress
        stress = max(uw % n(c), TINY)

        prod_and_coriolis = p13

        phi_ij_2   = phi_ij_2_13
        phi_ij_1_w = phi_ij_1_w_13
        phi_ij_2_w = phi_ij_2_w_13

        eps_h = eps_h_13

        b(c) = b(c) + grid % vol(c) * ( &
          - c_1 * eps % n(c)*TWO_THIRDS ) ! part of Phi_ij_1 to discard Kroneker
      end if
      !---------------!
      !   vw stress   !
      !---------------!
      if(name_phi .eq. 'VW') then
        ! limited stress
        stress = max(vw % n(c), TINY)

        prod_and_coriolis = p23

        phi_ij_2   = phi_ij_2_23
        phi_ij_1_w = phi_ij_1_w_23
        phi_ij_2_w = phi_ij_2_w_23

        eps_h = eps_h_23

        b(c) = b(c) + grid % vol(c) * ( &
          - c_1 * eps % n(c)*TWO_THIRDS ) ! part of Phi_ij_1 to discard Kroneker
      end if

      !-------------------------------------!
      !   repeating part for all stresses   !
      !-------------------------------------!
      b(c) = b(c) + grid % vol(c) *  ( &
          c_1 * eps % n(c)*TWO_THIRDS  & ! part of Phi_ij_1
        + max(prod_and_coriolis, 0.)   & ! P_ij + G_ij, if > 0
        + max(phi_ij_2, 0.)            & ! Phi_ij_2   , if > 0
        + max(phi_ij_1_w, 0.)          & ! Phi_ij_1_w , if > 0
        + max(phi_ij_2_w, 0.)          ) ! Phi_ij_2_w , if > 0

      A % val(A % dia(c)) = A % val(A % dia(c)) + grid % vol(c) * ( &
          c_1 * eps_2_kin            & ! part of Phi_ij_1
        + eps_h                      & ! eps_ij^h / u_iu_j
        - min(prod_and_coriolis, 0.) & ! (P_ij + G_ij) / u_iu_j, if < 0
        - min(phi_ij_2, 0.)          & ! (Phi_ij_2)    / u_iu_j, if < 0
        - min(phi_ij_1_w, 0.)        & ! (Phi_ij_1_w)  / u_iu_j, if < 0
        - min(phi_ij_2_w, 0.)        & ! (Phi_ij_2_w)  / u_iu_j, if < 0
        ) / stress
  !----------------------!
  !   Epsilon equation   !
  !----------------------!

    else ! it is eps eq.
    ! page 165 second term
      dissipation_2_term = - 2 * viscosity * (                                 &
                                                                               !
        uu % x(c) * U_xx(c) + uv % x(c) * V_xx(c) + uw % x(c) * W_xx(c)        &
      + uv % x(c) * U_xy(c) + vv % x(c) * V_xy(c) + vw % x(c) * W_xy(c)        &
      + uw % x(c) * U_xz(c) + vw % x(c) * V_xz(c) + ww % x(c) * W_xz(c)        &
                                                                               !
      + uu % y(c) * U_xy(c) + uv % y(c) * V_xy(c) + uw % y(c) * W_xy(c)        &
      + uv % y(c) * U_yy(c) + vv % y(c) * V_yy(c) + vw % y(c) * W_yy(c)        &
      + uw % y(c) * U_yz(c) + vw % y(c) * V_yz(c) + ww % y(c) * W_yz(c)        &
                                                                               !
      + uu % z(c) * U_xz(c) + uv % z(c) * V_xz(c) + uw % z(c) * W_xz(c)        &
      + uv % z(c) * U_yz(c) + vv % z(c) * V_yz(c) + vw % z(c) * W_yz(c)        &
      + uw % z(c) * U_zz(c) + vw % z(c) * V_zz(c) + ww % z(c) * W_zz(c)        &
                                                                               !
      + c_3e * kin % n(c) / eps_lim(c) * (                                     &
        uu % x(c)*( U % x(c)*U_xx(c) + V % x(c)*V_xx(c) + W % x(c)*W_xx(c) )   &
      + uu % y(c)*( U % x(c)*U_xy(c) + V % x(c)*V_xy(c) + W % x(c)*W_xy(c) )   &
      + uu % z(c)*( U % x(c)*U_xz(c) + V % x(c)*V_xz(c) + W % x(c)*W_xz(c) )   &
                                                                               !
      + uv % x(c)*( U % y(c)*U_xx(c) + V % y(c)*V_xx(c) + W % y(c)*W_xx(c) )   &
      + uv % y(c)*( U % y(c)*U_xy(c) + V % y(c)*V_xy(c) + W % y(c)*W_xy(c) )   &
      + uv % z(c)*( U % y(c)*U_xz(c) + V % y(c)*V_xz(c) + W % y(c)*W_xz(c) )   &
                                                                               !
      + uw % x(c)*( U % z(c)*U_xx(c) + V % z(c)*V_xx(c) + W % z(c)*W_xx(c) )   &
      + uw % y(c)*( U % z(c)*U_xy(c) + V % z(c)*V_xy(c) + W % z(c)*W_xy(c) )   &
      + uw % z(c)*( U % z(c)*U_xz(c) + V % z(c)*V_xz(c) + W % z(c)*W_xz(c) )   &
                                                                               !
      + uv % x(c)*( U % x(c)*U_xy(c) + V % x(c)*V_xy(c) + W % x(c)*W_xy(c) )   &
      + uv % y(c)*( U % x(c)*U_yy(c) + V % x(c)*V_yy(c) + W % x(c)*W_yy(c) )   &
      + uv % z(c)*( U % x(c)*U_yz(c) + V % x(c)*V_yz(c) + W % x(c)*W_yz(c) )   &
                                                                               !
      + vv % x(c)*( U % y(c)*U_xy(c) + V % y(c)*V_xy(c) + W % y(c)*W_xy(c) )   &
      + vv % y(c)*( U % y(c)*U_yy(c) + V % y(c)*V_yy(c) + W % y(c)*W_yy(c) )   &
      + vv % z(c)*( U % y(c)*U_yz(c) + V % y(c)*V_yz(c) + W % y(c)*W_yz(c) )   &
                                                                               !
      + vw % x(c)*( U % z(c)*U_xy(c) + V % z(c)*V_xy(c) + W % z(c)*W_xy(c) )   &
      + vw % y(c)*( U % z(c)*U_yy(c) + V % z(c)*V_yy(c) + W % z(c)*W_yy(c) )   &
      + vw % z(c)*( U % z(c)*U_yz(c) + V % z(c)*V_yz(c) + W % z(c)*W_yz(c) )   &
                                                                               !
      + uw % x(c)*( U % x(c)*U_xz(c) + V % x(c)*V_xz(c) + W % x(c)*W_xz(c) )   &
      + uw % y(c)*( U % x(c)*U_yz(c) + V % x(c)*V_yz(c) + W % x(c)*W_yz(c) )   &
      + uw % z(c)*( U % x(c)*U_zz(c) + V % x(c)*V_zz(c) + W % x(c)*W_zz(c) )   &
                                                                               !
      + vw % x(c)*( U % y(c)*U_xz(c) + V % y(c)*V_xz(c) + W % y(c)*W_xz(c) )   &
      + vw % y(c)*( U % y(c)*U_yz(c) + V % y(c)*V_yz(c) + W % y(c)*W_yz(c) )   &
      + vw % z(c)*( U % y(c)*U_zz(c) + V % y(c)*V_zz(c) + W % y(c)*W_zz(c) )   &
                                                                               !
      + ww % x(c)*( U % z(c)*U_xz(c) + V % z(c)*V_xz(c) + W % z(c)*W_xz(c) )   &
      + ww % y(c)*( U % z(c)*U_yz(c) + V % z(c)*V_yz(c) + W % z(c)*W_yz(c) )   &
      + ww % z(c)*( U % z(c)*U_zz(c) + V % z(c)*V_zz(c) + W % z(c)*W_zz(c) )   &
                                         )                                     &
                                                                               )
      ! page 165 f_eps
      f_eps = 1. - (1. - 1.4/c_2e)*exp(-(re_t(c)/6.)**2.)

      ! formula 2.19
      eps_2_kin_wave = eps % n(c) - &
        viscosity * (kin_x(c)**2. + kin_y(c)**2. + kin_z(c)**2.)

      ! page 165, diss. eq., second term
      b(c) = b(c) + grid % vol(c) * dissipation_2_term

      A % val(A % dia(c)) = A % val(A % dia(c)) + grid % vol(c) * ( &
        - c_1e * p_kin(c) / kin_lim(c)  & ! page 165, diss. eq., first term
        + c_2e * f_eps * eps_2_kin_wave ) ! page 165, diss. eq., third last
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

            eps % n(c2) = viscosity*(kin_x(c2)**2 + kin_y(c2)**2 + kin_z(c2)**2)
          end if ! end if of BC=wall
        end if ! end if of c2<0
      end if ! end if of c2<0
    end do
  end if

  end subroutine