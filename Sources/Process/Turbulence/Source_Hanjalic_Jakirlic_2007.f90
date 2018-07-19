!==============================================================================!
  subroutine Source_Hanjalic_Jakirlic_2007(grid, name_phi)
!------------------------------------------------------------------------------!
!   Calculate source terms for transport equations for Re stresses and         !
!   dissipation for Hanjalic-Jakirlic model.                                   !
!   Following paper "Near-wall, Reynolds-stress model calculations of          !
!   transonic flow configurations relevant to aircraft aerodynamics"           !
!   DOI: doi:10.1016/j.ijheatfluidflow.2007.04.001                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Work_Mod, only: l_sc_x               => r_cell_01,  &
                      l_sc_y               => r_cell_02,  &
                      l_sc_z               => r_cell_03,  &
                      kin_x                => r_cell_04,  &
                      kin_y                => r_cell_05,  &
                      kin_z                => r_cell_06,  &
                      kin_xx               => r_cell_07,  &
                      kin_yy               => r_cell_08,  &
                      kin_zz               => r_cell_09,  &
                      u_xx                 => r_cell_10,  &
                      u_yy                 => r_cell_11,  &
                      u_zz                 => r_cell_12,  &
                      u_xy                 => r_cell_13,  &
                      u_xz                 => r_cell_14,  &
                      u_yz                 => r_cell_15,  &
                      kin_sq               => r_cell_16,  &
                      eps_lim              => r_cell_17,  &
                      kin_lim              => r_cell_18,  &
                      re_t                 => r_cell_19,  &
                      dissipation_2_term   => r_cell_20
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2, iterator
  real    :: mag
  real    :: p11, p22, p33, p12, p13, p23
  real    :: a11, a22, a33, a12, a13, a23, a_pq_s_pq
  real    :: e11, e22, e33, e12, e13, e23
  real    :: s11, s22, s33, s12, s13, s23
  real    :: v12, v13, v23
  real    :: eps_h_11, eps_h_22, eps_h_33, eps_h_12, eps_h_13, eps_h_23
  real    :: eps_2_kin
  real    :: f_s,c_1,c_2
  real    :: phi_2_IP_11, phi_2_IP_22, phi_2_IP_33
  real    :: phi_2_IP_12, phi_2_IP_13, phi_2_IP_23, phi_km_2_IP_n_k_n_m
  real    :: phi_ij_w
  real    :: phi_ij_2_w
  real    :: u_k_u_m_n_k_n_m, phi_ij_1_second_term
  real    :: prod_and_coriolis, phi_ij, phi_ij_2, eps_h, stress
  real    :: a_, a_2, a_3
  real    :: e_, e_2, e_3
  real    :: c_, c_1_w, c_2_w, c_1_prime, c_3, c_4, c_5, f_, f_w
  real    :: n1n1, n2n2, n3n3, n1n2, n1n3, n2n3
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

  dissipation_2_term = 0.

  call Time_And_Length_Scale(grid) ! duplicated in Main_Pro

  !--------------------------------------------!
  !   preparations before filling source term  !
  !--------------------------------------------!

  do c = 1, grid % n_cells
    eps_lim(c) = max(eps % n(c), TINY) ! limited eps % n
    kin_lim(c) = max(kin % n(c), TINY) ! limited kin % n
    ! page 142 Re_t
    re_t(c)  = kin % n(c) * t_scale(c) / viscosity
  end do

  if(name_phi .ne. 'EPS') then ! it is not eps eq.

    ! needed for n_i_n_j
    call Grad_Mod_For_Phi(grid, l_scale, 1, l_sc_x,.true.)
    call Grad_Mod_For_Phi(grid, l_scale, 2, l_sc_y,.true.)
    call Grad_Mod_For_Phi(grid, l_scale, 3, l_sc_z,.true.)

  else ! it is eps eq.

    ! page 165 for term at -2*\nu(...)
    do iterator = 1, 3
      if(iterator == 1) then
        call Grad_Mod_For_Phi(grid, u % x, 1, u_xx, .true.) ! d2U/dxdx
        call Grad_Mod_For_Phi(grid, u % y, 2, u_yy, .true.) ! d2U/dydy
        call Grad_Mod_For_Phi(grid, u % z, 3, u_zz, .true.) ! d2U/dzdz
        call Grad_Mod_For_Phi(grid, u % x, 2, u_xy, .true.) ! d2U/dxdy
        call Grad_Mod_For_Phi(grid, u % x, 3, u_xz, .true.) ! d2U/dxdz
        call Grad_Mod_For_Phi(grid, u % y, 3, u_yz, .true.) ! d2U/dydz
      end if
      if(iterator == 2) then
        call Grad_Mod_For_Phi(grid, v % x, 1, u_xx, .true.) ! d2V/dxdx
        call Grad_Mod_For_Phi(grid, v % y, 2, u_yy, .true.) ! d2V/dydy
        call Grad_Mod_For_Phi(grid, v % z, 3, u_zz, .true.) ! d2V/dzdz
        call Grad_Mod_For_Phi(grid, v % x, 2, u_xy, .true.) ! d2V/dxdy
        call Grad_Mod_For_Phi(grid, v % x, 3, u_xz, .true.) ! d2V/dxdz
        call Grad_Mod_For_Phi(grid, v % y, 3, u_yz, .true.) ! d2V/dydz
      end if
      if(iterator == 3) then
        call Grad_Mod_For_Phi(grid, w % x, 1, u_xx, .true.) ! d2W/dxdx
        call Grad_Mod_For_Phi(grid, w % y, 2, u_yy, .true.) ! d2W/dydy
        call Grad_Mod_For_Phi(grid, w % z, 3, u_zz, .true.) ! d2W/dzdz
        call Grad_Mod_For_Phi(grid, w % x, 2, u_xy, .true.) ! d2W/dxdy
        call Grad_Mod_For_Phi(grid, w % x, 3, u_xz, .true.) ! d2W/dxdz
        call Grad_Mod_For_Phi(grid, w % y, 3, u_yz, .true.) ! d2W/dydz
      end if

      do  c = 1, grid % n_cells
        dissipation_2_term(c) = dissipation_2_term(c)                       +  &
          uu % n(c)*( u_xx(c)*u_xx(c) + u_xy(c)*u_xy(c) + u_xz(c)*u_xz(c) ) +  &
          uv % n(c)*( u_xx(c)*u_xy(c) + u_xy(c)*u_yy(c) + u_xz(c)*u_yz(c) ) +  &
          uw % n(c)*( u_xx(c)*u_xz(c) + u_xy(c)*u_yz(c) + u_xz(c)*u_zz(c) ) +  &
          uv % n(c)*( u_xy(c)*u_xx(c) + u_yy(c)*u_xy(c) + u_yz(c)*u_xz(c) ) +  &
          vv % n(c)*( u_xy(c)*u_xy(c) + u_yy(c)*u_yy(c) + u_yz(c)*u_yz(c) ) +  &
          vw % n(c)*( u_xy(c)*u_xz(c) + u_yy(c)*u_yz(c) + u_yz(c)*u_zz(c) ) +  &
          uw % n(c)*( u_xz(c)*u_xx(c) + u_yz(c)*u_xy(c) + u_zz(c)*u_xz(c) ) +  &
          vw % n(c)*( u_xz(c)*u_xy(c) + u_yz(c)*u_yy(c) + u_zz(c)*u_yz(c) ) +  &
          ww % n(c)*( u_xz(c)*u_xz(c) + u_yz(c)*u_yz(c) + u_zz(c)*u_zz(c) )
      end do ! c
    end do ! i

    dissipation_2_term(:) = c_3e*viscosity * t_scale(c) * dissipation_2_term(:)

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

    ! formula 2.4 | page 604
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

    !--------------------------------------------------!
    !   iterative procedure to find f_s and eps_h_ij   !
    !--------------------------------------------------!
    ! initial value for e_ = a_
    e_ = a_
    f_s = 1. - sqrt(a_) * e_**2.

    do iterator = 1, 6
      ! formula 2.14
      eps_h_11 = f_s * uu % n(c)*t_scale(c) + (1.-f_s) * TWO_THIRDS*eps % n(c)
      eps_h_22 = f_s * vv % n(c)*t_scale(c) + (1.-f_s) * TWO_THIRDS*eps % n(c)
      eps_h_33 = f_s * ww % n(c)*t_scale(c) + (1.-f_s) * TWO_THIRDS*eps % n(c)
      eps_h_12 = f_s * uv % n(c)*t_scale(c)
      eps_h_13 = f_s * uw % n(c)*t_scale(c)
      eps_h_23 = f_s * vw % n(c)*t_scale(c)

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

      ! for formula 10 Phi_ij^w
      u_k_u_m_n_k_n_m = uu % n(c)*n1n1 +    vv % n(c)*n2n2 +    ww % n(c)*n3n3 &
                   + 2.*uv % n(c)*n1n2 + 2.*uw % n(c)*n1n3 + 2.*vw % n(c)*n2n3

      ! formula 12: f
      f_  = min((re_t(c)/150)**1.5, 1.)
      ! formula 12: C
      c_  = 2.5*a_*min(0.6, a_2)**0.25*f_
      ! formula 12: C_1
      c_1  = c_ + sqrt(a_)*e_**2.
      ! formula 12: C_1_prime
      c_1_prime  = max(0.8*a_2,0.5)*c_1
      ! formula 12: C_2
      c_2  = 0.6*sqrt(a_)
      ! formula 12: C_3
      c_3  = (4./3.)*c_2
      ! formula 14
      c_1_w  = max(0.9 - 0.7*c_, 0.3)
      ! formula 13
      c_2_w  = min(a_,0.3)
      ! formula 14
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

      s11 = u % x(c) 
      s22 = v % y(c) 
      s33 = w % z(c) 
      s12 = 0.5*(u % y(c) + v % x(c))
      s13 = 0.5*(u % z(c) + w % x(c))
      s23 = 0.5*(v % z(c) + w % y(c))

      ! formula C.6
      v12 = 0.5*(u % y(c) - v % x(c))
      !v21 = -v12
      v13 = 0.5*(u % z(c) - w % x(c))
      !v31 = -v13
      v23 = 0.5*(v % z(c) - w % y(c))
      !v32 = -v23

      ! formula 11: Phi_ij,2^IP
      phi_2_IP_11 = - c_2 * (p11 - TWO_THIRDS*p_kin(c) )
      phi_2_IP_22 = - c_2 * (p22 - TWO_THIRDS*p_kin(c) )
      phi_2_IP_33 = - c_2 * (p33 - TWO_THIRDS*p_kin(c) )
      phi_2_IP_12 = - c_2 *  p12
      phi_2_IP_13 = - c_2 *  p13
      phi_2_IP_23 = - c_2 *  p23

      if(name_phi .eq. 'UU' .or. &
         name_phi .eq. 'VV' .or. &
         name_phi .eq. 'WW') then

        ! for formula 9
        a_pq_s_pq = a11*s11 + a22*s22 + a33*s33 &
             + 2.*(a12*s12 + a13*s13 + a23*s23)

        ! for formula 10: Phi_km,2^IP * n_k * n_m
        phi_km_2_IP_n_k_n_m =    phi_2_IP_11*n1n1  &
                            +    phi_2_IP_33*n3n3  &
                            +    phi_2_IP_22*n2n2  &
                            + 2.*phi_2_IP_12*n1n2  &
                            + 2.*phi_2_IP_13*n1n3  &
                            + 2.*phi_2_IP_23*n2n3
      end if

      !---------------!
      !   uu stress   !
      !---------------!
      if(name_phi .eq. 'UU') then
        ! limited stress
        stress = uu % n(c)

      ! formula 9: Phi_ij,1, second term
        phi_ij_1_second_term = c_1_prime * eps % n(c) * &
          ( a11*a11 + a12*a12 + a13*a13 - ONE_THIRD*a_2 )
        
        ! formula 9: Phi_ij,2
        phi_ij_2 = kin % n(c)*( c_3*s11                                     &
          + c_4*( 2.*( a11*s11 + a12*s12 + a13*s13) - TWO_THIRDS*a_pq_s_pq )&
          - c_5*  2.*(           a12*v12 + a13*v13)                         )

        ! formula 10, Phi_ij^w
        phi_ij_w = f_w * eps_2_kin * (                                    &
          c_1_w * ( u_k_u_m_n_k_n_m                                       &
          -3.*(uu % n(c)  *n1n1 + uv % n(c)  *n1n2 + uw % n(c)  *n1n3  ) )&
          + c_2_w * ( phi_km_2_IP_n_k_n_m                                 &
          -3.*(phi_2_IP_11*n1n1 + phi_2_IP_12*n1n2 + phi_2_IP_13*n1n3) )  )

        prod_and_coriolis = p11

        eps_h = eps_h_11
      end if
      !---------------!
      !   vv stress   !
      !---------------!
      if(name_phi .eq. 'VV') then
        ! limited stress
        stress = vv % n(c)

        prod_and_coriolis = p22

        ! formula 9: Phi_ij,1, second term
        phi_ij_1_second_term = c_1_prime * eps % n(c) * &
          ( a12*a12 + a22*a22 + a23*a23 - ONE_THIRD*a_2 )

        ! formula 9: Phi_ij,2
        phi_ij_2 = kin % n(c)*( c_3*s22                                     &
          + c_4*( 2.*( a12*s12 + a22*s22 + a23*s23) - TWO_THIRDS*a_pq_s_pq )&
          - c_5*  2.*(-a12*v12           + a23*v23)                         )

        ! formula 10, Phi_ij^w
        phi_ij_w = f_w * eps_2_kin * (                                    &
            c_1_w * ( u_k_u_m_n_k_n_m                                     &
          -3.*(uv % n(c)  *n1n2 + vv % n(c)  *n2n2 + vw % n(c)  *n2n3  ) )&
          + c_2_w * ( phi_km_2_IP_n_k_n_m                                 &
          -3.*(phi_2_IP_12*n1n2 + phi_2_IP_22*n2n2 + phi_2_IP_23*n2n3) )  )

        eps_h = eps_h_22
      end if
      !---------------!
      !   ww stress   !
      !---------------!
      if(name_phi .eq. 'WW') then
        ! limited stress
        stress = ww % n(c)

        prod_and_coriolis = p33

        ! formula 9: Phi_ij,1, second term
        phi_ij_1_second_term = c_1_prime * eps % n(c) * &
          ( a13*a13 + a23*a23 + a33*a33 - ONE_THIRD*a_2 )
        
        ! formula 9: Phi_ij,2
        phi_ij_2 = kin % n(c)*( c_3*s33                                     &
          + c_4*( 2.*( a13*s13 + a23*s23 + a33*s33) - TWO_THIRDS*a_pq_s_pq )&
          - c_5*  2.*(-a13*v13 - a23*v23          )                         )

        ! formula 10, Phi_ij^w
        phi_ij_w = f_w * eps_2_kin * (                                    &
            c_1_w * ( u_k_u_m_n_k_n_m                                     &
          -3.*(uw % n(c)  *n1n3 + vw % n(c)  *n2n3 + ww % n(c)  *n3n3  ) )&
          + c_2_w * ( phi_km_2_IP_n_k_n_m                                 &
          -3.*(phi_2_IP_13*n1n3 + phi_2_IP_23*n2n3 + phi_2_IP_33*n3n3) )  )

        eps_h = eps_h_33
      end if
      !---------------!
      !   uv stress   !
      !---------------!
      if(name_phi .eq. 'UV') then
        ! limited stress
        stress = uv % n(c)

        prod_and_coriolis = p12

        ! formula 9: Phi_ij,1, second term
        phi_ij_1_second_term = c_1_prime * eps % n(c) * &
          ( a11*a12 + a12*a22 + a13*a23 )

        ! formula 9: Phi_ij,2
        phi_ij_2 = kin % n(c)*( c_3*s12                                        &
          + c_4*( a11*s12 + a12*s22 + a13*s23 + a12*s11 + a22*s12 + a23*s13 )  &
          - c_5*(-a11*v12           + a13*v23 +           a22*v12 + a23*v13 )  )

        ! formula 10, Phi_ij^w
        phi_ij_w = - 1.5 * f_w * eps_2_kin * (                                 &
            c_1_w * ( uu % n(c)  *n1n2 + uv % n(c)  *n2n2 + uw % n(c)  *n2n3 + &
                      uv % n(c)  *n1n1 + vv % n(c)  *n1n2 + vw % n(c)  *n1n3 ) &
          + c_2_w * ( phi_2_IP_11*n1n2 + phi_2_IP_12*n2n2 + phi_2_IP_13*n2n3 + &
                      phi_2_IP_12*n1n1 + phi_2_IP_22*n1n2 + phi_2_IP_23*n1n3 ) )

        eps_h = eps_h_12

      end if
      !---------------!
      !   uw stress   !
      !---------------!
      if(name_phi .eq. 'UW') then
        ! limited stress
        stress = uw % n(c)

        prod_and_coriolis = p13

        ! formula 9: Phi_ij,1, second term
        phi_ij_1_second_term = c_1_prime * eps % n(c) * &
          ( a11*a13 + a12*a23 + a13*a33 )

        ! formula 9: Phi_ij,2
        phi_ij_2 = kin % n(c)*( c_3*s13                                        &
          + c_4*( a11*s13 + a12*s23 + a13*s33 + a13*s11 + a23*s12 + a33*s13 )  &
          - c_5*(-a11*v13 - a12*v23 +                     a23*v12 + a33*v13 )  )

        ! formula 10, Phi_ij^w
        phi_ij_w = - 1.5 * f_w * eps_2_kin * (                                 &
            c_1_w * ( uu % n(c)  *n1n3 + uv % n(c)  *n2n3 + uw % n(c)  *n3n3 + &
                      uw % n(c)  *n1n1 + vw % n(c)  *n1n2 + ww % n(c)  *n1n3 ) &
          + c_2_w * ( phi_2_IP_11*n1n3 + phi_2_IP_12*n2n3 + phi_2_IP_13*n3n3 + &
                      phi_2_IP_13*n1n1 + phi_2_IP_23*n1n2 + phi_2_IP_33*n1n3)  )

        eps_h = eps_h_13

      end if
      !---------------!
      !   vw stress   !
      !---------------!
      if(name_phi .eq. 'VW') then
        ! limited stress
        stress = vw % n(c)

        prod_and_coriolis = p23

        ! formula 9: Phi_ij,1, second term
        phi_ij_1_second_term = c_1_prime * eps % n(c) * &
          ( a12*a13 + a22*a23 + a23*a33 )

        ! formula 9: Phi_ij,2
        phi_ij_2 = kin % n(c)*( c_3*s23                                        &
          + c_4*( a12*s13 + a22*s23 + a23*s33 + a13*s12 + a23*s22 + a33*s23 )  &
          - c_5*(-a12*v13 - a22*v23            -a13*v12 +           a33*v23 )  )

        ! formula 10, Phi_ij^w
        phi_ij_w = - 1.5 * f_w * eps_2_kin * (                                 &
            c_1_w * ( uv % n(c)  *n1n3 + vv % n(c)  *n2n3 + vw % n(c)  *n3n3 + &
                      uw % n(c)  *n1n2 + vw % n(c)  *n2n2 + ww % n(c)  *n2n3 ) &
          + c_2_w * ( phi_2_IP_12*n1n3 + phi_2_IP_22*n2n3 + phi_2_IP_23*n3n3 + &
                      phi_2_IP_13*n1n2 + phi_2_IP_23*n2n2 + phi_2_IP_33*n2n3 ) )

        eps_h = eps_h_23

      end if

      !-------------------------------------!
      !   repeating part for all stresses   !
      !-------------------------------------!
      ! limit stress
      if (stress .lt. 0) stress = min(stress,-TINY)
      if (stress .ge. 0) stress = max(stress, TINY)

      phi_ij   = phi_ij_1_second_term + phi_ij_2 + phi_ij_w

      ! term depending on Kronecker symbol 
      if(name_phi .eq. 'UU' .or. &
         name_phi .eq. 'VV' .or. &
         name_phi .eq. 'WW') then
        b(c) = b(c) + grid % vol(c) *  (   &
          c_1 * eps % n(c)*TWO_THIRDS      & ! first term in Phi_ij_1
          - (1.-f_s)*TWO_THIRDS*eps % n(c) ) ! part of eps_ij^h
      end if

      b(c) = b(c) + grid % vol(c) *  (     &
        + max(prod_and_coriolis, 0.)       & ! P_ij + G_ij, if > 0
        + max(phi_ij, 0.)                  ) ! Phi_ij     , if > 0

      A % val(A % dia(c)) = A % val(A % dia(c)) + grid % vol(c) * ( &
          c_1 * eps_2_kin                  & ! first term in Phi_ij_1
        + f_s * t_scale(c)                 & ! part of eps_ij^h / u_iu_j
        + (                                &
        - min(prod_and_coriolis, 0.)       & ! (P_ij + G_ij) / u_iu_j, if < 0
        - min(phi_ij, 0.)                  & ! (Phi_ij) / u_iu_j, if < 0
          ) / stress                       )
  !----------------------!
  !   Epsilon equation   !
  !----------------------!
    else ! it is eps eq.

      ! page 165, diss. eq., second term
      b(c) = b(c) + grid % vol(c) * dissipation_2_term(c)

      A % val(A % dia(c)) = A % val(A % dia(c)) + grid % vol(c) * ( &
          c_1e * p_kin(c) / kin_lim(c)   & ! page 165, diss. eq., first term
        - c_2e * eps % n(c) / kin_lim(c) ) ! diss. eq., third last
    end if
  end do

  end subroutine