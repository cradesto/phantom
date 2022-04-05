!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module eos_idealpluspoly
!
! Ideal gas equation of state plus radiation pressure, assumes
!               inputs are in cgs units
!
! :References: Stellar Structure and Evolution (2nd Edition) (Kippenhahn,
!              Weigert, Weiss)
!
! :Owner: Potashov Marat
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
  use physcon,  only: Rg, radconst

  implicit none

  real, parameter :: tolerance = 1e-15

  public :: get_idealpluspoly_temp,&
            get_idealpluspoly_press_over_rho,&
            get_idealpluspoly_spsoundi,&
            get_idealpluspoly_temp_from_pres,&
            get_idealpluspoly_en_from_temp

  private

contains
!----------------------------------------------------------------
!  Solve for temperature as a function of (poly+gas) internal energy
!  per unit mass (eni [unit_ergg]) and density (rhoi [unit_density])
!  u_{\mathrm{poly}} = K \frac{\rho^{(\gamma - 1)}}{(\gamma - 1)}
!  % [unit_{ergg}]
!  \\
!  u_{\mathrm{therm}} = \frac{R T}{\mu (\tilde \gamma - 1)}
!  % [erg/g]
!  \\
!  u = u_{\mathrm{poly}} + u_{\mathrm{therm}}{unit_{ergg}}^{-1}
!  \\
!  T = (u - u_{\mathrm{poly}})\frac{\mu (\tilde \gamma - 1)}{R} unit_{ergg}
!  \\
!  \tilde \gamma = 5/3
!  \\
!  T = \frac{2}{3}(u - u_{\mathrm{poly}})\frac{\mu}{R} unit_{ergg}
!----------------------------------------------------------------
  subroutine get_idealpluspoly_temp(rhoi,eni,mu,polyk,gamma,tempi)

    use units,     only: unit_ergg

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg]
    real, intent(in)  :: mu, polyk, gamma

    real, intent(out) :: tempi ! [K]

    real              :: u_poly ! [unit_ergg]

    u_poly = polyk/(gamma-1.)*rhoi**(gamma-1.)
    tempi = (2.0/3.0*(eni - u_poly)*mu/Rg)*unit_ergg

  end subroutine get_idealpluspoly_temp
!----------------------------------------------------------------
!  u_{\mathrm{poly}} = K \frac{\rho^{(\gamma - 1)}}{(\gamma - 1)}
!  % [unit_{ergg}]
!  \\
!  \frac{P}{\rho} = (\gamma - 1) u_{\mathrm{poly}} + (\tilde \gamma - 1) (u - u_{\mathrm{poly}})
!  \\
!  \tilde \gamma = 5/3
!  \\
!  \frac{P}{\rho} = (\gamma - \frac{5}{3}) u_{\mathrm{poly}} + \frac{2}{3} u
  subroutine get_idealpluspoly_press_over_rho(rhoi,eni,polyk,gamma,ponrhoi)

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg]
    real, intent(in)  :: polyk, gamma

    real, intent(out) :: ponrhoi

    real              :: u_poly ! [unit_ergg]

    u_poly = polyk/(gamma-1.)*rhoi**(gamma-1.)
    ponrhoi = (gamma-5.0/3.0)*u_poly + 2.0/3.0*eni

  end subroutine get_idealpluspoly_press_over_rho
!----------------------------------------------------------------
  subroutine get_idealpluspoly_spsoundi(rhoi,presi,eni,spsoundi)

    real, intent(in)  :: rhoi,presi,eni

    real, intent(out) :: spsoundi

    real              :: gamma

    gamma = 1. + presi/(eni*rhoi)
    spsoundi = sqrt(gamma*presi/rhoi)

  end subroutine get_idealpluspoly_spsoundi
!----------------------------------------------------------------
!+
!  Calculates temperature from pressure and density
!+
!----------------------------------------------------------------
  subroutine get_idealpluspoly_temp_from_pres(presi,rhoi,mu,tempi)
    real, intent(in)    :: rhoi,presi,mu
    real, intent(inout) :: tempi
    real                :: imu,numerator,denominator,correction,temp_new
    integer             :: iter
    integer, parameter  :: iter_max = 1000

    iter = 0
    correction = huge(0.)
    imu = 1./mu
    tempi = min((3.*presi/radconst)**0.25, presi*mu/(rhoi*Rg))
    do while (abs(correction) > tolerance*tempi .and. iter < iter_max)
      numerator   = presi - rhoi*Rg*tempi*imu - radconst*tempi**4 /3.
      denominator =  - rhoi*Rg*imu - 4./3.*radconst*tempi**3
      correction  = numerator/denominator
      temp_new = tempi - correction
      if (temp_new > 1.2 * tempi) then
        tempi = 1.2 * tempi
      elseif (temp_new < 0.8 * tempi) then
        tempi = 0.8 * tempi
      else
        tempi = temp_new
      endif
      iter = iter + 1
    enddo

  end subroutine get_idealpluspoly_temp_from_pres
!----------------------------------------------------------------
!+
!  Calculates internal energy per unit mass from density
!  and temperature
!+
!----------------------------------------------------------------
  subroutine get_idealpluspoly_en_from_temp(densi,tempi,mu,gamma,eni)
    real, intent(in)  :: densi,tempi,mu,gamma
    real, intent(out) :: eni

    eni = Rg*tempi/((gamma-1.)*mu) + radconst*tempi**4/densi

  end subroutine get_idealpluspoly_en_from_temp

end module eos_idealpluspoly
