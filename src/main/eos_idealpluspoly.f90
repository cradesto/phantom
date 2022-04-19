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

    real, intent(out) :: ponrhoi ! [unit_ergg]

    real              :: u_poly ! [unit_ergg]

    u_poly = polyk/(gamma-1.)*rhoi**(gamma-1.)
    ponrhoi = (gamma-5.0/3.0)*u_poly + 2.0/3.0*eni

  end subroutine get_idealpluspoly_press_over_rho
!----------------------------------------------------------------
! % Pressure
! P = K \rho^{\gamma} + \frac{\rho R T}{\mu}
! \\ % Energy
! u = K \frac{\rho^{\gamma - 1}}{\gamma - 1} + \frac{3}{2}\frac{R T}{\mu}
! \\ % Entropy -- Sackur--Tetrode equation
! S \sim \ln\left(\frac{T^{\;\frac{3}{2}}}{\rho}\right)
! \\ % Speed of sound
! c^2_S = \left(\frac{\partial P}{\partial \rho}\right)_S
!   = \frac{\partial(P,S)}{\partial(\rho,S)}
!   =
!   \\
!   = \left(
!       \left(\frac{\partial P}{\partial \rho}\right)_T
!       \left(\frac{\partial S}{\partial T}\right)_\rho
!       -
!       \left(\frac{\partial P}{\partial T}\right)_\rho
!       \left(\frac{\partial S}{\partial \rho}\right)_T
!       \right)
!   \\
!     \left(
!       \left(\frac{\partial \rho}{\partial \rho}\right)_T
!       \left(\frac{\partial S}{\partial T}\right)_\rho
!       -
!       \left(\frac{\partial \rho}{\partial T}\right)_\rho
!       \left(\frac{\partial S}{\partial \rho}\right)_T
!       \right)^{-1}
!   =
!   \\
!   = \gamma K \rho^{\gamma - 1} + \frac{5}{3}\frac{R T}{\mu}
!   =
!   \\
!   = \frac{10}{9} u
!     + K \rho^{\gamma - 1} \left(\gamma - \frac{10}{9 (\gamma - 1)}\right)
  subroutine get_idealpluspoly_spsoundi(rhoi,eni,polyk,gamma,spsoundi)

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg]
    real, intent(in)  :: polyk, gamma

    real, intent(out) :: spsoundi ! [unit_velocity]

    ! from temperature
    ! spsoundi = gamma*polyk*(rhoi**(gamma-1.)) + 5.0/3.0*Rg*tempi/mu
    spsoundi = 10.0/9.0*eni + polyk*(rhoi**(gamma-1.))*(gamma - 10.0/(9.0*(gamma - 1.)))
    spsoundi = sqrt(spsoundi)

  end subroutine get_idealpluspoly_spsoundi
!----------------------------------------------------------------
!+
!  Calculates temperature from pressure and density
!+
!----------------------------------------------------------------
  subroutine get_idealpluspoly_temp_from_pres(presi,rhoi,mu,polyk,gamma,tempi)

    real, intent(in)  :: presi ! [unit_pressure]
    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: mu, polyk, gamma

    real, intent(out) :: tempi ! [K]

    tempi = presi - polyk*rhoi**gamma
    tempi = mu*tempi/(rhoi*Rg)

  end subroutine get_idealpluspoly_temp_from_pres
!----------------------------------------------------------------
!+
!  Calculates internal energy per unit mass from density
!  and temperature
!+
!----------------------------------------------------------------
  subroutine get_idealpluspoly_en_from_temp(rhoi,tempi,mu,polyk,gamma,eni)

    use units,     only: unit_ergg

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: tempi ! [K]
    real, intent(in)  :: mu, polyk, gamma

    real, intent(out) :: eni ! [unit_ergg]

    real              :: u_poly ! [unit_ergg]
    real              :: u_therm ! [unit_ergg]

    u_poly = polyk/(gamma-1.)*rhoi**(gamma-1.)
    u_therm = 3.0/2.0*Rg*tempi/mu
    eni = u_poly + u_therm/unit_ergg

  end subroutine get_idealpluspoly_en_from_temp

end module eos_idealpluspoly
