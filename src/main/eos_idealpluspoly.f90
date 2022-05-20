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

  public :: get_idealpluspoly_temp,&
    get_idealpluspoly_press_over_rho,&
    get_idealpluspoly_spsoundi,&
    get_idealpluspoly_temp_from_pres,&
    get_idealpluspoly_en_from_temp

  private

  interface get_idealpluspoly_temp
    module procedure &
      get_idealpluspoly_temp_from_full_en,&
      get_idealpluspoly_temp_from_therm_en
  end interface

  interface get_idealpluspoly_press_over_rho
    module procedure &
      ! get_idealpluspoly_full_press_over_rho_from_full_en,&
      get_idealpluspoly_full_press_over_rho_from_therm_en,&
      get_idealpluspoly_full_press_over_rho_from_temp,&
      get_idealpluspoly_therm_press_over_rho_from_therm_en
  end interface

  interface get_idealpluspoly_spsoundi
    module procedure &
      ! get_idealpluspoly_spsoundi_from_full_en,&
      get_idealpluspoly_spsoundi_from_therm_en,&
      get_idealpluspoly_spsoundi_from_temp
  end interface

  interface get_idealpluspoly_en_from_temp
  module procedure &
    get_idealpluspoly_full_en_from_temp,&
    get_idealpluspoly_therm_en_from_temp
  end interface

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
  subroutine get_idealpluspoly_temp_from_full_en(rhoi,eni,polyk,gamma,mu,tempi)

    use units,     only: unit_ergg

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg]
    real, intent(in)  :: polyk, gamma, mu

    real, intent(out) :: tempi ! [K]

    real              :: u_poly ! [unit_ergg]

    u_poly = polyk/(gamma-1.)*rhoi**(gamma-1.)
    if (eni < u_poly) then
      tempi = 0.
      ! write(*,*) 'get_idealpluspoly_temp: T < 0'
      ! write(*,'(2(a,es13.6))')&
      !   'eni = ', eni,&
      !   ',  u_poly = ', u_poly
      return
    endif

    tempi = (2.0/3.0*(eni - u_poly)*mu/Rg)*unit_ergg

  end subroutine get_idealpluspoly_temp_from_full_en
!----------------------------------------------------------------
  subroutine get_idealpluspoly_temp_from_therm_en(eni,mu,tempi)

    use units,     only: unit_ergg

    real, intent(in)  :: eni ! [unit_ergg] -- here only therm part of energy
    real, intent(in)  :: mu

    real, intent(out) :: tempi ! [K]

    tempi = (2.0/3.0*eni*mu/Rg)*unit_ergg

  end subroutine get_idealpluspoly_temp_from_therm_en
!----------------------------------------------------------------
!  u_{\mathrm{poly}} = K \frac{\rho^{(\gamma - 1)}}{(\gamma - 1)}
!  % [unit_{ergg}]
!  \\
!  \frac{P}{\rho} = (\gamma - 1) u_{\mathrm{poly}} + (\tilde \gamma - 1) (u - u_{\mathrm{poly}})
!  \\
!  \tilde \gamma = 5/3
!  \\
!  \frac{P}{\rho} = (\gamma - \frac{5}{3}) u_{\mathrm{poly}} + \frac{2}{3} u
!  \frac{P}{\rho} = (\gamma - 1) u_{\mathrm{poly}} + \frac{2}{3} u_{\mathrm{therm}}
  subroutine get_idealpluspoly_full_press_over_rho_from_full_en(rhoi,eni,polyk,gamma,ponrhoi)

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg]
    real, intent(in)  :: polyk, gamma

    real, intent(out) :: ponrhoi ! [unit_ergg]

    real              :: u_poly ! [unit_ergg]

    u_poly = polyk/(gamma-1.)*rhoi**(gamma-1.)
    ponrhoi = (gamma-5.0/3.0)*u_poly + 2.0/3.0*eni

  end subroutine get_idealpluspoly_full_press_over_rho_from_full_en
!----------------------------------------------------------------
  subroutine get_idealpluspoly_full_press_over_rho_from_therm_en(rhoi,eni,polyk,gamma,ponrhoi)

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg] -- here only therm part of energy
    real, intent(in)  :: polyk, gamma

    real, intent(out) :: ponrhoi ! [unit_ergg]

    ponrhoi = polyk*rhoi**(gamma-1.) + 2.0/3.0*eni

  end subroutine get_idealpluspoly_full_press_over_rho_from_therm_en
!----------------------------------------------------------------
  subroutine get_idealpluspoly_full_press_over_rho_from_temp(rhoi,polyk,gamma,mu,tempi,ponrhoi)

    use units,     only: unit_ergg

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: polyk, gamma, mu
    real, intent(in)  :: tempi ! [K]

    real, intent(out) :: ponrhoi ! [unit_ergg]

    ponrhoi = polyk*rhoi**(gamma-1.) + Rg*tempi/mu/unit_ergg

  end subroutine get_idealpluspoly_full_press_over_rho_from_temp
!----------------------------------------------------------------
  subroutine get_idealpluspoly_therm_press_over_rho_from_therm_en(eni,ponrhoi)

    real, intent(in)  :: eni ! [unit_ergg] -- here only therm part of energy

    real, intent(out) :: ponrhoi ! [unit_ergg]

    ponrhoi = 2.0/3.0*eni

  end subroutine get_idealpluspoly_therm_press_over_rho_from_therm_en
!----------------------------------------------------------------
!  % Pressure
!  P = K \rho^{\gamma} + \frac{\rho R T}{\mu}
!  \\ % Energy
!  u = K \frac{\rho^{\gamma - 1}}{\gamma - 1} + \frac{3}{2}\frac{R T}{\mu}
!  \\ % Entropy -- Sackur--Tetrode equation
!  S \sim \ln\left(\frac{T^{\;\frac{3}{2}}}{\rho}\right)
!  \\ % Speed of sound
!  c^2_S = \left(\frac{\partial P}{\partial \rho}\right)_S
!    = \frac{\partial(P,S)}{\partial(\rho,S)}
!    =
!    \\
!    = \left(
!        \left(\frac{\partial P}{\partial \rho}\right)_T
!        \left(\frac{\partial S}{\partial T}\right)_\rho
!        -
!        \left(\frac{\partial P}{\partial T}\right)_\rho
!        \left(\frac{\partial S}{\partial \rho}\right)_T
!        \right)
!    \\
!      \left(
!        \left(\frac{\partial \rho}{\partial \rho}\right)_T
!        \left(\frac{\partial S}{\partial T}\right)_\rho
!        -
!        \left(\frac{\partial \rho}{\partial T}\right)_\rho
!        \left(\frac{\partial S}{\partial \rho}\right)_T
!        \right)^{-1}
!    =
!    \\
!    = \gamma K \rho^{\gamma - 1} + \frac{5}{3}\frac{R T}{\mu}
!    =
!    \\
!    = \gamma(\gamma - 1) u_{\mathrm{poly}} + \frac{10}{9} u_{\mathrm{therm}}
!    =
!    \\
!    = \frac{10}{9} u
!      + K \rho^{\gamma - 1} \left(\gamma - \frac{10}{9 (\gamma - 1)}\right)
  subroutine get_idealpluspoly_spsoundi_from_full_en(rhoi,eni,polyk,gamma,spsoundi)

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg]
    real, intent(in)  :: polyk, gamma

    real, intent(out) :: spsoundi ! [unit_velocity]

    spsoundi = 10.0/9.0*eni + polyk*(rhoi**(gamma-1.))*(gamma - 10.0/(9.0*(gamma - 1.)))
    spsoundi = sqrt(spsoundi)

  end subroutine get_idealpluspoly_spsoundi_from_full_en
!----------------------------------------------------------------
  subroutine get_idealpluspoly_spsoundi_from_therm_en(rhoi,eni,polyk,gamma,spsoundi)

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: eni ! [unit_ergg]
    real, intent(in)  :: polyk, gamma

    real, intent(out) :: spsoundi ! [unit_velocity]

    spsoundi = gamma*(gamma-1.)*polyk*(rhoi**(gamma-1.)) + 10.0/9.0*eni
    spsoundi = sqrt(spsoundi)

  end subroutine get_idealpluspoly_spsoundi_from_therm_en
!----------------------------------------------------------------
  subroutine get_idealpluspoly_spsoundi_from_temp(rhoi,polyk,gamma,mu,tempi,spsoundi)

    use units,     only: unit_ergg

    real, intent(in)  :: rhoi ! [unit_density]
    real, intent(in)  :: polyk, gamma, mu
    real, intent(in)  :: tempi ! [K]

    real, intent(out) :: spsoundi ! [unit_velocity]

    spsoundi = gamma*polyk*rhoi**(gamma-1.) + 5.0/3.0*Rg*tempi/mu/unit_ergg
    spsoundi = sqrt(spsoundi)

  end subroutine get_idealpluspoly_spsoundi_from_temp
!----------------------------------------------------------------
!+
!  Calculates temperature from pressure and density
!+
!----------------------------------------------------------------
  subroutine get_idealpluspoly_temp_from_pres(presi_cgs,rhoi_cgs,polyk,gamma,mu,tempi)

    use units,     only: unit_density, unit_pressure

    real, intent(in)  :: presi_cgs ! [g/(cm*s**2)]
    real, intent(in)  :: rhoi_cgs ! [g/cm**3]
    real, intent(in)  :: polyk, gamma, mu

    real, intent(out) :: tempi ! [K]

    real              :: rhoi ! [unit_density]
    real              :: presi_poly_cgs ! [g/(cm*s**2)]

    rhoi = rhoi_cgs/unit_density
    presi_poly_cgs = (polyk*rhoi**gamma)*unit_pressure
    if (presi_cgs < presi_poly_cgs) then
      tempi = 0.
      write(*,*) 'get_idealpluspoly_temp_from_pres: T < 0'
      write(*,'(2(a,3x,es13.6))')&
        'presi_cgs = ', presi_cgs,&
        ',  presi_poly_cgs = ', presi_poly_cgs
      return
    endif

    tempi = presi_cgs - presi_poly_cgs
    tempi = mu*tempi/(rhoi_cgs*Rg)

    ! NB:
    tempi = 0.0

  end subroutine get_idealpluspoly_temp_from_pres
!----------------------------------------------------------------
!+
!  Calculates internal energy per unit mass from density
!  and temperature
!+
!----------------------------------------------------------------
  subroutine get_idealpluspoly_full_en_from_temp(rhoi_cgs,polyk,gamma,mu,tempi,eni_cgs)

    use units,     only: unit_density, unit_ergg

    real, intent(in)  :: rhoi_cgs ! [g/cm**3]
    real, intent(in)  :: polyk, gamma, mu
    real, intent(in)  :: tempi ! [K]

    real, intent(out) :: eni_cgs ! [erg/g]

    real              :: rhoi ! [unit_density]

    real              :: u_poly_cgs ! [erg/g]
    real              :: u_therm_cgs ! [erg/g]

    rhoi = rhoi_cgs/unit_density
    u_poly_cgs = (polyk/(gamma-1.)*rhoi**(gamma-1.))*unit_ergg
    u_therm_cgs = 3.0/2.0*Rg*tempi/mu
    eni_cgs = u_poly_cgs + u_therm_cgs

  end subroutine get_idealpluspoly_full_en_from_temp
!----------------------------------------------------------------
  subroutine get_idealpluspoly_therm_en_from_temp(mu,tempi,eni_cgs)

    real, intent(in)  :: mu
    real, intent(in)  :: tempi ! [K]

    real, intent(out) :: eni_cgs ! [erg/g] -- here only therm part of energy

    eni_cgs = 3.0/2.0*Rg*tempi/mu

  end subroutine get_idealpluspoly_therm_en_from_temp

end module eos_idealpluspoly
