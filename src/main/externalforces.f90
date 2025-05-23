!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module externalforces
!
! Routines dealing with external forces/ potentials
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - accradius1      : *soft accretion radius of central object*
!   - accradius1_hard : *hard accretion radius of central object*
!   - eps_soft        : *softening length (Plummer) for central potential in code units*
!   - mass1           : *mass of central object in code units*
!
! :Dependencies: dump_utils, extern_Bfield, extern_binary, extern_corotate,
!   extern_densprofile, extern_geopot, extern_gnewton, extern_gwinspiral,
!   extern_lensethirring, extern_prdrag, extern_spiral, extern_staticsine,
!   infile_utils, io, part, units
!
 use extern_binary,        only:accradius1,mass1,accretedmass1,accretedmass2
 use extern_corotate,      only:omega_corotate  ! so public from this module
 use extern_lensethirring, only:a=>blackhole_spin
 implicit none

 private
 public :: externalforce,externalforce_vdependent
 public :: accrete_particles,was_accreted
 public :: accradius1,omega_corotate,accretedmass1,accretedmass2
 public :: write_options_externalforces,read_options_externalforces
 public :: initialise_externalforces,is_velocity_dependent
 public :: update_vdependent_extforce
 public :: update_externalforce
 public :: write_headeropts_extern,read_headeropts_extern

 real, public :: eps_soft = 0.d0
 real, private :: eps2_soft = 0.d0
 real, public :: Rdisc = 5.

 real, public :: accradius1_hard = 0.
 logical, public :: extract_iextern_from_hdr = .false.

 public :: mass1,a

 !
 ! enumerated list of external forces
 !
 integer, parameter, public :: &
   iext_star          = 1, &
   iext_corotate      = 2, &
   iext_binary        = 3, &
   iext_prdrag        = 4, &
   iext_torus         = 5, &
   iext_toystar       = 6, &
   iext_externB       = 7, &
   iext_spiral        = 8, &
   iext_lensethirring = 9, &
   iext_densprofile   = 10, &
   iext_einsteinprec  = 11, &
   iext_gnewton       = 12, &
   iext_staticsine    = 13, &
   iext_gwinspiral    = 14, &
   iext_discgravity   = 15, &
   iext_corot_binary  = 16, &
   iext_geopot        = 17

 !
 ! Human-readable labels for these
 !
 integer, parameter, public  :: iexternalforce_max = 17
 character(len=*), parameter, public :: externalforcetype(iexternalforce_max) = (/ &
    'star                 ', &
    'corotate             ', &
    'binary               ', &
    'prdrag               ', &
    'torus                ', &
    'toystar              ', &
    'external B field     ', &
    'spiral               ', &
    'Lense-Thirring       ', &
    'density profile      ', &
    'Einstein-prec        ', &
    'generalised Newtonian', &
    'static sinusoid      ', &
    'grav. wave inspiral  ', &
    'disc gravity         ', &
    'corotating binary    ', &
    'geopotential model   '/)

contains
!-----------------------------------------------------------------------
!+
!  Computes external (body) forces on a particle given its co-ordinates
!+
!-----------------------------------------------------------------------
subroutine externalforce(iexternalforce,xi,yi,zi,hi,ti,fextxi,fextyi,fextzi,phi,dtf,ii)
 use extern_corotate,  only:get_centrifugal_force,get_companion_force,icompanion_grav
 use extern_binary,    only:binary_force
 use extern_prdrag,    only:get_prdrag_spatial_force
 use extern_gnewton,   only:get_gnewton_spatial_force
 use extern_spiral,    only:s_potential,schmidt_potential,&
  pichardo_potential,Wang_bar,LogDisc,&
  MNDisc,KFDiscSp,PlumBul,HernBul,HubbBul,COhalo,Flathalo,AMhalo,KBhalo,LMXbar,&
  LMTbar,Orthog_basisbar,DehnenBar,VogtSbar,BINReadPot3D,NFWhalo,&
  ibar,idisk,ihalo,ibulg,iarms,iread,Wadabar
 use extern_densprofile, only:densityprofile_force
 use extern_Bfield,      only:get_externalB_force
 use extern_staticsine,  only:staticsine_force
 use extern_gwinspiral,  only:get_gw_force_i
 use extern_geopot,      only:get_geopot_force,J2,spinvec
 use units,              only:get_G_code
 use io,                 only:fatal
 use part,               only:rhoh,massoftype,igas
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: xi,yi,zi,hi,ti
 real,    intent(out) :: fextxi,fextyi,fextzi,phi
 real,    intent(out), optional :: dtf
 integer, intent(in),  optional :: ii ! NOTE: index-base physics can be dangerous; treat with caution!
 real            :: r2,dr,dr3,r,d2,f2i
 real            :: rcyl2,rcyl,rsph,rsph3,v2onr,dtf1,dtf2
 real            :: phii,gcode,R_g,factor,rhoi
 real, parameter :: Rtorus = 1.0
 real,dimension(3) :: pos
!-----------------------------------------------------------------------
!
!--set external force to zero
!
 fextxi = 0.
 fextyi = 0.
 fextzi = 0.
 phi    = 0.

 select case(iexternalforce)

 case(iext_star,iext_lensethirring,iext_geopot)
!
!--1/r^2 force from central point mass
!
    r2 = xi*xi + yi*yi + zi*zi + eps2_soft

    if (r2 > epsilon(r2)) then
       dr = 1./sqrt(r2)
       dr3 = mass1*dr**3
       fextxi = fextxi - xi*dr3
       fextyi = fextyi - yi*dr3
       fextzi = fextzi - zi*dr3
       phi    = -mass1*dr
    endif

    if (iexternalforce==iext_geopot) then
       call get_geopot_force(xi,yi,zi,dr,dr3,accradius1,J2,spinvec,fextxi,fextyi,fextzi,phi)
    endif

 case(iext_corotate)
!
!--spatial part of forces in corotating frame, i.e. centrifugal force
!
    pos = (/xi,yi,zi/)
    call get_centrifugal_force(pos,fextxi,fextyi,fextzi,phi)
    if ( (icompanion_grav == 1) .or. (icompanion_grav == 2) ) then
       call get_companion_force(pos,fextxi,fextyi,fextzi,phi)
    endif

 case(iext_binary)
!
!--gravitational force from central binary
!
    call binary_force(xi,yi,zi,ti,fextxi,fextyi,fextzi,phi)

 case(iext_prdrag)
!
!--Spatial component of force due to luminous central object
!
    call get_prdrag_spatial_force(xi,yi,zi,mass1,fextxi,fextyi,fextzi,phi)

 case(iext_torus)
!
!--effective potential for equilibrium torus (Stone, Pringle, Begelman)
!  centripedal force balances pressure gradient and gravity of central point mass
!
    rcyl2 = xi*xi + yi*yi
    rcyl  = sqrt(rcyl2)
    rsph  = sqrt(rcyl2 + zi*zi)
    rsph3 = rsph**3
    v2onr = (Rtorus/(rcyl*rcyl2) - rcyl/rsph3)

    fextxi = v2onr*xi/rcyl
    fextyi = v2onr*yi/rcyl
    fextzi = -zi/rsph3   ! vertical component of gravity

    phi = - 1./rsph + 0.5*Rtorus/rcyl2

 case(iext_toystar)
!
!--Toy star force, centred on the origin
!
    r2 = xi**2 + yi**2 + zi**2
    if (r2 <= 0.25) then
       fextxi = -xi
       fextyi = -yi
       fextzi = -zi
       phi    = 0.5*r2
    endif

 case(iext_externB)
!
!--External force due to an assumed external B field (with non-zero curl)
!
    rhoi = rhoh(hi,massoftype(igas))
    call get_externalB_force(xi,yi,zi,hi,rhoi,fextxi,fextyi,fextzi)
    phi = 0.

 case(iext_spiral)

!--Spiral/galactic potential references in extern_spiral.
!--Calc r^2 (3D) and d^2 (2D), phi (x-y), theta (xy-z)
    gcode = get_G_code() !Needed to convert gg in codes natural units
    d2   = (xi*xi + yi*yi)       !or  r    = sqrt(d2+zi*zi)
    dr   = 1./sqrt(d2+zi*zi)    !1/r
    r    = 1./dr
    phii = atan2(yi,xi)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=DISCS=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case(idisk)
    case(0)
       !--No potential
    case(1)
       !--Binney&Tremmaine~Logdisc:
       call LogDisc(xi,yi,zi,d2,phi,fextxi,fextyi,fextzi)
    case(2)
       !--Miyamoto&Niagi~inc.Plummer/Kuzmin:
       call MNDisc(xi,yi,zi,d2,phi,fextxi,fextyi,fextzi)
    case(3)
       !--Khoperskov:Freeman disc+TWA spirals+bar:
       call KFDiscSp(xi,yi,zi,d2,r,phii,ti,phi,fextxi,fextyi,fextzi)
    end select
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=BULGES=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case(ibulg)
    case(0)
       !--No potential
    case(1)
       !--Plummer bulge:
       call PlumBul(xi,yi,zi,r,phi,fextxi,fextyi,fextzi)
    case(2)
       !--Hernequist bulge:
       call HernBul(xi,yi,zi,r,phi,fextxi,fextyi,fextzi)
    case(3)
       !--Hubble luminosity profile bulge:
       call HubbBul(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
    end select
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=HALO/CORONA-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case(ihalo)
    case(0)
       !--No potential
    case(1)
       !--Caldwell&Ostriker~ln/atan:
       call COhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
    case(2)
       !--Flat~ln:
       call flathalo(xi,yi,zi,r,phi,fextxi,fextyi,fextzi)
    case(3)
       !--Allen&Martos~gamma~power~law:
       call AMhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
    case(4)
       !--Khoperskov/Begeman~Isothermal:
       call KBhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
    case(5)
       !--Navarro/Frenk/White~CDMhalo:
       call NFWhalo(xi,yi,zi,r,dr,phi,fextxi,fextyi,fextzi)
    end select
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=SPIRALS-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case(iarms)
    case(0)
       !--No potential
    case(1)
       !--Cox&Gomez~~Altered TWA spiral pattern
       call s_potential(xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
    case(2,3)
       !--Pichardo&Martros~superposition of Schmidt Spheroids
       !--rho(a)=p0+p1*a, linear drop in density inside spheroids
       call pichardo_potential(xi,yi,zi,d2,ti,phii,phi,fextxi,fextyi,fextzi)
    case(4)
       !--Pichardo&Martros~superposition of Schmidt Spheroids
       !--rho(a)=p0+p1/a, inverse drop in density inside spheroids
       call schmidt_potential(xi,yi,zi,d2,ti,phii,phi,fextxi,fextyi,fextzi)
    end select

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=BARS-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case(ibar)
    case(0)
       !--No potential
    case(1)
       !--Long&Murali~smoothed biaxial bar
       call LMXbar(0.,1.,1.,xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
       !call LMXbar(90.*3.14159/180.,1.0,1.0,xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
    case(2)
       !--Long&Murali~smoothed triaxial bar
       call LMTbar(xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
    case(3)
       !--Wang~DwekG2+G3 bar+bulge, but able to do any distribution.
       call Orthog_basisbar(xi,yi,zi,r,dr,ti,hi,phi,fextxi,fextyi,fextzi)
    case(4)
       !--Dehnen~cosine_R^1/3 in/out
       call DehnenBar(xi,yi,d2,phii,ti,phi,fextxi,fextyi)
    case(5)
       !--Vogt&Letelier~smoothed S-shape
       call VogtSbar(xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
    case(6)
       !--Wada bar~Sinusoidal bar similar to Dehnen bar:
       call WadaBar(xi,yi,d2,phii,ti,hi,phi,fextxi,fextyi)
    end select

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=READIN-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case(iread)
    case(0)
       !--No potential
    case(1)
       !--Read in the potential from some gridded file.
       call BINReadPot3D(xi,yi,zi,ti,phi,fextxi,fextyi,fextzi)
    end select

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

 case(iext_densprofile)
    ! gravitational potential from tabulated density vs r profile
    call densityprofile_force(xi,yi,zi,fextxi,fextyi,fextzi,phi)

 case(iext_einsteinprec)
!
!--Potential given by (8) of NP2000, "Einstein precession"
!
    r2 = xi*xi + yi*yi + zi*zi + eps2_soft

    if (r2 > epsilon(r2)) then
       dr = 1./sqrt(r2)
       R_g = 1.0
       factor = 1. + 6.*R_g*dr
       dr3 = mass1*dr**3
       fextxi = fextxi - xi*dr3*factor
       fextyi = fextyi - yi*dr3*factor
       fextzi = fextzi - zi*dr3*factor
       phi    = -mass1*dr*(1. + 3.*R_g*dr)
    endif


 case(iext_gnewton)
!
!--Spatial component of the generalized Newtonian force
!
    call get_gnewton_spatial_force(xi,yi,zi,mass1,fextxi,fextyi,fextzi,phi)

 case(iext_staticsine)

!
!--Sinusoidal force, phi = A cos(k(x+B))
!

    call staticsine_force(xi,yi,fextxi,fextyi,fextzi,phi)


 case(iext_gwinspiral)
!
!--Gravitational wave inspiral
!
    if (present(ii)) then
       call get_gw_force_i(ii,fextxi,fextyi,fextzi,phi)
    else
       ! This will return 0 force, but this should not happen
       ! if initialised properly
       call get_gw_force_i(0, fextxi,fextyi,fextzi,phi)
    endif

 case(iext_discgravity)
!
!--vertical gravity in disc section
!
    phi = -mass1/sqrt(Rdisc**2 + yi**2)
    fextyi = -yi*mass1/sqrt(Rdisc**2 + yi**2)**3

 case(iext_corot_binary)
    !
    !--gravitational force from central binary
    !
    call binary_force(xi,yi,zi,ti,fextxi,fextyi,fextzi,phi)

    !
    !--spatial part of forces in corotating frame, i.e. centrifugal force
    !
    pos = (/xi,yi,zi/)
    call get_centrifugal_force(pos,fextxi,fextyi,fextzi,phi)

 case default
!
!--external forces should not be called if iexternalforce = 0
!
    call fatal('externalforces','external force not implemented',&
               var='iexternalforce',ival=iexternalforce)
 end select

!
!--return a timestep based only on the external force
!  so that we can do substeps with only the external force call
!
 if (present(dtf)) then
    f2i = fextxi*fextxi + fextyi*fextyi + fextzi*fextzi
    if (abs(f2i) > epsilon(f2i)) then
       !
       !--external force timestep based on sqrt(h/accel)
       !
       if (hi > epsilon(hi)) then
          dtf1 = sqrt(hi/sqrt(f2i))
       else
          dtf1 = huge(dtf1)
       endif
       !
       !--external force timestep based on sqrt(phi)/accel
       !
       if (abs(phi) > epsilon(phi)) then
          dtf2 = sqrt(abs(phi)/f2i)
       else
          dtf2 = huge(dtf2)
       endif
       dtf  = min(dtf1,dtf2)
       !if (dtf2 < dtf1) print*,' phi timestep = ',dtf2,' h/a = ',dtf1, ' ratio = ',dtf2/dtf1
    else
       dtf = huge(dtf)
    endif
 endif

 return
end subroutine externalforce

!-----------------------------------------------------------------------
!+
!  Query function to determine whether or not an external force
!  has velocity dependent parts
!+
!-----------------------------------------------------------------------
logical function is_velocity_dependent(iexternalforce)
 integer, intent(in) :: iexternalforce

 select case(iexternalforce)
 case(iext_corotate,iext_corot_binary,iext_prdrag,iext_lensethirring,iext_einsteinprec,iext_gnewton)
    is_velocity_dependent = .true.
 case default
    is_velocity_dependent = .false.
 end select

end function is_velocity_dependent

!-----------------------------------------------------------------------
!+
!  Returns the part of the external force that is velocity dependent
!  (these are treated implicitly in the leapfrog integrator)
!  This routine returns an explicit evaluation
!+
!-----------------------------------------------------------------------
subroutine externalforce_vdependent(iexternalforce,xyzi,veli,fexti,poti,densi,ui)
 use extern_corotate,      only:get_coriolis_force
 use extern_prdrag,        only:get_prdrag_vdependent_force
 use extern_lensethirring, only:get_lense_thirring_force
 use extern_gnewton,       only:get_gnewton_vdependent_force
 integer, intent(in)  :: iexternalforce
 real,    intent(in)  :: xyzi(3),veli(3)
 real,    intent(out) :: fexti(3)
 real,    intent(inout) :: poti
 real,    intent(in), optional :: densi,ui ! Needed for compatibility with gr

 select case(iexternalforce)
 case(iext_corotate,iext_corot_binary)
    call get_coriolis_force(xyzi,veli,fexti,poti)
 case(iext_prdrag)
    call get_prdrag_vdependent_force(xyzi,veli,mass1,fexti)
 case(iext_lensethirring,iext_einsteinprec)
    call get_lense_thirring_force(xyzi,veli,mass1,fexti)
 case(iext_gnewton)
    call get_gnewton_vdependent_force(xyzi,veli,mass1,fexti)
 case default
    fexti(:) = 0.
 end select

end subroutine externalforce_vdependent

!-----------------------------------------------------------------------
!+
!  Solves for velocity-dependent part of external force via an
!  implicit inversion of v^1 = v^1/2 + 1/2*dt*f0 + 1/2*dt*f1(x^1,v^1)
!  necessary for using v-dependent forces in leapfrog
!+
!-----------------------------------------------------------------------
subroutine update_vdependent_extforce(iexternalforce, &
           vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dkdt,xi,yi,zi,densi,ui)
 use extern_corotate,      only:update_coriolis
 use extern_prdrag,        only:update_prdrag
 use extern_lensethirring, only:update_ltforce
 use extern_gnewton,       only:update_gnewton
 integer, intent(in)    :: iexternalforce
 real,    intent(in)    :: dkdt,xi,yi,zi
 real,    intent(in)    :: vhalfx,vhalfy,vhalfz
 real,    intent(inout) :: fxi,fyi,fzi
 real,    intent(out)   :: fexti(3)
 real,    intent(in), optional :: densi,ui ! Needed for compatibility with gr

 select case(iexternalforce)
 case(iext_corotate,iext_corot_binary)
    call update_coriolis(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dkdt)
 case(iext_prdrag)
    call update_prdrag(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dkdt,xi,yi,zi,mass1)
 case(iext_lensethirring,iext_einsteinprec)
    call update_ltforce(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dkdt,xi,yi,zi,mass1)
 case(iext_gnewton)
    call update_gnewton(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,fexti,dkdt,xi,yi,zi,mass1)
 end select

end subroutine update_vdependent_extforce

!-----------------------------------------------------------------------
!+
!  update time-dependent external forces where necessary
!+
!-----------------------------------------------------------------------
subroutine update_externalforce(iexternalforce,ti,dmdt)
 use io,                only:warn
 use part,              only:xyzh,vxyzu,igas,npart,nptmass,&
                             xyzmh_ptmass,vxyz_ptmass
 use extern_gwinspiral, only:gw_still_inspiralling,get_gw_force
 use extern_binary,     only:update_binary
 integer, intent(in) :: iexternalforce
 real,    intent(in) :: ti,dmdt
 logical             :: stopped_now

 select case(iexternalforce)
 case(iext_binary,iext_corot_binary)
    call update_binary(ti)
 case(iext_gwinspiral)
    call gw_still_inspiralling(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,stopped_now)
    call get_gw_force()
    if (stopped_now) call warn('externalforces','Stars have merged. Disabling GW inspiral',2)
 end select

end subroutine update_externalforce

!-----------------------------------------------------------------------
!+
!  Performs checks to see whether or not a particle should be
!  accreted by crossing a boundary of the external potential
!  (at the moment this is just a hard boundary, but could in principle
!   add checks to see if particle is bound etc. here)
!+
!-----------------------------------------------------------------------
subroutine accrete_particles(iexternalforce,xi,yi,zi,hi,mi,ti,accreted,i)
 use extern_binary, only:binary_accreted,accradius1
 integer, intent(in)    :: iexternalforce
 real,    intent(in)    :: xi,yi,zi,mi,ti
 real,    intent(inout) :: hi
 logical, intent(out)   :: accreted
 integer, intent(in), optional :: i  ! for compatibility with GR routine
 real :: r2

 accreted = .false.
 select case(iexternalforce)
 case(iext_star,iext_prdrag,iext_lensethirring,iext_einsteinprec,iext_gnewton)

    r2 = xi*xi + yi*yi + zi*zi
    if (r2 < (accradius1)**2) accreted = .true.

 case(iext_binary,iext_corot_binary)

    accreted = binary_accreted(xi,yi,zi,mi,ti)

 end select

 if (accreted) then
    hi = -abs(hi)
 endif

end subroutine accrete_particles

!-----------------------------------------------------------------------
!+
!  query function for particles that have already been accreted
!+
!-----------------------------------------------------------------------
pure logical function was_accreted(iexternalforce,hi)
 integer, intent(in) :: iexternalforce
 real,    intent(in) :: hi

 select case(iexternalforce)
 case(iext_star,iext_binary,iext_corot_binary,iext_prdrag,&
      iext_lensethirring,iext_einsteinprec,iext_gnewton)
    ! An accreted particle is indicated by h < 0.
    ! Note less than, but not equal.
    ! (h=0 indicates dead MPI particle)
    was_accreted = (hi < 0.)
 case default
    was_accreted = .false.
 end select

end function was_accreted

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_externalforces(iunit,iexternalforce)
 use infile_utils,         only:write_inopt,get_optstring
 use extern_corotate,      only:write_options_corotate
 use extern_binary,        only:write_options_externbinary
 use extern_prdrag,        only:write_options_prdrag
 use extern_lensethirring, only:write_options_ltforce
 use extern_spiral,        only:write_options_spiral
 use extern_Bfield,        only:write_options_externB
 use extern_staticsine,    only:write_options_staticsine
 use extern_gwinspiral,    only:write_options_gwinspiral
 use extern_geopot,        only:write_options_geopot
 integer, intent(in) :: iunit,iexternalforce
 character(len=80) :: string

 write(iunit,"(/,a)") '# options relating to external forces'

 call get_optstring(iexternalforce_max,externalforcetype,string,4)
 call write_inopt(iexternalforce,'iexternalforce',trim(string),iunit)

 select case(iexternalforce)
 case(iext_star,iext_prdrag,iext_lensethirring,iext_einsteinprec,iext_gnewton,iext_geopot)
    call write_inopt(mass1,'mass1','mass of central object in code units',iunit)
    if (accradius1_hard < tiny(0.)) accradius1_hard = accradius1
    call write_inopt(accradius1,'accradius1','soft accretion radius of central object',iunit)
    call write_inopt(accradius1_hard,'accradius1_hard','hard accretion radius of central object',iunit)
 end select

 select case(iexternalforce)
 case(iext_star,iext_lensethirring,iext_einsteinprec,iext_gnewton)
    call write_inopt(eps_soft,'eps_soft','softening length (Plummer) for central potential in code units',iunit)
 end select

 select case(iexternalforce)
 case(iext_corotate)
    call write_options_corotate(iunit)
 case(iext_corot_binary)
    call write_options_corotate(iunit)
    call write_options_externbinary(iunit)
 case(iext_binary)
    call write_options_externbinary(iunit)
 case(iext_prdrag)
    call write_options_prdrag(iunit)
 case(iext_externB)
    call write_options_externB(iunit)
 case(iext_spiral)
    call write_options_spiral(iunit)
 case(iext_lensethirring,iext_einsteinprec)
    call write_options_ltforce(iunit)
 case(iext_staticsine)
    call write_options_staticsine(iunit)
 case(iext_gwinspiral)
    call write_options_gwinspiral(iunit)
 case(iext_geopot)
    call write_options_geopot(iunit)
 end select

end subroutine write_options_externalforces

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_extern(iexternalforce,hdr,time,ierr)
 use dump_utils,        only:dump_h,add_to_rheader
 use extern_binary,     only:write_headeropts_externbinary
 use extern_gwinspiral, only:write_headeropts_gwinspiral
 integer,      intent(in)    :: iexternalforce
 type(dump_h), intent(inout) :: hdr
 real,         intent(in)    :: time
 integer,      intent(out)   :: ierr

 select case(iexternalforce)
 case(iext_gwinspiral)
    call write_headeropts_gwinspiral(hdr,ierr)
 case(iext_binary,iext_corot_binary)
    call write_headeropts_externbinary(hdr,time,ierr)
 end select

end subroutine write_headeropts_extern

!-----------------------------------------------------------------------
!+
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_extern(iexternalforce,hdr,nptmass,ierr)
 use dump_utils,        only:dump_h
 use extern_binary,     only:read_headeropts_externbinary
 use extern_gwinspiral, only:read_headeropts_gwinspiral
 integer,      intent(in)  :: iexternalforce,nptmass
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr = 0
 select case(iexternalforce)
 case(iext_gwinspiral)
    call read_headeropts_gwinspiral(hdr,nptmass,ierr)
 case(iext_binary,iext_corot_binary)
    call read_headeropts_externbinary(hdr,ierr)
 end select

end subroutine read_headeropts_extern

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_externalforces(name,valstring,imatch,igotall,ierr,iexternalforce)
 use io,                   only:fatal,warn
 use extern_corotate,      only:read_options_corotate
 use extern_binary,        only:read_options_externbinary
 use extern_prdrag,        only:read_options_prdrag
 use extern_spiral,        only:read_options_spiral
 use extern_lensethirring, only:read_options_ltforce
 use extern_Bfield,        only:read_options_externB
 use extern_staticsine,    only:read_options_staticsine
 use extern_gwinspiral,    only:read_options_gwinspiral
 use extern_geopot,        only:read_options_geopot
 character(len=*), intent(in)    :: name,valstring
 logical,          intent(out)   :: imatch,igotall
 integer,          intent(out)   :: ierr
 integer,          intent(inout) :: iexternalforce
 integer, save :: ngot = 0
 logical :: igotallcorotate,igotallbinary,igotallprdrag
 logical :: igotallltforce,igotallspiral,igotallexternB
 logical :: igotallstaticsine,igotallgwinspiral,igotallgeopot
 character(len=30), parameter :: tag = 'externalforces'

 imatch            = .true.
 igotall           = .false.
 igotallcorotate   = .true.
 igotallbinary     = .true.
 igotallprdrag     = .true.
 igotallexternB    = .true.
 igotallspiral     = .true.
 igotallltforce    = .true.
 igotallstaticsine = .true.
 igotallgwinspiral = .true.
 igotallgeopot     = .true.

 !call read_inopt(db,'iexternalforce',iexternalforce,min=0,max=9,required=true)
 !if (imatch) ngot = ngot + 1

 select case(trim(name))
 case('iexternalforce')
    read(valstring,*,iostat=ierr) iexternalforce
    if (iexternalforce < 0) call fatal(tag,'silly choice of iexternalforce, use 0')
    ngot = ngot + 1
 case('mass1')
    read(valstring,*,iostat=ierr) mass1
    if (mass1 < 0)           call fatal(tag,'mass of central object cannot be -ve')
    if (mass1 < tiny(mass1)) call warn(tag,'mass of central object is zero')
    ngot = ngot + 1
 case('accradius1')
    read(valstring,*,iostat=ierr) accradius1
    if (iexternalforce <= 0) call warn(tag,'no external forces: ignoring accradius1 value')
    if (accradius1 < 0.)    call fatal(tag,'negative accretion radius')
 case('accradius1_hard')
    read(valstring,*,iostat=ierr) accradius1_hard
    if (iexternalforce <= 0) call warn(tag,'no external forces: ignoring accradius1_hard value')
    if (accradius1_hard > accradius1) call fatal(tag,'hard accretion boundary must be within soft accretion boundary')
 case('eps_soft')
    read(valstring,*,iostat=ierr) eps_soft
    if (iexternalforce <= 0) call warn(tag,'no external forces: ignoring accradius1 value')
    if (eps_soft < 0.)       call fatal(tag,'negative softening parameter',var='eps_soft',val=eps_soft)
    eps2_soft = eps_soft*eps_soft
 case default
    imatch = .false.
    select case(iexternalforce)
    case(iext_corotate)
       call read_options_corotate(name,valstring,imatch,igotallcorotate,ierr)
    case(iext_corot_binary)
       call read_options_corotate(name,valstring,imatch,igotallcorotate,ierr)
       call read_options_externbinary(name,valstring,imatch,igotallbinary,ierr)
    case(iext_binary)
       call read_options_externbinary(name,valstring,imatch,igotallbinary,ierr)
    case(iext_prdrag)
       call read_options_prdrag(name,valstring,imatch,igotallprdrag,ierr)
    case(iext_externB)
       call read_options_externB(name,valstring,imatch,igotallexternB,ierr)
    case(iext_spiral)
       call read_options_spiral(name,valstring,imatch,igotallspiral,ierr)
    case(iext_lensethirring,iext_einsteinprec)
       call read_options_ltforce(name,valstring,imatch,igotallltforce,ierr)
    case(iext_staticsine)
       call read_options_staticsine(name,valstring,imatch,igotallstaticsine,ierr)
    case(iext_gwinspiral)
       call read_options_gwinspiral(name,valstring,imatch,igotallgwinspiral,ierr)
    case(iext_geopot)
       call read_options_geopot(name,valstring,imatch,igotallgwinspiral,ierr)
    end select
 end select
 igotall = (ngot >= 1      .and. igotallcorotate   .and. &
            igotallbinary  .and. igotallprdrag     .and. &
            igotallspiral  .and. igotallltforce    .and. &
            igotallexternB .and. igotallstaticsine .and. &
            igotallgwinspiral .and. igotallgeopot)

 !--make sure mass is read where relevant
 select case(iexternalforce)
 case(iext_star,iext_lensethirring,iext_einsteinprec,iext_gnewton,iext_geopot)
    igotall = igotall .and. (ngot >= 2)
 end select

end subroutine read_options_externalforces

!-----------------------------------------------------------------------
!+
!  interface to initialisation/checking routines for external forces
!+
!-----------------------------------------------------------------------
subroutine initialise_externalforces(iexternalforce,ierr)
 use io,                   only:error
 use extern_lensethirring, only:check_lense_thirring_settings
 use extern_spiral,        only:initialise_spiral
 use extern_densprofile,   only:load_extern_densityprofile
 use extern_Bfield,        only:check_externB_settings
 use extern_gwinspiral,    only:initialise_gwinspiral
 use units,                only:G_is_unity,c_is_unity,get_G_code,get_c_code
 use part,                 only:npart,nptmass
 integer, intent(in)  :: iexternalforce
 integer, intent(out) :: ierr

 ierr = 0
 select case(iexternalforce)
 case(iext_lensethirring,iext_einsteinprec)
    call check_lense_thirring_settings(ierr,accradius1)
 case(iext_externB)
    call check_externB_settings(ierr)
 case(iext_spiral)
    call initialise_spiral(ierr)
 case(iext_densprofile)
    call load_extern_densityprofile(ierr)
 case (iext_gwinspiral)
    call initialise_gwinspiral(npart,nptmass,ierr)
    if (ierr > 0) then
       call error('externalforces','Require number of particles per star or sink particles',&
                  var='iexternalforce',ival=iexternalforce)
       ierr = ierr + 1
    endif
 case default
    if (iexternalforce <= 0 .or. iexternalforce > iexternalforce_max) then
       call error('externalforces','externalforce not implemented',var='iexternalforce',ival=iexternalforce)
       ierr = ierr + 1
    endif
 end select

 select case(iexternalforce)
 case(iext_star,iext_binary,iext_corot_binary,iext_prdrag,iext_spiral,iext_lensethirring,iext_einsteinprec,iext_gnewton)
    !
    !--check that G=1 in code units
    !
    if (.not.G_is_unity()) then
       call error('units',trim(externalforcetype(iexternalforce))//&
                  ' external force assumes G=1 in code units but we have',var='G',val=real(get_G_code()))
       ierr = ierr + 1
    endif
 end select

 select case(iexternalforce)
 case(iext_prdrag)
    !
    !--check that c=1 in code units for prdrag only
    !
    if (.not.c_is_unity()) then
       call error('units',trim(externalforcetype(iexternalforce))//&
                  ' external force assumes c=1 in code units but we have',var='c',val=real(get_c_code()))
       ierr = ierr + 1
    endif
 end select

end subroutine initialise_externalforces

end module externalforces
