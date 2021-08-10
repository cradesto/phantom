!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!

module analysis
!
! analysis programmes used to analys a neutron star merger.
!  This is an interactive routine that includes multiple analysis options
!  Note: all outputs are in code units unless explicitly stated
!  Author: Bernard Field & Madeline Marshall (supervisors: James Wurster & Paul Lasky)
!  Changes for stripping were done by Marat Potashov
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, extern_gwinspiral, io, part, physcon,
!   prompting, readwrite_dumps, units
!
  use io,                only: fatal
  use part,              only: rhoh
  use physcon,           only: pi
  use centreofmass,      only: get_centreofmass
  use readwrite_dumps,   only: opened_full_dump
  use extern_gwinspiral, only: Nstar,&
                                get_momentofinertia,&
                                calculate_omega,&
                                com, vcom,&
                                evector_old,&
                                omega_old,&
                                time_old

  implicit none

  private

  character(len=20), parameter, public :: analysistype = 'Stripping'
  !
  ! integer, parameter :: nana_opts          = 20
  real               :: density_cutoff_cgs = 5.0d-5     ! The density threshold in code units (opts 2,3,4)
  real               :: thickness          = 2.0
  real               :: dtheta             = 5.*pi/180. ! width of slice in the directions of the minor and major axes (opt 4)

  logical            :: firstcall          = .true.
  logical            :: iexist
  character(len=200) :: fileout
  ! character(len=200) :: analysis_opt(nana_opts)

  ! integer,            :: choice

  ! for L1 calculation
  real               :: L1_projection
  real, pointer      :: xyzh_(:,:), vxyzu_(:,:)
  real               :: particlemass_
  integer            :: npart_

  public             :: do_analysis

contains
!--------------------------------------------------------------------------
  subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

    use centreofmass, only: get_centreofmass, reset_centreofmass
    use prompting,    only: prompt
    use units,        only: unit_density

    character(len=*), intent(in)    :: dumpfile
    integer,          intent(in)    :: num, npart, iunit
    real, target,     intent(inout) :: xyzh(:,:), vxyzu(:,:)
    real,             intent(in)    :: particlemass, time

    real                            :: density_cutoff

    if(firstcall) then
      evector_old = 0.
      omega_old = 0.
      time_old = 0.
    endif

    xyzh_ => xyzh
    vxyzu_ => vxyzu
    particlemass_ = particlemass
    npart_ = npart

    !
    !--Determine which analysis to run if this is the first step
    !

    !--Reset centre of mass
    call reset_centreofmass(npart, xyzh, vxyzu)

    !--Calculate the centre of mass and velocity of the system
    call get_centreofmass(com, vcom, npart, xyzh, vxyzu)

    density_cutoff = density_cutoff_cgs / unit_density

    !--Calculate the moment of inertia tensor
    call calculate_I(dumpfile, xyzh, vxyzu, time, npart, density_cutoff, iunit, particlemass)

    call trace_com(dumpfile, xyzh, vxyzu, time, npart, iunit, particlemass)

    ! if(firstcall) then
    !   analysis_opt(:) = 'none'
    !   analysis_opt(1) = 'Trace centre of mass of each star and the system'
    !   analysis_opt(2) = 'Brute force calculation of T/W'
    !   analysis_opt(3) = 'Calculate the moment of inertia tensor'
    !   analysis_opt(4) = 'Calculate the radial profile of a slice through the midplane'
    !   write(*,"(a)") 'Analysis options: '
    !   do i = 1, nana_opts
    !     if(trim(analysis_opt(i)) /= 'none') write(*,"(a5,i2,1x,a60)") 'Case ', i, analysis_opt(i)
    !   enddo
    !   choice = 1
    !   call prompt('Enter analysis choice',choice,1,nana_opts)
    !   !
    !   !--Get the index range for each star
    !   ! (Note from DJP: This should now be automatically read from the dump header)
    !   if(Nstar(1) <= 0) call fatal('analysis_NSmerger','Require Nstar(1) > 0 in header of dump file')
    !   if(Nstar(2) <= 0) call fatal('analysis_NSmerger','Require Nstar(2) > 0 in header of dump file')
    !   !
    !   !--Prompt for the density_cut off
    !   if(choice > 1) then
    !     call prompt('Enter cutoff density (cgs):',density_cutoff_cgs,0.)
    !     density_cutoff = density_cutoff_cgs / unit_density
    !   endif
    !   !--Prompt for thickness
    !   if(choice == 4) then
    !     call prompt('Enter the thickness of the mid-plane slice (code units):',thickness,0.)
    !   endif
    ! endif
    ! !
    ! !--Calculate the centre of mass and velocity of the system
    ! !
    ! call get_centreofmass(com,vcom,npart,xyzh,vxyzu)
    ! !
    ! !--Run the analysis
    ! !
    ! select case(choice)
    !  case(1)
    !   !--Trace the centre of masses
    !   call trace_com(dumpfile,xyzh,vxyzu,time,npart,iunit)
    !  case(2)
    !   !--Calculate T/W using brute force
    !   call calculate_TW(dumpfile,xyzh,vxyzu,time,npart,iunit,particlemass)
    !  case(3)
    !   !--Calculate the moment of inertia tensor
    !   call calculate_I(dumpfile,xyzh,time,npart,iunit,particlemass)
    !  case(4)
    !   !--Calculate the radial profile of a slice through the midplane'
    !   call calculate_midplane_profile(dumpfile,xyzh,vxyzu,npart,iunit,particlemass)
    ! end select
    ! !
    ! close(iunit)
    ! firstcall = .false. ! done here since this logical is required for opening files

    if(firstcall) firstcall = .false. ! performed here since multiple commands require this knowledge

  end subroutine do_analysis
!-----------------------------------------------------------------------
!+
! Trace centre of mass of each star and the system
!+
!-----------------------------------------------------------------------
  subroutine trace_com(dumpfile,xyzh,vxyzu,time,npart,iunit,particlemass)

    use dim,          only: maxp, maxvxyzu
    use centreofmass, only: get_centreofmass, get_total_angular_momentum

    character(len=*), intent(in) :: dumpfile
    integer,          intent(in) :: npart,iunit
    real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
    real,             intent(in) :: time
    real,             intent(in) :: particlemass

    integer                      :: i, iA, iB

    real(kind=8)                 :: xposA(3), xposB(3), vpos(3), sep
    integer                      :: npartA, npartB
    real                         :: xyzhA(4,maxp), vxyzuA(maxvxyzu, maxp)
    real                         :: xyzhB(4,maxp), vxyzuB(maxvxyzu, maxp)
    real                         :: L_A(3), L_B(3), L_tot(3)

    !--Open file (appendif exists)
    fileout = trim(dumpfile(1:index(dumpfile,'_')-1))//'_centres.dat'

    inquire(file=fileout,exist=iexist)

    if(firstcall .or. .not.iexist) then
      open(iunit,file=fileout,status='replace')
      write(iunit,"('#',26(1x,'[',i2.2,1x,a11,']',2x))") &
        1, 'time',   &
        2, 'xcom',   &
        3, 'ycom',   &
        4, 'zcom',   &
        5, 'x_A',    &
        6, 'y_A',    &
        7, 'z_A',    &
        8, 'x_B',    &
        9, 'y_B',    &
        10,'z_B',    &
        11,'d',      &
        12,'m_A',    &
        13,'m_B',    &
        14,'dm',     &
        15,'L_tot,1',&
        16,'L_tot,2',&
        17,'L_tot,3',&
        18,'n L_tot',&
        19,'L_A,1',  &
        20,'L_A,2',  &
        21,'L_A,3',  &
        22,'n L_A',  &
        23,'L_B,1',  &
        24,'L_B,2',  &
        25,'L_B,3',  &
        26,'n L_B'
    else
      open(iunit,file=fileout,position='append')
    endif

!--Determine the centre of mass of each star
!--NB: In this version two center masses of star should be lie on x axe!

    iA = 1
    iB = 1
    do i = 1, npart
      if(dot_product(xyzh(1:3,i), evector_old) >= L1_projection) then
        xyzhA(:,iA) = xyzh(:,i)
        vxyzuA(:,iA) = vxyzu(:,i)
        iA = iA + 1
      else
        xyzhB(:,iB) = xyzh(:,i)
        vxyzuB(:,iB) = vxyzu(:,i)
        iB = iB + 1
      endif
    enddo

    npartA = iA - 1
    npartB = iB - 1

    call get_centreofmass(xposA, vpos, npartA, xyzhA, vxyzuA)
    call get_centreofmass(xposB, vpos, npartB, xyzhB, vxyzuB)
    sep = sqrt(dot_product(xposA - xposB, xposA - xposB))

    call get_total_angular_momentum(xyzhA, vxyzuA, npartA, L_A)
    call get_total_angular_momentum(xyzhB, vxyzuB, npartB, L_B)
    call get_total_angular_momentum(xyzh, vxyzu, npart, L_tot)

    write(iunit,'(26(es18.10,1x))')&
      time,&
      com,&
      xposA,&
      xposB,&
      sep,&
      npartA*particlemass,&
      npartB*particlemass,&
      abs(npartA*particlemass - npartB*particlemass),&
      L_tot,&
      norm2(L_tot),&
      L_A,&
      norm2(L_A),&
      L_B,&
      norm2(L_B)

    close(iunit)
    !
  end subroutine trace_com
!-----------------------------------------------------------------------
!+
! Calculate T/W using brute force
!+
!-----------------------------------------------------------------------
  subroutine calculate_TW(dumpfile,xyzh,vxyzu,time,npart,density_cutoff,iunit,particlemass)

    character(len=*), intent(in) :: dumpfile
    integer,          intent(in) :: npart,iunit
    real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
    real,             intent(in) :: density_cutoff,particlemass,time
    integer                      :: i,j,npartmeasured
    real                         :: rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2,rad2,rmax2
    real                         :: erot,erotx,eroty,erotz,grav
    real                         :: r(3),v(3)
    !
    !--Skip if not a full dump
    if(.not.opened_full_dump) return
    !
    !--Open file (appendif exists)
    fileout = trim(dumpfile(1:index(dumpfile,'_')-1))//'_TW.dat'
    inquire(file=fileout,exist=iexist)
    if(firstcall .or. .not.iexist) then
      open(iunit,file=fileout,status='replace')
      write(iunit,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'time',       &
        2,'T/W',     &
        3,'T',     &
        4,'W',    &
        5,'Star mass', &
        6,'Star radius'
    else
      open(iunit,file=fileout,position='append')
    endif
    !
    !--Calculate rotational kinetic energy
    npartmeasured = 0
    rmax2 = 0.
    erot  = 0.
    erotx = 0.
    eroty = 0.
    erotz = 0.
!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,particlemass,com,vcom,density_cutoff) &
!$omp private(i,r,v,rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2,rad2) &
!$omp reduction(+:erotx,eroty,erotz,npartmeasured) &
!$omp reduction(max:rmax2)
!$omp do
    do i=1,npart
      if(rhoh(xyzh(4,i),particlemass) > density_cutoff) then
        r = xyzh(1:3,i)  - com
        v = vxyzu(1:3,i) - vcom
        ! r cross v
        rcrossvx = (r(2)*v(3) - r(3)*v(2))
        rcrossvy = (r(3)*v(1) - r(1)*v(3))
        rcrossvz = (r(1)*v(2) - r(2)*v(1))
        ! rotational energy around each axis through the origin
        radxy2 = r(1)*r(1) + r(2)*r(2)
        radyz2 = r(3)*r(3) + r(2)*r(2)
        radxz2 = r(1)*r(1) + r(3)*r(3)
        if(radyz2 > 0.) erotx = erotx + particlemass*rcrossvx*rcrossvx/radyz2
        if(radxz2 > 0.) eroty = eroty + particlemass*rcrossvy*rcrossvy/radxz2
        if(radxy2 > 0.) erotz = erotz + particlemass*rcrossvz*rcrossvz/radxy2
        !
        ! size of the star
        rad2 = dot_product(r,r)
        rmax2 = max(rmax2,rad2)
        !
        ! additional bookkeeping
        npartmeasured = npartmeasured + 1
      endif
    enddo
!$omp enddo
!$omp end parallel
    erotx = 0.5*erotx
    eroty = 0.5*eroty
    erotz = 0.5*erotz
    erot  = sqrt(erotx**2 + eroty**2 + erotz**2)
    !
    !--Calculate gravitational potential energy
    grav = 0.
!$omp parallel default(none) &
!$omp shared(npart,xyzh,particlemass,density_cutoff) &
!$omp private(i,j,r,rad2) &
!$omp reduction(+:grav)
!$omp do
    do i=1,npart
      if(rhoh(xyzh(4,i),particlemass) > density_cutoff) then
        do j=i+1,npart
          if(rhoh(xyzh(4,j),particlemass) > density_cutoff) then
            r    = xyzh(1:3,i) - xyzh(1:3,j)
            rad2 = dot_product(r,r)
            if(rad2 > 0.0) grav = grav + 1.0/sqrt(rad2)
          endif
        enddo
      endif
    enddo
!$omp enddo
!$omp end parallel
    grav = -grav*particlemass*particlemass
    !
    !--Write results
    print *, "time:",time," T/|W|:", erot/abs(grav), " T:", erot, " W:", grav, "ignored particles:", npart - npartmeasured, &
      "star mass:", npartmeasured*particlemass, "rstar:",sqrt(rmax2)
    write(iunit,'(6(es18.10,1x))') time,erot/abs(grav),erot,grav,npartmeasured*particlemass,sqrt(rmax2)
    !
  end subroutine calculate_TW
!-----------------------------------------------------------------------
!+
!  Determines the moment of inertia tensor, in diagonalised form.
!  Can exclude particles which are below a specified cut-off density.
!  Will output the mass and radius of the measured area.
!  Will output the eigenvalues and eigenvectors of inertia tensor
!   and L1 point coordinates
!+
!-----------------------------------------------------------------------
  subroutine calculate_I(dumpfile,xyzh,vxyzu,time,npart,density_cutoff,iunit,particlemass)

    character(len=*), intent(in) :: dumpfile
    real,             intent(in) :: time
    integer,          intent(in) :: npart
    real,             intent(in) :: xyzh(:,:), vxyzu(:,:)
    real,             intent(in) :: density_cutoff
    integer,          intent(in) :: iunit
    real,             intent(in) :: particlemass

    integer                      :: npartused
    real                         :: rmax, smallI, medI, bigI
    integer                      :: smallIIndex, middleIIndex, bigIIndex
    real                         :: principle(3), evectors(3,3), ellipticity(2)
    real                         :: omega(3)
    real                         :: omega_mean(3)
    real                         :: L1(3)

    !--Open file (appendif exists)
    fileout = trim(dumpfile(1:index(dumpfile,'_')-1))//'_inertia.dat'
    inquire(file=fileout,exist=iexist)
    if(firstcall .or. .not.iexist) then
      open(iunit,file=fileout,status='replace')
      write(iunit,"('#',25(1x,'[',i2.2,1x,a11,']',2x))") &
        1, 'time',           &
        2, 'I1',             &
        3, 'I2',             &
        4, 'I3',             &
        5, 'e1',             &
        6, 'e2',             &
        7, 'v1,1',           &
        8, 'v1,2',           &
        9, 'v1,3',           &
        10,'v2,1',           &
        11,'v2,2',           &
        12,'v2,3',           &
        13,'v3,1',           &
        14,'v3,2',           &
        15,'v3,3',           &
        16,'ang v1',         &
        17,'ang v2',         &
        18,'ang v3',         &
        19,'ang v',          &
        20,'L1,1',           &
        21,'L1,2',           &
        22,'L1,3',           &
        23,'excluded parts', &
        24,'mstar',          &
        25,'rstar'
    else
      open(iunit,file=fileout,position='append')
    endif
    !
    !--Calculate the tensor
    call get_momentofinertia(xyzh, npart, density_cutoff, particlemass, npartused, principle, evectors, rmax)
    !
    !--Sort the principle moments, since ellipticity depends on it.
    smallIIndex = minloc(principle, dim=1)
    bigIIndex = maxloc(principle, dim=1)
    middleIIndex = 6 - smallIIndex - bigIIndex

    smallI = principle(smallIIndex)
    medI   = principle(middleIIndex)
    bigI   = principle(bigIIndex)

    ellipticity(1) = sqrt(2.0*(bigI-smallI)/smallI)
    ellipticity(2) = sqrt(2.0*(bigI-medI)/medI)

    ! Calculate omega
    omega = calculate_omega(evectors(:, smallIIndex), evector_old, time, time_old, omega_old)
    write(*,*) "Omega coords = ", omega
    write(*,*) "Omega norm2 = ", norm2(omega)

    evector_old = evectors(:, smallIIndex)
    omega_old = omega
    time_old = time

    ! NB: only for full dumpfile
    omega_mean = calculate_mean_omega(npart, xyzh, vxyzu)
    write(*,*) "Mean Omega coords = ", omega_mean
    write(*,*) "Mean Omega norm2 = ", norm2(omega_mean)

    call L1_point(2, xyzh, particlemass, npart, L1_projection, L1)

    !--Write to file
    write(iunit,'(25(es18.10,1x))') &
      time, principle(smallIIndex), principle(middleIIndex), principle(bigIIndex),&
      ellipticity(1), ellipticity(2),&
      evectors(1,smallIIndex), evectors(2,smallIIndex), evectors(3,smallIIndex),&
      evectors(1,middleIIndex), evectors(2,middleIIndex), evectors(3,middleIIndex),&
      evectors(1,bigIIndex), evectors(2,bigIIndex), evectors(3,bigIIndex),&
      omega(1), omega(2), omega(3), norm2(omega),&
      L1(1), L1(2), L1(3),&
      real(npart-npartused), npartused*particlemass, rmax

  end subroutine calculate_I
!-------------------------------------------------------------
! Function to find mean omega vector from a list
! of positions and velocities
!-------------------------------------------------------------
  function calculate_mean_omega(n,xyz,vxyz) result(l)

    use vectorutils, only:cross_product3D

    integer, intent(in) :: n
    real,    intent(in) :: xyz(:,:),vxyz(:,:)
    real    :: l(3)

    real    :: li(3)
    integer :: i

    l = 0.
    do i = 1,n
      call cross_product3D(xyz(:,i),vxyz(:,i),li)
      li = li/norm2(xyz(:,i))**2
      l = l + li
    enddo
    l = l/real(n)

  end function calculate_mean_omega
!-----------------------------------------------------------------------
!+
!  Determines the radial profile of the midplane (of a given thickness)
!+
!-----------------------------------------------------------------------
  subroutine calculate_midplane_profile(dumpfile,xyzh,vxyzu,npart,density_cutoff,iunit,particlemass)

    use part, only: alphaind

    character(len=*), intent(in) :: dumpfile
    integer,          intent(in) :: npart,iunit
    real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
    real,             intent(in) :: density_cutoff
    real,             intent(in) :: particlemass
    integer,          parameter  :: nbins = 800
    integer                      :: i,j,npartused,zloc,majorloc,minorloc
    real                         :: rmax,thetamajor,thetaminor
    integer                      :: bincountmaj(nbins),bincountmin(nbins),bincountavg(nbins)
    real                         :: rtocm(npart),theta(npart),vtheta(npart),angv(npart)
    real                         :: avvinbinmaj(nbins),avvinbinmin(nbins),avvinbinavg(nbins)
    real                         :: vinbinmaj(nbins),  vinbinmin(nbins),  vinbinavg(nbins)
    real                         :: alphabinmaj(nbins),alphabinmin(nbins),alphabinavg(nbins)
    real                         :: partdensmaj(nbins),partdensmin(nbins),partdensavg(nbins)
    real                         :: radbin(nbins),vol(nbins)
    real                         :: principle(3),principlenew(2),evectors(3,3),evectorsnew(3,2)
    real                         :: major(3),minor(3)
    !
    !--Skip if not a full dump
    if(.not.opened_full_dump)return
    !
    !--Open file
    fileout = trim(dumpfile(1:index(dumpfile,'_')-1))//'_rotataxesprofile'//trim(dumpfile(index(dumpfile,'_'):))//'.dat'
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',10(1x,'[',i2.2,1x,a11,']',2x))") &
      1,'outer bin rad',&
      2,'density_min',  &
      3,'ang v_minor',  &
      4,'alpha_minor',  &
      5,'density_maj',  &
      6,'ang v_major',  &
      7,'alpha_major',  &
      8,'density_tot',  &
      9,'ang v_tot',    &
      10,'alpha_tot'
    !
    !--Initialise variables
    !
    theta       = 0.0
    rtocm       = 0.0
    vtheta      = 0.0
    angv        = 0.0
    vol         = 0.0
    bincountmaj = 0
    vinbinmaj   = 0.0
    avvinbinmaj = 0.0
    alphabinmaj = 0.0
    partdensmaj = 0.0
    bincountmin = 0
    vinbinmin   = 0.0
    avvinbinmin = 0.0
    alphabinmin = 0.0
    partdensmin = 0.0
    bincountavg = 0
    vinbinavg   = 0.0
    avvinbinavg = 0.0
    alphabinavg = 0.0
    partdensavg = 0.0
    !
    !--Calculate radius and angle of particle from CoM coordinates in x-y plane
!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,com,vcom,rtocm,theta,vtheta,angv) &
!$omp private(i)
!$omp do
    do i=1,npart
      rtocm (i) = sqrt((xyzh(1,i)-com(1))**2 + (xyzh(2,i)-com(2))**2)
      theta (i) = atan2((xyzh(2,i)-com(2)),(xyzh(1,i)-com(1)))
      vtheta(i) = -(vxyzu(1,i)-vcom(1))*sin(theta(i)) + (vxyzu(2,i)-vcom(2))*cos(theta(i))
      angv  (i) = abs(vtheta(i))/rtocm(i)
    enddo
!$omp enddo
!$omp end parallel
    !
    !--Calculate moment of inertia
    call get_momentofinertia(xyzh, npart, density_cutoff, particlemass, npartused, principle, evectors, rmax)
    !
    !--Find location of major and minor axes
    zloc = maxloc(evectors(3,:),1)
    j = 1
    do i = 1,3
      if(i/=zloc) then
        principlenew( j) =principle( i)
        evectorsnew(:,j) =evectors(:,i)
        j = j+1
      endif
    enddo
    majorloc   = minloc(principlenew,1)
    minorloc   = maxloc(principlenew,1)
    major      = evectorsnew(:,majorloc)
    minor      = evectorsnew(:,minorloc)
    thetamajor = atan2(major(2),major(1))
    thetaminor = atan2(minor(2),minor(1))
    write(*,*) 'evectors',evectors
    write(*,*) 'evectorsnew',evectorsnew
    write(*,*) 'major',major
    write(*,*) 'minor',minor
    write(*,*) 'thetamajor',thetamajor
    write(*,*) 'thetaminor',thetaminor
    !
    !Set radii and calculate volume of slice bins
    rmax = maxval(rtocm)
    do i = 1,nbins
      radbin(i) = rmax*float(i)/float(nbins)
      if(i==1) then
        vol(i) = thickness*dtheta*radbin(1)**2
      else
        vol(i) = thickness*dtheta*(radbin(i)**2-radbin(i-1)**2)
      endif
    enddo
    !
    !--Sort particles into bins, find total angular velocity and total alpha in each
    do i=1,npart
      do j=1,nbins
        if(rtocm(i)<radbin(j)) then
          if(abs(xyzh(3,i)-com(3)) < 0.5*thickness) then
            bincountavg(j) = bincountavg(j) + 1
            vinbinavg  (j) = vinbinavg  (j) + angv(i)
            alphabinavg(j) = alphabinavg(j) + alphaind(1,i)
            if(theta(i) < (thetamajor+dtheta) .and. theta(i) > (thetamajor-dtheta)) then
              bincountmaj(j) = bincountmaj(j) + 1
              vinbinmaj  (j) = vinbinmaj  (j) + angv(i)
              alphabinmaj(j) = alphabinmaj(j) + alphaind(1,i)
            elseif (theta(i) < (thetaminor+dtheta) .and. theta(i) > (thetaminor-dtheta)) then
              bincountmin(j) = bincountmin(j) + 1
              vinbinmin  (j) = vinbinmin  (j) + angv(i)
              alphabinmin(j) = alphabinmin(j) + alphaind(1,i)
            endif
          endif
          exit
        endif
      enddo
    enddo
    !
    !--Convert totals to averages for each bin
    do i = 1,nbins
      if(bincountmaj(i) > 0) then
        avvinbinmaj(i) = vinbinmaj(i)  /float(bincountmaj(i))
        alphabinmaj(i) = alphabinmaj(i)/float(bincountmaj(i))
        partdensmaj(i) = float(bincountmaj(i))*particlemass/vol(i)
      endif
      if(bincountmin(i) > 0) then
        avvinbinmin(i) = vinbinmin(i)  /float(bincountmin(i))
        alphabinmin(i) = alphabinmin(i)/float(bincountmin(i))
        partdensmin(i) = float(bincountmin(i))*particlemass/vol(i)
        print*, partdensmin(i) ,float(bincountmin(i)),particlemass,vol(i)
      endif
      if(bincountavg(i) > 0) then
        avvinbinavg(i) = vinbinavg(i)  /float(bincountavg(i))
        alphabinavg(i) = alphabinavg(i)/float(bincountavg(i))
        partdensavg(i) = float(bincountavg(i))*particlemass/(vol(i)*pi/dtheta)
      endif
    enddo
    !
    !--Write results to file
    do i = 1,nbins
      write(iunit,'(10(es18.10,1x))')  radbin(i), partdensmin(i), avvinbinmin(i), alphabinmin(i), &
        partdensmaj(i), avvinbinmaj(i), alphabinmaj(i), partdensavg(i), avvinbinavg(i), alphabinavg(i)
    enddo
!
  end subroutine calculate_midplane_profile
!-----------------------------------------------------------------------
! line - is a unit vector in the direction of the line
! point - vector to the point
!-----------------------------------------------------------------------
  real function distance_to_line(line, point) result(d)

    ! use vectorutils, only:cross_product3D

    real, intent(in) :: line(3), point(3)

    real :: p, vd(3)

    p = dot_product(-point, line)
    vd = -point - p*line
    d = norm2(vd)

    ! or

    ! call cross_product3D(-point,line,vd)
    ! d = norm2(vd)/norm2(line)

  end function distance_to_line
!-----------------------------------------------------------------------
! Finding the value of the gravitational potential by Lagrange interpolation method
! NB: Due to the factorial growth in the numerator and denominator this method apears unsuccessful!
!-----------------------------------------------------------------------
  real function interpolate_potential_lagrange(p, x, y, npoints) result(phi)

    real,    intent(in) :: p
    real,    intent(in) :: x(:)
    real,    intent(in) :: y(:)
    integer, intent(in) :: npoints

    integer             :: i, k
    real                :: g

    real                :: g1
    real                :: g2

    phi = 0.
    do i = 1, npoints
      g = 1.
      do k = 1, npoints
        if(k /= i) then
          g1 = p - x(k)
          g2 = x(i) - x(k)
          g = g*g1/g2
        endif
      enddo
      phi = phi + y(i)*g
    enddo

  end function interpolate_potential_lagrange
!-----------------------------------------------------------------------
  subroutine gravitational_potential_wrapper(p, potential)

    real, intent (in)  :: p
    real, intent (out) :: potential

    call gravitational_potential(p, xyzh_, particlemass_, npart_, potential)

  end subroutine gravitational_potential_wrapper
!-----------------------------------------------------------------------
  subroutine gravitational_potential(p, xyzh, particlemass, npart, potential)

    real,           intent (in)  :: p
    integer,        intent (in)  :: npart
    real,           intent (in)  :: xyzh(4,npart)
    real,           intent (in)  :: particlemass
    real,           intent (out) :: potential

    integer                      :: i
    real(kind=8)                 :: point(3)
    real(kind=8)                 :: dpoint(3)
    real(kind=8)                 :: dr ! 1/sqrt(r^2)

    potential = 0.

    point = p*evector_old

    do i = 1, npart

      dpoint = point - xyzh(1:3,i)
      dr = 1./norm2(dpoint)
      potential = potential - dr

    enddo

    potential = potential*particlemass

  end subroutine gravitational_potential
!-----------------------------------------------------------------------
  subroutine roche_potential_wrapper(p, potential)

    real, intent (in)  :: p
    real, intent (out) :: potential

    real               :: omega

    omega = norm2(omega_old)
    call roche_potential(p, xyzh_, particlemass_, npart_, omega, potential)

  end subroutine roche_potential_wrapper
!-----------------------------------------------------------------------
  subroutine roche_potential(p, xyzh, particlemass, npart, omega, potential)

    real,    intent (in)  :: p
    integer, intent (in)  :: npart
    real,    intent (in)  :: xyzh(4,npart)
    real,    intent (in)  :: particlemass
    real,    intent (in)  :: omega
    real,    intent (out) :: potential

    potential = 0.

    call gravitational_potential(p, xyzh, particlemass, npart, potential)

    potential = potential - 0.5*(omega*omega)*(p*p)

  end subroutine roche_potential
!-----------------------------------------------------------------------
  subroutine gravitational_force_wrapper(p, force, dforce)

    real, intent (in)  :: p
    real, intent (out) :: force
    real, intent (out) :: dforce

    call gravitational_force(p, xyzh_, particlemass_, npart_, force, dforce)

  end subroutine gravitational_force_wrapper
!-----------------------------------------------------------------------
  subroutine gravitational_force(p, xyzh, particlemass, npart, force, dforce)

    real,           intent (in)  :: p
    integer,        intent (in)  :: npart
    real,           intent (in)  :: xyzh(4,npart)
    real,           intent (in)  :: particlemass
    real,           intent (out) :: force
    real,           intent (out) :: dforce

    integer                      :: i
    real(kind=8)                 :: point(3)
    real(kind=8)                 :: dpoint(3)
    real(kind=8)                 :: f(3)
    real(kind=8)                 :: df(6)
    real(kind=8)                 :: dr  ! 1/sqrt(r^2)
    real(kind=8)                 :: dr3 ! 1/sqrt(r^2)^3
    real(kind=8)                 :: dr5 ! 1/sqrt(r^2)^5

    force = 0.
    dforce = 0.

    f = 0.
    df = 0.

    point = p*evector_old

    do i = 1, npart

      dpoint = point - xyzh(1:3,i)
      dr = 1./norm2(dpoint)
      dr3 = dr*dr*dr
      dr5 = dr3*dr*dr

      f = f - dpoint*dr3

      ! NB: Check for correctness
      df(1) = df(1) + dr5*(3.*dpoint(1)*dpoint(1) - 1.) ! dfx/dx
      df(2) = df(2) + dr5*(3.*dpoint(1)*dpoint(2))      ! dfx/dy = dfy/dx
      df(3) = df(3) + dr5*(3.*dpoint(1)*dpoint(3))      ! dfx/dz = dfz/dx
      df(4) = df(4) + dr5*(3.*dpoint(2)*dpoint(2) - 1.) ! dfy/dy
      df(5) = df(5) + dr5*(3.*dpoint(2)*dpoint(3))      ! dfy/dz = dfz/dy
      df(6) = df(6) + dr5*(3.*dpoint(3)*dpoint(3) - 1.) ! dfz/dz

    enddo

    force = norm2(f)*particlemass
    dforce = norm2(df)*particlemass

  end subroutine gravitational_force
!-----------------------------------------------------------------------
! Finding the zero of gravitational force by Newton method
! Due to the noise in the function this method is unsuccessful
!-----------------------------------------------------------------------
  subroutine newton_method(p, xyzh, particlemass, pNew, fNew, residual, npart, func)

    real,    intent(in)  :: p
    integer, intent(in)  :: npart
    real,    intent(in)  :: xyzh(4,npart)
    real,    intent(in)  :: particlemass
    real,    intent(out) :: fNew, residual

    real                 :: pNew, f, df

    interface
      subroutine func(point, force, dforce)
        real, intent (in)  :: point
        real, intent (out) :: force
        real, intent (out) :: dforce
      end subroutine func
    end interface

    ! compute function value evaluated at x
    call func(p, f, df)

    ! numerical second derivative
    ! write(*,*) p, f, df
    ! pNew = p + 1.e-4
    ! call func(pNew, fNew, df)
    ! df = (fNew - f)/(pNew - p)
    ! write(*,*) pNew, f, df
    ! stop

    ! Exit if f' is near or become zero
    if(abs(df) < 1.e-12) then
      print *, '[Error: newton_method] Function derivative becomes very close to zero or zero.'
      print *, 'f=',f, 'df/dp =',df
      print *, 'Aborting now in order to avoid division by zero.'
      stop
    end if

    ! Algorithm
    pNew = p - f/df
    fNew = f

    ! Search fails if a newly updated value x is out of the search domain
    ! if((pNew < pBeg) .or. (pNew > pEnd)) then
    !   print *, '[Error: newton_method] pNew',pNew, 'is out of domain.'
    !   print *, 'Failed in search. Aborting now.'
    !   stop
    ! end if

    ! Calculate a new residual
    residual = abs(pNew - p)

  end subroutine newton_method
!-----------------------------------------------------------------------
! The golden-section search is a technique for finding an extremum
!   (minimum or maximum) of a function inside a specified interval.
! The implementation is based on the more robust approach described in
!   V.G. Karmanov Mathematical programming, Moscow: FML, 2008, pp. 134-142.
!-----------------------------------------------------------------------
  real function golden_section_search_method(a, b, xyzh, particlemass, npart,&
    eps, err, extr, maxIter, Nest, iter, func)

    real,    intent(in)  :: a, b        ! left and right boundaries
    ! of the extremum search interval
    integer, intent(in)  :: npart
    real,    intent(in)  :: xyzh(4,npart)
    real,    intent(in)  :: particlemass
    real,    intent(in)  :: eps        ! specified accuracy

    real,    intent(out) :: err        ! achieved accuracy
    integer, intent(out) :: extr       ! extremum type: 1 - minimum; -1 - maximum

    integer, intent(in)  :: maxIter
    integer, intent(out) :: iter, Nest ! number of iterations actually made
    ! and their lower bound

    real                 :: q = 0.5d0*(sqrt(5.0d0)-1.0d0),&
      alpha, fy0, fz0,&
      a0, a1, b0, b1, y0,&
      y1, z0, z1, d0, d1, d2, d3, d10
    integer              :: nfail

    interface
      subroutine func(point, potential)
        real, intent (in)  :: point
        real, intent (out) :: potential
      end subroutine func
    end interface

    golden_section_search_method = 0.

    iter = 1
    err = 1.0d0
    nfail = 0
    Nest = int(log(eps/(b-a))/log(q))
    alpha = 0.8d0

    fy0 = 0.; fz0 = 0.; a0 = 0.; a1 = 0.; b0 = 0.; b1 = 0.
    y0 = 0.; y1 = 0.; z0 = 0.; z1 = 0.; d0 = 0.
    d1 = 0.; d2 = 0.; d3 = 0.; d10 = 0.

    call func(a, fy0)
    call func(a+eps, fz0)
    if(fy0 > fz0) then
      extr = 1         ! looking for a minimum
    else
      extr = -1        ! looking for a maximum
    end if
! step 1
    a0 = a
    b0 = b
! step 2
    do
      d10 = d1
      d0 = b0-a0
      d1 = q*d0
      d2 = d0-d1
! checking if precision is achieved (algorithm loops)
      if(abs(d1-d10) <= epsilon(1.0d0) .and. iter > maxIter) then
        nfail = nfail+1
! exit after two consecutive non-decreasing precision
        if(nfail >= 2) then
          err = d1
          golden_section_search_method = 0.5d0*(a1+b1)
          return
        end if
      else
        nfail = 0
      end if

      y0 = a0+d2
      z0 = b0-d2

      call func(y0, fy0)
      call func(z0, fz0)
! step 3
      do
        iter = iter+1
        d3 = d1-d2

        if(extr*fy0 <= extr*fz0) then
          a1 = a0
          b1 = z0
          z1 = y0
          y1 = a1+d3
          fz0 = fy0
          call func(y1, fy0)

          if(y1 >= z1) then
            a0 = a1
            b0 = b1
            exit   ! to step 2
          end if
        else
          a1 = y0
          b1 = b0
          y1 = z0
          z1 = b1-d3
          fy0 = fz0
          call func(z1, fz0)

          if(z1 <= y1) then
            a0 = a1
            b0 = b1
            exit   ! to step 2
          end if

        end if
! step 4
        if(d1 <= eps) then
          golden_section_search_method = 0.5d0*(a1+b1)
          err = d1
          return
        else

          if(d1 <= alpha*d0) then
            a0 = a1
            b0 = b1
            y0 = y1
            z0 = z1
            d1 = d2
            d2 = d3
            cycle ! to step 3
          else  ! d1 > eps .and. d1 > alpha**d0)
            a0 = a1
            b0 = b1
            exit ! to step 2
          end if

        end if

      end do

    end do

  end function golden_section_search_method
!-----------------------------------------------------------------------
  subroutine L1_point(method, xyzh, particlemass, npart, L1_proj, L1)

    integer, intent(in)  :: method ! 1 - Newton, 2 - Golden Section Search
    integer, intent(in)  :: npart
    real,    intent(in)  :: xyzh(4,npart)
    real,    intent(in)  :: particlemass
    real,    intent(out) :: L1_proj
    real,    intent(out) :: L1(3)

    integer              :: nIter
    integer, parameter   :: maxIter = 10
    real,    parameter   :: threshold = 1e-6

    ! for Newton's method
    real                 :: p, pNew, residual, f

    ! for Golden section search method
    real                 :: p1, p2
    integer              :: nest, extr
    real                 :: eps, err

    L1_proj = 0.
    L1 = 0.

    ! Initial values for residual and number of iteration
    residual = 1.e10
    nIter = 1

    p = 0.

    if(method == 1) then
      ! Keep search iteration until
      ! (a) residual is bigger then a user-defined threshold value, and
      ! (b) iteration number is less than a user-defined maximum iteration number.

      do while ((residual > threshold) .and. (nIter < maxIter))

        ! Search using conventional Newton's method
        call newton_method(p, xyzh, particlemass, pNew, f,&
          residual, npart, gravitational_force_wrapper)

        ! Save for the next search iteration
        p = pNew

        ! Update iteration number
        nIter = nIter + 1

        write(*,*) nIter, pNew, residual
      end do

      L1 = p*evector_old
      write(*,*) "L1 by Newton's method - ", L1, p

    else

      ! NB: find correct p1 and p2
      p1 = -20.
      p2 = 0.
      eps = threshold

      write(*,*) 'Golden section search method - Interval of extremum: p1=', p1, ' p2=', p2
      p = golden_section_search_method(p1, p2, xyzh, particlemass,&
        npart, eps, err, extr,&
        maxIter, Nest, nIter,&
        roche_potential_wrapper)

      L1 = p*evector_old

      if(extr == 1)  write(*,*) 'Minimum'
      if(extr == -1) write(*,*) 'Maximum'

      write(*,*) 'Iterations done: ', nIter, ', accuracу achieved is ', err
      write(*,*) 'L1 by Golden section search method - ', L1, p

    endif

    L1_proj = p

  end subroutine L1_point
!-----------------------------------------------------------------------

end module analysis
