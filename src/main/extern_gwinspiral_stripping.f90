!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!

! #define LAPACK

module extern_gwinspiral
!
! Simulates the inspiral of two stars in a circular orbit caused by gravitational wave
!   radiation.
!   Author: Bernard Field (supervisor: James Wurster & Paul Lasky)
!   Changes for stripping were done by Marat Potashov
!
! :References: e.g. Tong (2015) classical dynamics lecture notes
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - stop_ratio : *ratio of particles crossing CoM to indicate a merger*
!
! :Dependencies: centreofmass, dump_utils, infile_utils, io, physcon, units
!
  implicit none

  !
  ! Runtime parameters
  !
  real,       public :: stopratio = 0.005
  real,       public :: density_cutoff_cgs = 5.0d-5 ! The density threshold in code units

  integer,    public :: nstar(2) = 0 ! give default value in case dump header not read

  real,       public :: com(3)
  real,       public :: vcom(3)

  ! save from previous call (of get_gw_force)
  real, save, public :: evector_old(3)
  real, save, public :: time_old(2)
  real, save, public :: omega_old(3)

  real,       public :: fstar1(3), fstar2(3)
  real,       public :: fstar_tensor(3,3)
  real,       public :: comstar1(3), comstar2(3)
  real,       public :: vcomstar1(3), vcomstar2(3)

  !
  ! local variables
  !
  integer,   private :: n_threshhold
  real,      private :: m_density_cutoff
  real,      private :: fstar1_coef, fstar2_coef

  logical,   private :: isseparate = .true.

  !
  ! subroutines
  !
  public :: initialise_gwinspiral, gw_still_inspiralling
  public :: get_gw_force, get_gw_force_i
  ! public :: get_gwinspiral_vdependent_force, update_gwinspiral_leapfrog
  public :: get_momentofinertia, calculate_omega, correct_evector
  public :: read_options_gwinspiral, write_options_gwinspiral
  public :: read_headeropts_gwinspiral, write_headeropts_gwinspiral
  private

contains

!-----------------------------------------------------------------------
!+
!  Initialises external force to determine how many particle are in
!  each star
!+
!-----------------------------------------------------------------------
  subroutine initialise_gwinspiral(npart,nptmass,ierr)

    use units,        only: unit_density

    integer, intent(in)  :: npart,nptmass
    integer, intent(out) :: ierr
    integer              :: nerr

    !
    ! Calculate the number of particle required to 'merge' before the
    !  entire star is considered to be merged
    !
    n_threshhold = max(int(npart*stopratio),1)
    m_density_cutoff = density_cutoff_cgs / unit_density

    ierr = 0
    nerr = 0
    if(nstar(1) > 0) then
      write(*,"(2(a,i8))") ' Initialising inspiral two stars scenario, Nstar_1 = ',nstar(1),' Nstar_2 = ', nstar(2)
    elseif (nstar(1)==0 .and. nstar(2)==0 .and. nptmass==2) then
      write(*,"(2(a,i8))") ' Initialising inspiral on two sink particles'
    elseif (nstar(2)==0 .and. nptmass==1) then
      write(*,"(2(a,i8))") ' Initialising inspiral star-particle scenario, Nstar_1 = ',nstar(1)
    else
      ierr = 1
      isseparate = .false.
      return
    endif
    if(nerr > 0) ierr = 1

    evector_old = 0.d0
    time_old = 0.d0
    omega_old = 0.d0
    fstar_tensor = 0.d0

  end subroutine initialise_gwinspiral
!-----------------------------------------------------------------------
!+
!  Determines if the neutron stars are still inspirialing
!  (i.e. do we need to calculate the GW energy)
!+
!-----------------------------------------------------------------------
  subroutine gw_still_inspiralling(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,stopped_now)

    use physcon,      only: c
    use units,        only: unit_velocity
    use centreofmass, only: get_centreofmass

    integer, intent(in)  :: npart,nptmass
    real,    intent(in)  :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
    logical, intent(out) :: stopped_now

    integer              :: i,k1,k2
    real                 :: dir,dx,dy,dz
    real                 :: mstar1,mstar2
    real                 :: dirstar1(3)
    real                 :: c_code,sep

    !
    ! Calculate force coefficients
    !
    c_code      = c/unit_velocity
    fstar1_coef = 0.
    fstar2_coef = 0.
    stopped_now = .false.
    if( isseparate ) then
      if(nstar(1) == 0 .and. nptmass==2) then
        comstar1  = xyzmh_ptmass(1:3,2)
        mstar1    = xyzmh_ptmass(4,2)
        vcomstar1 = vxyz_ptmass(1:3,2)
      else
        call get_centreofmass(comstar1,vcomstar1,nstar(1),xyzh(:,1:nstar(1)),&
          vxyzu(:,1:nstar(1)),mass=mstar1)
      endif
      if(nstar(2) == 0 .and. nptmass>=1) then
        comstar2  = xyzmh_ptmass(1:3,1)
        mstar2    = xyzmh_ptmass(4,1)
        vcomstar2 = vxyz_ptmass(1:3,1)
      else
        call get_centreofmass(comstar2,vcomstar2,nstar(2),xyzh(:,nstar(1)+1:npart),&
          vxyzu(:,nstar(1)+1:npart),mass=mstar2)
      endif
      com = (comstar1*mstar1 + comstar2*mstar2)/(mstar1+mstar2)
      dirstar1 = comstar1 - com ! The directional vector to star 1
      vcom = (vcomstar1*mstar1 + vcomstar2*mstar2)/(mstar1+mstar2)
      !
      ! Determine how many particle are 'in' star 2 when they should be 'in' star 1
      ! as determined by their CoM location
      !
      k1 = 0
!$omp parallel default(none) &
!$omp shared(nstar,xyzh,com,dirstar1) &
!$omp private(i,dx,dy,dz,dir) &
!$omp reduction(+:k1)
!$omp do
      do i=1,nstar(1)
        dx  = xyzh(1,i) - com(1)
        dy  = xyzh(2,i) - com(2)
        dz  = xyzh(3,i) - com(3)
        dir = dirstar1(1)*dx + dirstar1(2)*dy + dirstar1(3)*dz
        if( dir < 0.0 ) k1 = k1 + 1
      enddo
!$omp enddo
!$omp end parallel
      !
      ! Determine how many particle are 'in' star 1 when they should be 'in' star 2
      ! as determined by their CoM location
      !
      k2 = 0
!$omp parallel default(none) &
!$omp shared(nstar,npart,xyzh,com,dirstar1) &
!$omp private(i,dx,dy,dz,dir) &
!$omp reduction(+:k2)
!$omp do
      do i=nstar(1)+1,nstar(1)+nstar(2)
        dx  = xyzh(1,i) - com(1)
        dy  = xyzh(2,i) - com(2)
        dz  = xyzh(3,i) - com(3)
        dir = dirstar1(1)*dx + dirstar1(2)*dy + dirstar1(3)*dz
        if( dir > 0.0 ) k2 = k2 + 1
      enddo
!$omp enddo
!$omp end parallel

      fstar1_coef = -32.0/5.0*mstar1*mstar2**3/c_code**5
      fstar2_coef = -32.0/5.0*mstar2*mstar1**3/c_code**5
      ! TODO:
      ! fstar1_coef = fstar1_coef*1.0e+5
      ! fstar2_coef = fstar2_coef*1.0e+5

      if(k1+k2 >= n_threshhold) then
        ! The stars have merged!
        isseparate = .false.    ! Stop emitting gravitational waves
        stopped_now = .true.    ! to trigger a print statement
      elseif (nptmass == 2) then
        ! stop the merger if two sinks are within each others accretion radii
        sep = sqrt(dot_product(comstar1 - comstar2,comstar1 - comstar2))
        print*,' SEP is ',sep,mstar1,mstar2,fstar1_coef,fstar2_coef
        if(sep <= xyzmh_ptmass(5,1) + xyzmh_ptmass(5,2)) then
          isseparate = .false.    ! Stop emitting gravitational waves
          stopped_now = .true.    ! to trigger a print statement
        endif
      endif
    endif

  end subroutine gw_still_inspiralling
!-----------------------------------------------------------------------
!+
!  Calculate the loss of energy (per star) from gravitational waves
!  This energy loss is in the form of a uniform reduction in force
!  Note: this is called immediately after gw_still_inspiralling, thus
!        the CoM's have just been calculated and stored
!+
!-----------------------------------------------------------------------
  subroutine get_gw_force(time,npart,xyzh,vxyzu,particlemass,nptmass,xyzmh_ptmass,vxyz_ptmass)

    use physcon,      only: c
    use units,        only: unit_velocity

    real,    intent(in) :: time
    integer, intent(in) :: npart, nptmass
    real,    intent(in) :: particlemass
    real,    intent(in) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)

    ! old version
    real                :: dx,dy,dz,separation,vstar1sq,vstar2sq
    ! new version
    integer             :: npartused
    real                :: rmax
    integer             :: smallIIndex
    real                :: inertia(3,3)
    real                :: principle(3), evectors(3,3)
    real                :: omega(3)
    real                :: d5q(3,3)
    real                :: c_code

    !
    ! Calculate force coefficients
    !

    ! the old version of drag force calculation

    if( isseparate ) then
      dx = comstar1(1) - comstar2(1)
      dy = comstar1(2) - comstar2(2)
      dz = comstar1(3) - comstar2(3)
      separation = sqrt(dx*dx + dy*dy + dz*dz)
      !
      ! Determine the speed of the two stars
      !
      vstar1sq = dot_product(vcomstar1,vcomstar1)
      vstar2sq = dot_product(vcomstar2,vcomstar2)
      !
      ! Compute the drag force vectors for each star
      !
      fstar1 = fstar1_coef * vcomstar1 / (vstar1sq*separation**5)
      fstar2 = fstar2_coef * vcomstar2 / (vstar2sq*separation**5)
    endif

    ! the new one version of gw drag force
    ! NB: only for two stars in stripping scenario

    c_code = c/unit_velocity
    if(nstar(1) > 0) then
      ! two stars in stripping scenario

      ! if(time > time_old(2)) then

      ! Calculate the inertia tensor with omega
      call get_momentofinertia(xyzh, vxyzu, com, vcom, npart, m_density_cutoff, particlemass,&
        npartused, inertia, principle, evectors, rmax, omega)

      smallIIndex = minloc(principle, dim=1)
      call correct_evector(evectors(:, smallIIndex), evector_old, time_old(2), time_old(1), omega_old)
      ! omega = calculate_omega(evectors(:, smallIIndex), evector_old, time_old(2), time_old(1), omega_old)
      evector_old = evectors(:, smallIIndex)

      time_old(1) = time_old(2)
      time_old(2) = time
      omega_old = omega

      call dquadrupole5(npart, xyzh, omega, particlemass, d5q)
      fstar_tensor = -2.d0/5.d0*d5q/c_code**5
      ! endif

    endif

  end subroutine get_gw_force
!-----------------------------------------------------------------------
!+
!  Calculate the loss of energy (per particle) from gravitational waves
!    i.e. determine if the particle is in star 1, or star 2, and use the
!    required force
!+
!-----------------------------------------------------------------------
  subroutine get_gw_force_i(i,xi,yi,zi,fextxi,fextyi,fextzi,phi)

    integer, intent(in)    :: i
    real,    intent(in)    :: xi,yi,zi
    real,    intent(inout) :: fextxi,fextyi,fextzi,phi

    fextxi = 0.0
    fextyi = 0.0
    fextzi = 0.0

    if(i <= 0) then
      return
    endif

    ! the old version of drag force calculation

    if(isseparate) then
      if(i <= nstar(1)) then
        fextxi = fstar1(1)
        fextyi = fstar1(2)
        fextzi = fstar1(3)
      elseif (nstar(2) > 0) then
        fextxi = fstar2(1)
        fextyi = fstar2(2)
        fextzi = fstar2(3)
      endif
    elseif (i == -1 .and. nstar(2)==0 .and. isseparate) then
      fextxi = fstar2(1) ! acceleration applied to sink particle 1 (star 2)
      fextyi = fstar2(2)
      fextzi = fstar2(3)
    elseif (i == -2 .and. nstar(1)==0 .and. isseparate) then
      fextxi = fstar1(1) ! acceleration applied to sink particle 2 (star 1)
      fextyi = fstar1(2)
      fextzi = fstar1(3)
    endif

    ! if(d1 <= 0.2) then
    ! write(*,*) 'f1 = ', i, fextxi, fextyi, fextzi
    ! endif

    ! the new one version of gw drag force
    ! NB: only for two stars in stripping scenario

    ! TODO: turn on new type of the forces
    if(nstar(1) > 0) then
      fextxi = fstar_tensor(1,1)*xi&
        + fstar_tensor(1,2)*yi&
        + fstar_tensor(1,3)*zi
      fextyi = fstar_tensor(2,1)*xi&
        + fstar_tensor(2,2)*yi&
        + fstar_tensor(2,3)*zi
      fextzi = fstar_tensor(3,1)*xi&
        + fstar_tensor(3,2)*yi&
        + fstar_tensor(3,3)*zi
    endif

    ! if(d1 <= 0.2) then
    ! write(*,*) 'f2 = ', i, fextxi, fextyi, fextzi
    ! endif

  end subroutine get_gw_force_i
!---------------------------------------------------------------
!+
!+
!---------------------------------------------------------------
  ! subroutine get_gwinspiral_vdependent_force(r,vel,bh_mass,vcrossomega)

  !   use vectorutils, only: cross_product3D

  !   real, intent(in)  :: r(3),vel(3)
  !   real, intent(in)  :: bh_mass
  !   real, intent(out) :: vcrossomega(3)

  !   vcrossomega = 0.d0

  ! end subroutine get_gwinspiral_vdependent_force
!---------------------------------------------------------------
!+
!+
!---------------------------------------------------------------
  ! subroutine update_gwinspiral_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,&
  !   vcrossomega,dt,xi,yi,zi,bh_mass)

  !   use vectorutils, only : cross_product3D,matrixinvert3D
  !   use io,          only : fatal,warning

  !   real, intent(in)    :: dt,xi,yi,zi,bh_mass
  !   real, intent(in)    :: vhalfx,vhalfy,vhalfz
  !   real, intent(inout) :: fxi,fyi,fzi
  !   real, intent(out)   :: vcrossomega(3)

  !   vcrossomega = 0.d0

  ! end subroutine update_gwinspiral_leapfrog
!-----------------------------------------------------------------------
!+
! Calculates the moment of inertia
! This is done about the coordinate axes whose origin is at the
!   centre of mass
! Mechanics, Third Edition: Volume 1 (Course of Theoretical Physics)
!   L. Landau, and E. Lifshitz. eq. 32.6
!+
!-----------------------------------------------------------------------
  subroutine get_momentofinertia(xyzh,vxyzu,center_of_mass,vcenter_of_mass,npart,density_cutoff,particlemass,&
    npartused,inertia,principle,evectors,rmax,omega)

    use part, only: rhoh
    use vectorutils, only: cross_product3D

    real,             intent(in)  :: xyzh(:,:)
    real,             intent(in)  :: vxyzu(:,:)
    real,             intent(in)  :: center_of_mass(3)
    real,             intent(in)  :: vcenter_of_mass(3)
    integer,          intent(in)  :: npart
    real,             intent(in)  :: density_cutoff
    real,             intent(in)  :: particlemass
    integer,          intent(out) :: npartused
    real,             intent(out) :: inertia(3,3)
    real,             intent(out) :: principle(3), evectors(3,3)
    real,             intent(out) :: rmax
    real,   optional, intent(out) :: omega(3)

    integer                       :: i
    real                          :: dot_inertia(3,3)
    integer                       :: smallIIndex
    real                          :: smallI
    real                          :: smallIEvector(3)
    real                          :: c(3)
    real                          :: c1(3)
    real                          :: dRdt(3)
! #ifdef LAPACK
    ! real                          :: inertia2(3,3)
! #endif
    real                          :: x,y,z,vx,vy,vz,r2,rmax2

    inertia     = 0.0
    dot_inertia = 0.0
    npartused   = 0
    rmax2       = 0.0
    if (present(omega))&
      omega       = 0.0

!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,center_of_mass,vcenter_of_mass,particlemass,density_cutoff) &
!$omp private(i,x,y,z,vx,vy,vz,r2) &
!$omp reduction(+:inertia,dot_inertia,npartused) &
!$omp reduction(max:rmax2)
!$omp do
    do i = 1, npart
      if(rhoh(xyzh(4,i),particlemass) > density_cutoff) then
        x = xyzh(1,i) - center_of_mass(1)
        y = xyzh(2,i) - center_of_mass(2)
        z = xyzh(3,i) - center_of_mass(3)
        vx = vxyzu(1,i) - vcenter_of_mass(1)
        vy = vxyzu(2,i) - vcenter_of_mass(2)
        vz = vxyzu(3,i) - vcenter_of_mass(3)
        inertia(1,1) = inertia(1,1) + y**2 + z**2
        inertia(2,2) = inertia(2,2) + x**2 + z**2
        inertia(3,3) = inertia(3,3) + x**2 + y**2
        inertia(1,2) = inertia(1,2) - x*y
        inertia(1,3) = inertia(1,3) - x*z
        inertia(2,3) = inertia(2,3) - y*z
        dot_inertia(1,1) = dot_inertia(1,1) + 2.0*y*vy + 2.0*z*vz
        dot_inertia(2,2) = dot_inertia(2,2) + 2.0*x*vx + 2.0*z*vz
        dot_inertia(3,3) = dot_inertia(3,3) + 2.0*x*vx + 2.0*y*vy
        dot_inertia(1,2) = dot_inertia(1,2) - vx*y - x*vy
        dot_inertia(1,3) = dot_inertia(1,3) - vx*z - x*vz
        dot_inertia(2,3) = dot_inertia(2,3) - vy*z - y*vz
        ! Additional useful values
        npartused    = npartused + 1
        r2           = x*x + y*y + z*z
        rmax2        = max(rmax2, r2)
      endif
    enddo
!$omp enddo
!$omp end parallel
    rmax = sqrt(rmax2)
    !--The symmetric components
    inertia(2,1) = inertia(1,2)
    inertia(3,1) = inertia(1,3)
    inertia(3,2) = inertia(2,3)
    dot_inertia(2,1) = dot_inertia(1,2)
    dot_inertia(3,1) = dot_inertia(1,3)
    dot_inertia(3,2) = dot_inertia(2,3)
    !--Multiply in constant
    inertia      = inertia*particlemass
    dot_inertia  = dot_inertia*particlemass
    !
! #ifdef LAPACK
    ! inertia2 = inertia
! #endif
    !
    !--Find the eigenvectors
    !
#ifndef LAPACK
    !  note: i is a dummy out-integer that we don't care about
    call jacobi(inertia,3,3,principle,evectors,i)

    ! write(*,*) 'Eigenvalues JACOBI:'
    ! do i = 1, 3
    !   write(*,*) i, principle(i)
    ! enddo
    ! write(*,*)
    ! write(*,*) 'Eigenvectors JACOBI:'
    ! do i = 1, 3
    !   write(*,*) i, evectors(:,i)
    ! enddo
    ! write(*,*)
#else
    call eigensystem(inertia,3,principle)
    evectors = inertia
    ! call eigensystem(inertia2,3,principle)
    ! evectors = inertia2

    ! write(*,*) 'Eigenvalues LAPACK:'
    ! do i = 1, 3
    !   write(*,*) i, principle(i)
    ! enddo
    ! write(*,*)
    ! write(*,*) 'Eigenvectors LAPACK:'
    ! do i = 1, 3
    !   write(*,*) i, evectors(:,i)
    ! enddo
    ! write(*,*)
#endif
    !
    if (present(omega)) then
      ! \[
      ! \mathbf{\Omega}_i^\mathrm{orb} =
      !   \sum_{j \neq i}\frac{1}{\lambda_j - \lambda_i}
      !     \left[\Big((\mathbf{e}_j^T \mathbf{\dot{I}} \mathbf{e}_i)\mathbf{e}_j\Big)
      !       \times \mathbf{e}_i \right]
      ! \]
      smallIIndex = minloc(principle, dim=1)
      smallIEvector = evectors(:, smallIIndex)
      smallI = principle(smallIIndex)
      c = matmul(dot_inertia, smallIEvector)
      dRdt = 0.0
      do i = 1, 3
        if(i == smallIIndex) cycle
        c1 = evectors(:, i)
        dRdt = dRdt + (dot_product(c,c1)*c1)/(principle(i) - smallI)
      enddo
      call cross_product3D(dRdt, smallIEvector, omega)
    endif
    !
  end subroutine get_momentofinertia
!-----------------------------------------------------------------------
!+
! LAPACK: DSYEV computes the eigenvalues and, optionally,
!   the left and/or right eigenvectors for SY matrices
! Calls the LAPACK diagonalization subroutine DSYEV
! input:  a(n,n) = real symmetric matrix to be diag
!         n  = size of a
! output: a(n,n) = orthonormal eigenvectors of a
!         v(n) = eigenvalues of a in ascending order
!+
!-----------------------------------------------------------------------
#ifdef LAPACK
  subroutine eigensystem(a,n,v)

    integer, intent(in)    :: n
    real,    intent(inout) :: a(n,n)
    real,    intent(out) :: v(n)

    integer :: lda
    real(kind=8) :: work(3*n-1)
    integer :: lwork
    integer :: info
    integer :: i

    info = 0
    lda = n
    lwork = 3*n-1
    call dsyev('V','U',n,a,lda,v,work,lwork,info)
    if(info < 0) then
      write(*,'(a, i3, a)') "INFO = ", info,&
        " the i-th argument had an illegal value"
    else if(info < 0) then
      write(*,'(a, i3, a)') "INFO = ", info,&
        " the algorithm failed to converge;&
      & i off-diagonal elements of an intermediate tridiagonal&
      & form did not converge to zero."
    endif

  end subroutine eigensystem
#endif
!-----------------------------------------------------------------------
!+
! Calculates the Jacobian
! Source: http://www.fing.edu.uy/if/cursos/fiscomp/extras/numrec/book/f11.pdf
!+
!-----------------------------------------------------------------------
  subroutine jacobi(a,n,np,d,v,nrot)

    integer, intent(in)    :: n,np
    integer, intent(out)   :: nrot
    real,    intent(inout) :: a(np,np)
    real,    intent(out)   :: d(np),v(np,np)
    integer, parameter     :: nmax = 500
!
! Computes all eigenvalues and eigenvectors of a real symmetric matrix, a,
!   which is of size n by n, stored in a physical np by np array.
! On output, elements of a above the diagonal are destroyed.
! d returns the eigenvalues of a in its first n elements.
! v is a matrix with the same logical and physical dimensions as a,
!   whose columns contain, on output, the normalized eigenvectors of a.
! nrot returns the number of Jacobi rotations that were required.
!
    integer :: i,ip,iq,j
    real :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

    do 12, ip=1,n  !Initialize  to  the  identity  matrix.
      do 11, iq=1,n
        v(ip,iq)=0.
11    enddo
      v(ip,ip)=1.
12  enddo
    do 13, ip=1,n
      b(ip)=a(ip,ip)
!Initialize b and d to the diagonal of a.
      d(ip)=b(ip)
      z(ip)=0.
!This  vector  will  accumulate  terms  of  the  form tapq as  in equation  (11.1.14).
13  enddo

    nrot=0
    do 24,i=1,50
      sm=0.
      do 15,ip=1,n-1
!Sum  off-diagonal elements.
        do 14,iq=ip+1,n
          sm=sm+abs(a(ip,iq))
14      enddo
15    enddo
      if(sm==0.)&
        return
!The normal return, which relies on quadratic convergence to machine  underflow.
      if(i < 4) then
        tresh=0.2*sm/n**2
!...on the first  three sweeps.
      else
        tresh=0.
!...thereafter.
      endif
      do 22,ip=1,n-1
        do 21,iq=ip+1,n
          g=100.*abs(a(ip,iq))
!After four sweeps, skip the rotation if the off-diagonal element is small.
          if((i > 4).and.(abs(d(ip))+g==abs(d(ip))).and.(abs(d(iq))+g==abs(d(iq)))) then
            a(ip,iq)=0.
          elseif (abs(a(ip,iq)) > tresh) then
            h=d(iq)-d(ip)
            if(abs(h)+g==abs(h)) then
              t=a(ip,iq)/h
!t=1/(2(theta))
            else
              theta=0.5*h/a(ip,iq)
!Equation  (11.1.10).
              t=1./(abs(theta)+sqrt(1.+theta**2))
              if(theta < 0.)t=-t
            endif
            c=1./sqrt(1+t**2)
            s=t*c
            tau=s/(1.+c)
            h=t*a(ip,iq)
            z(ip)=z(ip)-h
            z(iq)=z(iq)+h
            d(ip)=d(ip)-h
            d(iq)=d(iq)+h
            a(ip,iq)=0.
            do 16,j=1,ip-1
!Case of rotations 1<=j<p.
              g=a(j,ip)
              h=a(j,iq)
              a(j,ip)=g-s*(h+g*tau)
              a(j,iq)=h+s*(g-h*tau)
16          enddo
            do 17,j=ip+1,iq-1
!Case of rotations p<j<q.
              g=a(ip,j)
              h=a(j,iq)
              a(ip,j)=g-s*(h+g*tau)
              a(j,iq)=h+s*(g-h*tau)
17          enddo
            do 18,j=iq+1,n
!Case of rotations q<j<=n.
              g=a(ip,j)
              h=a(iq,j)
              a(ip,j)=g-s*(h+g*tau)
              a(iq,j)=h+s*(g-h*tau)
18          enddo
            do 19,j=1,n
              g=v(j,ip)
              h=v(j,iq)
              v(j,ip)=g-s*(h+g*tau)
              v(j,iq)=h+s*(g-h*tau)
19          enddo
            nrot=nrot+1
          endif
21      enddo
22    enddo
      do 23,ip=1,n
        b(ip)=b(ip)+z(ip)
        d(ip)=b(ip)
!Update d with the  sum of tapq,
        z(ip)=0.
!and  reinitialize z.
23    enddo
24  enddo
    return
  end subroutine jacobi
!-----------------------------------------------------------------------
! Function for finding omega vector from changes of the eigenvector
!-----------------------------------------------------------------------
  function calculate_omega(evector,evector_prev,time,time_prev,omega_prev) result(omega)

    use vectorutils, only: cross_product3D

    real, intent(inout) :: evector(3)
    real, intent(in)    :: evector_prev(3)
    real, intent(in)    :: time
    real, intent(in)    :: time_prev
    real, intent(in)    :: omega_prev(3)

    real                :: omega(3)
    real                :: dRdt(3)

    ! Calculate omega
    dRdt = 0.
    omega = 0.

    if(time > time_prev) then
      ! NB: change sign of evector if need it
      call correct_evector(evector,evector_prev,time,time_prev,omega_prev)
      dRdt = (evector - evector_prev)/(time - time_prev)
      call cross_product3D(evector,dRdt,omega)
      omega = omega/(norm2(evector)**2)
    endif

  end function calculate_omega
!-----------------------------------------------------------------------
  subroutine correct_evector(evector,evector_prev,time,time_prev,omega_prev)

    use vectorutils, only: cross_product3D

    real, intent(inout) :: evector(3)
    real, intent(in)    :: evector_prev(3)
    real, intent(in)    :: time
    real, intent(in)    :: time_prev
    real, intent(in)    :: omega_prev(3)

    real                :: omega(3)
    real                :: dRdt(3)

    ! Calculate omega
    dRdt = 0.
    omega = 0.

    if(time > time_prev) then
      dRdt = (evector - evector_prev)/(time - time_prev)
      call cross_product3D(evector,dRdt,omega)
      ! NB: change sign of evector if need it
      if(omega_prev(3)*omega(3) < 0.0)&
        evector = -1.*evector
    endif

  end subroutine correct_evector
!-----------------------------------------------------------------------
! Calculate the fifth time derivative of the quadrupole moment
!-----------------------------------------------------------------------
  subroutine dquadrupole5(npart,xyzh,omega,particlemass,d5q)

    integer, intent(in)  :: npart
    real,    intent(in)  :: xyzh(:,:)
    real,    intent(in)  :: omega(3)
    real,    intent(in)  :: particlemass
    real,    intent(out) :: d5q(3,3)

    real                 :: omegasq
    real                 :: omegari
    integer              :: i, ia, ib, ii, ik
    real                 :: coeff(3)

    d5q = 0.d0

    omegasq = dot_product(omega, omega)

!$omp parallel default(none) &
!$omp shared(npart,xyzh,omega,omegasq) &
!$omp private(i,omegari,ii,ia,ib,coeff,ik) &
!$omp reduction(+:d5q)
!$omp do
    do i = 1, npart
      omegari = dot_product(omega, xyzh(1:3,i))

      coeff = 0.d0
      do ii = 1, 3
        do ia = 1, 3
          do ib = 1, 3
            coeff(ii) = coeff(ii) + levicivita(ii,ia,ib)*omega(ia)*xyzh(ib,i)
          enddo
        enddo
      enddo

      do ii = 1, 3
        do ik = 1, 3
          d5q(ii,ik) = d5q(ii,ik) +&
            (16.d0*xyzh(ik,i)*omegasq -&
            15.d0*omega(ik)*omegari)*coeff(ii) +&
            (16.d0*xyzh(ii,i)*omegasq -&
            15.d0*omega(ii)*omegari)*coeff(ik)
        enddo
      enddo
    enddo
!$omp enddo
!$omp end parallel

    d5q = d5q*omegasq*particlemass

  end subroutine dquadrupole5
!-----------------------------------------------------------------------
! Returns the Levi-Civita-Symbol (permutation symbol)
!-----------------------------------------------------------------------
  pure function levicivita(i,j,k) result(lc)

    integer, intent(in)           :: i, j, k
    real(8)                       :: lc

    lc = 0.5d0 * (i - j) * (j - k) * (k - i)

  end function levicivita
!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
  subroutine write_options_gwinspiral(iunit)

    use infile_utils, only: write_inopt

    integer, intent(in) :: iunit

    call write_inopt(stopratio,'stop_ratio','ratio of particles crossing CoM to indicate a merger',iunit)

  end subroutine write_options_gwinspiral
!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
  subroutine read_options_gwinspiral(name,valstring,imatch,igotall,ierr)

    use io, only: fatal

    character(len=*), intent(in)  :: name,valstring
    logical,          intent(out) :: imatch,igotall
    integer,          intent(out) :: ierr

    integer, save                 :: ngot = 0
    character(len=30), parameter  :: where = 'read_options_gwinspiral'

    imatch  = .true.
    igotall = .false.
    select case(trim(name))
     case('stop_ratio')
      read(valstring,*,iostat=ierr) stopratio
      if(stopratio < 0.)  call fatal(where,'Cannot have negative merger percentage of particle overlap')
      if(stopratio > 1.)  call fatal(where,'Cannot have ratio of particle overlap > 1')
      ngot = ngot + 1
    end select

    igotall = (ngot >= 1)

  end subroutine read_options_gwinspiral
!-----------------------------------------------------------------------
!+
!  writes relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
  subroutine write_headeropts_gwinspiral(hdr,ierr)

    use dump_utils, only: dump_h,add_to_header

    type(dump_h), intent(inout) :: hdr
    integer,      intent(out)   :: ierr

    ierr = 0
    call add_to_header(nstar(1),'Nstar_1',hdr,ierr)
    call add_to_header(nstar(2),'Nstar_2',hdr,ierr)

  end subroutine write_headeropts_gwinspiral
!-----------------------------------------------------------------------
!+
!  reads relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
  subroutine read_headeropts_gwinspiral(hdr,nptmass,ierr)

    use dump_utils, only: dump_h,extract

    type(dump_h), intent(in)  :: hdr
    integer,      intent(in)  :: nptmass
    integer,      intent(out) :: ierr
    integer :: ierr1,ierr2

    ierr  = 0
    call extract('Nstar_1',nstar(1),hdr,ierr1)
    call extract('Nstar_2',nstar(2),hdr,ierr2)

    if(ierr1 /= 0 .or. ierr2 /= 0) then
      ! if there are two sink particles and Nstar_1
      if(nptmass >= 2 .and. ierr1 /= 0 .and. ierr2 /= 0) then
        write(*,*) ' WARNING: NStar_1 and NStar_2 not found: applying GW emission to sink particles'
      elseif (nptmass > 1 .and. ierr2 /= 0) then
        write(*,*) ' WARNING: NStar_2 not found: applying GW emission to 1 sink particle + gas'
      else
        write(*,*) ' ERROR: NStar_1 and NStar_2 not present in dump file'
        ierr = 1
      endif
    endif

  end subroutine read_headeropts_gwinspiral

end module extern_gwinspiral
