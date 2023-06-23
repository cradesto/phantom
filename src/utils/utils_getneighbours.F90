!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module getneighbours
!
! A set of routines generate neighbour lists for all particles, and
!  read/write the list to file
!
! :References:
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, kdtree, kernel, linklist, part
!
  implicit none

  public :: generate_neighbour_lists, write_neighbours
  integer, public, allocatable, dimension(:)   :: sph_neighbours_count
  integer, public, allocatable, dimension(:,:) :: sph_neighbours
  integer, public, allocatable, dimension(:)   :: particles_in_cells_count
  integer, public, allocatable, dimension(:)   :: cells
  integer, public            :: min_neighbours_count, max_neighbours_count
  real,    public            :: mean_neighbours_count, sd_neighbours_count
  integer, public            :: min_particles_in_cells_count, max_particles_in_cells_count
  real,    public            :: mean_particles_in_cells_count, sd_particles_in_cells_count
  logical                    :: neigh_overload
  integer,         parameter :: maxcellcache =  50000
  integer,         parameter :: neighall     = 100000 ! maximum number of neighbours to test
  integer, public, parameter :: neighmax     =   2000 ! maximum number of neighbours to store

  private

contains

!-----------------------------------------------------------------------
!+
! Generate neighbour lists for all particles
!+
!-----------------------------------------------------------------------
  subroutine generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile,write_neighbour_list)

    use dim,      only:maxneigh,maxp,ncellsmax
    use kernel,   only:radkern2
    use linklist, only:ncells, ifirstincell, set_linklist, get_neighbour_list
    use part,     only:get_partinfo, igas, maxphase, iphase, iamboundary, iamtype
    use kdtree,   only:inodeparts,inoderange
#ifdef PERIODIC
    use boundary, only:dxbound,dybound,dzbound
#endif
    real,             intent(in)     :: xyzh(:,:),vxyzu(:,:)
    integer,          intent(in)     :: npart
    character(len=*), intent(in)     :: dumpfile
    logical,          intent(in)     :: write_neighbour_list
    real,allocatable, dimension(:,:) :: dumxyzh

    integer      :: i,j,k,p,ip,icell,ineigh,nneigh,dummynpart
    integer      :: ineigh_all(neighall)
    real         :: dx,dy,dz,rij2
    real         :: hi1,hj1,hi21,hj21,q2i,q2j
    integer,save :: listneigh(maxneigh)
    real,   save :: xyzcache(maxcellcache,4)
    real         :: rneigh_all(neighall)
    !$omp threadprivate(xyzcache,listneigh)
    character(len=100) :: neighbourfile

    !****************************************
    ! 1. Build kdtree and linklist
    ! --> global (shared) neighbour lists for all particles in tree cell
    !****************************************

    print*, 'Building kdtree and linklist: '
    allocate(dumxyzh(4,npart))
    dumxyzh = xyzh
    dummynpart = npart
    call set_linklist(dummynpart,npart,dumxyzh,vxyzu(:,1:npart))

    print*, 'Allocating arrays for neighbour storage : '
    if(.not. allocated(sph_neighbours_count))&
      allocate(sph_neighbours_count(npart))
    if(.not. allocated(sph_neighbours))&
      allocate(sph_neighbours(npart,neighmax))
    if(.not. allocated(particles_in_cells_count))&
      allocate(particles_in_cells_count(ncellsmax))
    if(.not. allocated(cells))&
      allocate(cells(npart))

    sph_neighbours_count(:)     = 0
    sph_neighbours(:,:)         = 0
    particles_in_cells_count(:) = 0
    cells(:)                    = 0

    print "(A,I5)", 'Maximum neighbour number allocated:  ', neighmax

    !***************************************
    ! 2. Assign neighbour lists to particles by searching shared list of host cell
    !***************************************

    print*, 'Creating neighbour lists for particles'

    !$omp parallel default(none) &
    !$omp shared(ncells,ifirstincell,npart,maxphase,maxp,inodeparts,inoderange) &
    !$omp shared(xyzh,vxyzu,iphase,sph_neighbours_count,sph_neighbours,cells) &
#ifdef PERIODIC
    !$omp shared(dxbound,dybound,dzbound) &
#endif
    !$omp private(icell,i,j,k,p,ip)&
    !$omp private(nneigh,ineigh_all,rneigh_all) &
    !$omp private(hi1,hi21,hj1,hj21,rij2,q2i,q2j) &
    !$omp private(dx,dy,dz)
    !$omp do schedule(runtime)
    over_cells: do icell=1,int(ncells)

      k = ifirstincell(icell)

      ! Skip empty/inactive cells
      if (k <= 0) cycle over_cells

      ! Get neighbour list for the cell
      call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)

      particles_in_cells_count(icell) = inoderange(2,icell) - inoderange(1,icell) + 1

      ! Loop over particles in the cell
      over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
        i = inodeparts(ip)
        cells(i) = icell
        if (maxphase==maxp) then
          if (iamboundary( iamtype(iphase(i)) )) cycle over_parts
        endif

        ! Fill neighbour list for this particle
        sph_neighbours_count(i) = 0
        ineigh_all    = 0
        rneigh_all    = 0.
        hi1  = 1.0/xyzh(4,i)
        hi21 = hi1*hi1

        over_neighbours: do ineigh = 1,nneigh
          j = abs(listneigh(ineigh))

          ! Skip self
          if (i==j) cycle over_neighbours

          dx = xyzh(1,i) - xyzh(1,j)
          dy = xyzh(2,i) - xyzh(2,j)
          dz = xyzh(3,i) - xyzh(3,j)
#ifdef PERIODIC
          if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
          if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
          if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
          rij2 = dx*dx + dy*dy + dz*dz
          q2i  = rij2*hi21

          hj1  = 1.0/xyzh(4,j)
          hj21 = hj1*hj1
          q2j  = rij2*hj21

          is_sph_neighbour: if (q2i < radkern2 .or. q2j < radkern2) then
            sph_neighbours_count(i) = sph_neighbours_count(i) + 1
            if (sph_neighbours_count(i) <= neighall) then
              ineigh_all(sph_neighbours_count(i)) = j
              rneigh_all(sph_neighbours_count(i)) = rij2
            else
              print*, 'neighbour finding.  neighcount > neighall.  aborting'
              stop
            endif
          endif is_sph_neighbour
        enddo over_neighbours ! End loop over neighbours
        ! Failsafe if too many neighbours
        if (sph_neighbours_count(i) <= neighmax) then
          sph_neighbours(i,1:sph_neighbours_count(i)) = ineigh_all(1:sph_neighbours_count(i))
        else
          print*, 'Neighbour finding: There are ',sph_neighbours_count(i),' neighbours for i = ',i
          print*, 'Neighbour finding: Keeping the ',neighmax,' closest neighbours.'
          sph_neighbours_count(i) = neighmax
          do p = 1,neighmax
            j = minloc(rneigh_all(1:sph_neighbours_count(i)),1)
            sph_neighbours(i,p) = ineigh_all(j)
            rneigh_all(j) = huge(rij2)
          enddo
        endif
      enddo over_parts         ! End loop over particles in the cell
    enddo over_cells           ! End loop over cells in the kd-tree
    !$omp enddo
    !$omp end parallel

    ! Do some simple stats on neighbour numbers
    print*, 'Statistics for neighbours'
    min_neighbours_count  = 0
    max_neighbours_count  = 0
    mean_neighbours_count = 0.0
    sd_neighbours_count   = 0.0
    call stats(npart, sph_neighbours_count, min_neighbours_count, max_neighbours_count,&
      mean_neighbours_count, sd_neighbours_count)

    print*, 'Statistics for cells'
    min_particles_in_cells_count  = 0
    max_particles_in_cells_count  = 0
    mean_particles_in_cells_count = 0.0
    sd_particles_in_cells_count   = 0.0
    i = count(particles_in_cells_count > 0.0)
    call stats(i, particles_in_cells_count, min_particles_in_cells_count, max_particles_in_cells_count,&
      mean_particles_in_cells_count, sd_particles_in_cells_count)

    !**************************************
    ! 3. Output neighbour lists to file (if requested; these files can become very big)
    !**************************************
    if (write_neighbour_list) then
      neighbourfile = 'neigh_'//TRIM(dumpfile)
      call write_neighbours(neighbourfile, xyzh, vxyzu, npart)
      print*, 'Neighbour finding complete for file ', TRIM(dumpfile)
    endif

    deallocate(dumxyzh)

  end subroutine generate_neighbour_lists
!-----------------------------------------------------------------------
!+
! Calculates the mean and standard deviation of the array number
! Also calculates a 5 sigma deviation from mean
! (This is principally used as a diagnostic aid for structure finding
! algorithms that rely on the nearest neighbours, like CLUMPFIND)
!+
!-----------------------------------------------------------------------
  subroutine stats(n, array, minimum, maximum, mean, sd)

    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: array
    integer, intent(out) :: minimum, maximum
    real, intent(out) :: mean, sd

    integer             :: i, s

    ! Calculate mean and standard deviation of array

    maximum = maxval(array)
    minimum = minval(array)

    print*, 'The maximum count is ', maximum
    print*, 'The minimum count is ', minimum

    if (maximum > neighmax) then
      print*, 'WARNING! Count too large for allocated arrays.  aborting'
      stop
    endif

    mean = sum(array)/REAL(n)
    sd   = 0.0

    s = size(array)
!$omp parallel default(none) &
!$omp shared(array,mean,s) &
!$omp private(i) &
!$omp reduction(+:sd)
!$omp do schedule(runtime)
    do i = 1, s
      sd = sd + (array(i) - mean)**2
    enddo
!$omp enddo
!$omp end parallel

    sd = sqrt(sd/REAL(n))

    print*, 'Mean number is ', mean
    print*, 'Standard Deviation: ', sd

  end subroutine stats
!-----------------------------------------------------------------------
!+
! Writes neighbour data to binary file
!+
!-----------------------------------------------------------------------
  subroutine write_neighbours(neighbourfile,xyzh,vxyzu,npart)

    use linklist, only: get_cell_location
    use kernel,   only: radkern
    use linklist, only: ncells

    real, intent(in)    :: xyzh(:,:),vxyzu(:,:)
    integer, intent(in) :: npart
    integer             :: i,j,icell
    character(len=100)  :: neighbourfile
    real                :: com(3), size, hmax
    real, parameter     :: tolerance = 2.0e0  ! A dummy parameter used to keep file format similar to other codes (Probably delete later)

    neigh_overload = .false.
    neighbourfile  = TRIM(neighbourfile)
    print*, 'Writing neighbours to file ', neighbourfile

    open(2,file='data.txt')
    do i=1,npart
      icell = cells(i)
      write(2,'(i7)', advance="no") icell
    enddo
    write(2,'()')

    open(2, file=neighbourfile)!, form='unformatted')
    write(2,*) neighmax, tolerance, mean_neighbours_count, sd_neighbours_count
    write(2,*) (sph_neighbours_count(i), i=1,npart)

    write(2,*) mean_particles_in_cells_count, sd_particles_in_cells_count
    do icell=1,ncells
      if(particles_in_cells_count(icell) > 0) then
        write(2,'(i12)', advance="no") icell
      endif
    enddo
    write(2,'()')

    do icell=1,ncells
      if(particles_in_cells_count(icell) > 0) then
        write(2,'(i12)', advance="no") particles_in_cells_count(icell)
      endif
    enddo
    write(2,'()')

    do i=1,npart
      icell = cells(i)
      call get_cell_location(icell, com, size, hmax)
      if (sph_neighbours_count(i) > neighmax) then
        neigh_overload = .true.
        print*, 'neighbour overload: ', sph_neighbours_count(i), neighmax
        write(2,*) xyzh(4,i)*radkern, icell, com, size, hmax, (sph_neighbours(i,j), j=1,neighmax)
      else
        write(2,*) xyzh(4,i)*radkern, icell, com, size, hmax, (sph_neighbours(i,j), j=1,sph_neighbours_count(i))
      endif
    enddo
    close(2)

    if (neigh_overload) then
      print*, 'WARNING! File write incomplete: neighbour count exceeds array size'
    else
      print*, 'File Write Complete'
    endif

  end subroutine write_neighbours
!-----------------------------------------------------------------------
end module getneighbours
