!This is the complete program to compute the potential anywhere outside the cube. 
program potential_parallel
    implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include <mpif.h>
#include <sprng_f.h>

    !Parameters for the diffusion 
    real*8 :: x, y, z
    real*8, parameter :: side = 1.0d0
    real*8, parameter :: radius = sqrt(3.0d0) / 2.0d0
    logical :: inf, inside_cube
    integer*8, parameter :: N = 10**5               !Number of constituent particles of the cube
    integer*8, parameter :: num_wos = 10**6         !Number of diffused particles 
    integer*8 :: particles_cube = 0
    integer*8 :: particles_inf = 0
    real*8 :: sigma, theta

    real*8, parameter :: lp_radius = 0.5d0
    real*8 :: x0, y0, z0

    !Parameters for the grid 
    integer*8 :: p, q 
    real*8, parameter :: a = -0.5d0, b = 0.5d0
    integer*8, parameter :: grid_divisions = CEILING(sqrt(REAL(N) / 6.0d0))
    real*8 :: spacing = (b - a) / grid_divisions

    !Arrays to store the charge density values from the 6 faces
    real*8, dimension(grid_divisions, grid_divisions) :: sigma_array             !TOP

    !Parameters for the potential calculation
    real*8 :: V, d 

    !The point where we want to calculate the potential 
    real*8, parameter :: x1 = 0.300d0
    real*8, parameter :: y1 = 0.300d0
    real*8, parameter :: z1 = 0.5d0

    integer :: ierr, rank, processN
    integer :: counti, countf, count_rate, count_max
    real(8) :: runtime
 
    integer :: seed, gentype
    SPRNG_POINTER :: stream
    common /rnd/ stream 

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, processN, ierr)

    seed = 2131257613
    gentype = 5
    call init_sprng(seed, SPRNG_DEFAULT, gentype, stream)

    call system_clock(counti, count_rate, count_max)

    real*8, dimension(:,:), allocatable :: sigma_local
    allocate(sigma_local(grid_divisions, grid_divisions))

    !Distribute workload among processes
    do p = rank + 1, grid_divisions, processN
        do q = 1, grid_divisions
            x0 = a + spacing * p
            y0 = a + spacing * q
            z0 = 0.5d0

            sigma = 0.0d0
            do i = 1, num_wos
                inf = .false.
                inside_cube = .false.

                call center(x0, y0, z0, x, y, z, lp_radius)
                call unif_samp(x, y, z, lp_radius, theta)

                do while ((inside_cube .eqv. .false.) .and. (inf .eqv. .false.))
                    call wop(x, y, z)
                    call in_the_cube(inside_cube, x, y, z)
                    call infinity(inf, x, y, z)

                    if (inside_cube) then
                        particles_cube = particles_cube + 1
                        exit
                    else if (inf) then
                        particles_inf = particles_inf + 1
                        call summ(theta, sigma, inf)
                        exit
                    else
                        call center(x0, y0, z0, x, y, z, lp_radius)
                        call bts(inf, inside_cube, x, y, z)
                    end if
                end do     

                sigma = sigma / real(num_wos)
                sigma_local(p, q) = sigma
            end do
        end do
    end do

    !Gather sigma_local from all processes to sigma_array
    call MPI_GATHER(sigma_local, grid_divisions * grid_divisions, MPI_DOUBLE_PRECISION, &
                    sigma_array, grid_divisions * grid_divisions, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !POTENTIAL CALCULATION
    V = 0.0d0

    !Distribute workload among processes for potential calculation
    do p = rank + 1, grid_divisions, processN
        do q = 1, grid_divisions
            x0 = a + spacing * p
            y0 = a + spacing * q
            z0 = 0.5d0

            call distance(x0, y0, z0, x1, y1, z1, d)
            V = V + sigma_array(p, q) / d
        end do
    end do

    !Sum up local potentials across all processes
    call MPI_REDUCE(V, V, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        print *, "Potential: ", V
    end if

    deallocate(sigma_local)

    call system_clock(countf, count_rate, count_max)

    write(*, "('Location: ', 3F10.3)") x1, y1, z1
    runtime = real(countf - counti) / real(count_rate)
    write(*, '(A, F6.3, A)') "Runtime: ", runtime, "s"

    call MPI_FINALIZE(ierr)

end program potential_parallel


!Coordinate of the center of the last passage sphere about any face.  
subroutine center(x,y,z,x_center,y_center,z_center,lp_radius)
    implicit none 
    real*8, intent(in) :: x,y,z,lp_radius
    real*8, intent(inout) :: x_center, y_center, z_center
    real*8, parameter :: side = 0.5d0

    if (z >= side .and. abs(x) <= side .and. abs(y) <= side) then
        !Upper x-y plane
        x_center = x
        y_center = y
        z_center = z+lp_radius
    else if (z <= -side .and. abs(x) <= side .and. abs(y) <= side) then
        !Lower x-y plane
        x_center = x
        y_center = y
        z_center = z-lp_radius
    else if (y >= side .and. abs(x) <= side .and. abs(z) <= side) then
        !Back x-z plane
        x_center = x
        y_center = y+lp_radius
        z_center = z
    else if (y <= -side .and. abs(x) <= side .and. abs(z) <= side) then
        !Front x-z plane
        x_center = x
        y_center = y-lp_radius
        z_center = z
    else if (x >= side .and. abs(y) <= side .and. abs(z) <= side) then
        !Right y-z plane
        x_center = x+lp_radius
        y_center = y
        z_center = z
    else
        !Left y-z plane
        x_center = x-lp_radius
        y_center = y
        z_center = z
    end if

end subroutine center 

!Uniform sample of points in the sphere
subroutine unif_samp(x, y, z, radius_inf,theta)
    implicit none
    real*8, intent(in):: radius_inf
    real*8, intent(inout):: x, y, z
    real*8:: urnd1, urnd2, phi, R
    real*8, parameter:: pi = dacos(-1.d0)
    real*8, intent(inout) :: theta

#define SIMPLE_SPRNG
#include <sprng_f.h>

    SPRNG_POINTER :: stream
    common /rnd/ stream 

    urnd1 = sprng(stream)
    urnd2 = sprng(stream)

    R = radius_inf
    theta = dacos(1.d0-2.d0 * urnd1)
    phi = 2.d0 * pi * urnd2

    x = x + R * dsin(theta) * dcos(phi)  
    y = x + R * dsin(theta) * dsin(phi) 
    z = z + R * dcos(theta) 
end subroutine unif_samp

!Distance from the point in the sphere to the nearest plane 
subroutine d_to_plane(x, y, z, d)
    implicit none
    real*8, intent(in) :: x, y, z
    real*8, intent(out) :: d
    real*8, parameter :: side = 0.5d0 

    d = max(abs(x),abs(y),abs(z))
    d = d-side   
end subroutine d_to_plane

!Distance from any two arbitrary points 
subroutine distance(x, y, z, x0, y0, z0, d)
    implicit none
    real*8, intent(in):: x, y, z, x0, y0, z0
    real*8, intent(out) :: d

    d = dsqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
end subroutine distance

!WOP
subroutine wop(x, y, z)
    implicit none
    real*8, intent(inout) :: x, y, z
    real*8 :: xxx, r, yyy, phi1, d
    real*8, parameter :: pi = dacos(-1.d0)
    double precision :: a,b

#define SIMPLE_SPRNG
#include <sprng_f.h>
    
    SPRNG_POINTER :: stream
    common /rnd/ stream 
  
    a = 0.5d0
    b = a * sqrt(3.0d0)
  
    xxx = sprng(stream)
    r = sqrt((1.d0 - xxx**2) / xxx**2) 
    yyy = sprng(stream)
    phi1 = 2.d0 * pi * yyy
  
    if ((abs(x) >= abs(y)) .and. (abs(x) >= abs(z))) then
      d = abs(x) - a
      if (x > 0) then
        x = a
        y = y + r * d * sin(phi1)
        z = z + r * d * cos(phi1)
      else
        x = -a
        y = y + r * d * sin(phi1)
        z = z + r * d * cos(phi1)
      endif
    elseif ((abs(y) >= abs(x)) .and. (abs(y) >= abs(z))) then
      d = abs(y) - a
      if (y > 0) then
        y = a
        x = x + r * d * sin(phi1)
        z = z + r * d * cos(phi1)
      else
        y = -a
        x = y + r * d * sin(phi1)
        z = z + r * d * cos(phi1)
      endif
    elseif ((abs(z) >= abs(x)) .and. (abs(z) >= abs(y))) then
      d = abs(z) - a
      if (z > 0) then
        z = a
        x = x + r * d * sin(phi1)
        y = y + r * d * cos(phi1)
      else
        z = -a
        x = x + r * d * sin(phi1)
        y = y + r * d * cos(phi1)
      endif
    endif
end subroutine wop

!Back to sphere 
subroutine bts(x,y,z)
    implicit none
  
    real*8, intent(inout) :: x, y, z
    real*8 :: r
    real*8, parameter :: pi = dacos(-1.d0)
    real*8 :: t, p2, costh, sinth, phi, cosph, sinph, xold, yold, zold, p, theta
    double precision :: a,b

#define SIMPLE_SPRNG
#include <sprng_f.h>
    
    SPRNG_POINTER :: stream
    common /rnd/ stream 
  
    a = 0.5d0
    b = a * sqrt(3.0d0)
  
    call distance(0.0d0,0.0d0,0.0d0,x,y,z,r)
  
    t = b / r 
    p2=sprng(stream)
  
    costh = -(1.d0 - t)**2 + 2.d0 * (1.d0 - t) * (1.d0 + t**2) * p2 + 2.d0 * t * (1.d0 + t**2) * p2**2
    costh = costh / (1.d0 - t + 2.d0 * t * p2)**2
    theta = dacos(costh)
    sinth = dsqrt(dabs(1.d0 - costh**2))
  
    phi=sprng(stream)
    phi = 2.d0 * pi * phi
    cosph = dcos(phi)
    sinph = dsin(phi)
  
    xold = x
    yold = y
    zold = z

    p = dsqrt(xold**2+yold**2) 
    x = b * (sinth * cosph * xold * zold / p / r - sinth * sinph * yold / p + costh * xold / r)
    y = b * (sinth * cosph * yold * zold / p / r + sinth * sinph * xold / p + costh * yold / r)
    z = b * (-sinth * cosph * p / r + costh * zold / r)

end subroutine bts 

!Infinity Check
subroutine infinity(inf,x,y,z)
    implicit none
    real*8, intent(in) :: x,y,z
    logical, intent(inout) :: inf
    real*8 :: ff
    double precision :: a,b,r

#define SIMPLE_SPRNG
#include <sprng_f.h>

    SPRNG_POINTER :: stream
    common /rnd/ stream 

    inf = .false. 
    a = 0.5d0
    b = a * sqrt(3.0d0)
  
    call distance(0.0d0, 0.0d0, 0.0d0, x,y,z,r)
  
    
    ff = sprng(stream)
    if (ff * r >= b) then
      inf = .true.
    end if

end subroutine infinity


!Inside the cube 
subroutine in_the_cube(inside_cube, x, y, z) 
    implicit none
    logical, intent(inout) :: inside_cube
    real*8, intent(in) :: x, y, z

    if ((abs(x) <= 0.5d0) .and. (abs(y) <= 0.5d0) .and. (abs(z) <= 0.5d0)) then
        inside_cube = .true.
    else
        inside_cube = .false.
    end if
end subroutine in_the_cube

!Integration
subroutine summ(theta,sigma,inf)
    implicit none 
    real*8, intent(in) :: theta
    logical, intent(in) :: inf
    real*8, intent(inout) :: sigma
    real*8 :: const 

    const = 16.d0 * dsqrt(2.d0)

    if (inf .eqv. .true.) then 
      sigma = sigma + dsin(theta)/(const*(1.d0+dcos(theta))**1.5)
    end if 

end subroutine summ

