
PROGRAM wolff
      IMPLICIT NONE
      INTEGER, PARAMETER:: length = 10, latsize = length * length
      INTEGER, PARAMETER:: eq = 100000, steps = 100, tsteps = 400
      REAL*8, PARAMETER:: tsize = 0.01
      INTEGER:: lattice(latsize)
      INTEGER:: i,j, mag
      REAL*8:: tau, beta,  magav, magsum, mag2sum, mag2av
      REAL*8:: ensum, enav, en2sum, en2av

      OPEN (UNIT=20,FILE="averages.dat",ACTION="write",STATUS="replace")

      !Initialize the seed
      CALL init_random_seed()

      DO j = 1, tsteps
        !Initialize variables
        tau = tsize * DFLOAT(j)
        beta = 1.0d0 / tau
        magsum = 0
        mag2sum = 0
        ensum = 0
        en2sum = 0

        !Initialize the lattice
        CALL init_lattice(lattice)


        CALL equilibrate(lattice)

        !Initialize magnetization measurements
        magsum = magsum + ABS(SUM(lattice))
        mag2sum = mag2sum + (ABS(SUM(lattice)) * ABS(SUM(lattice)))
        ensum = ensum + energy(lattice)
        en2sum = en2sum + energy(lattice) * energy(lattice)

        !Do the wolff algorithm and take measurements for the given number of
        !steps
        DO i = 1, steps
          !Flip a cluster
          CALL flip_cluster(lattice)

          !Measure magnetization
          magsum = magsum + ABS(SUM(lattice))
          mag2sum = mag2sum + (ABS(SUM(lattice)) * ABS(SUM(lattice)))
          ensum = ensum + energy(lattice)
          en2sum = en2sum + energy(lattice) * energy(lattice)

          END DO

        magav =  (magsum / DFLOAT(steps + 1)) /DFLOAT(latsize)
        mag2av =  (mag2sum / DFLOAT(steps + 1)) /DFLOAT(latsize)
        enav = ensum / DFLOAT(steps + 1)
        en2av = en2sum / DFLOAT(steps + 1)

        WRITE(20,*) tau, magav, mag2av, enav, en2av
        END DO
        CLOSE(20)

      CONTAINS
! initialize the seed for random numbers based on the cpu clock
        SUBROUTINE init_random_seed()
          INTEGER :: i, n, clock
          INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          CALL RANDOM_SEED(size = n)
          ALLOCATE(seed(n))
          CALL SYSTEM_CLOCK(COUNT=clock)
          seed = clock + 37 * (/ (i - 1, i = 1, n) /)
          CALL RANDOM_SEED(PUT = seed)
          DEALLOCATE(seed)
          END SUBROUTINE

! initialize the lattice randomly with spins of 1 or -1 at each site
        SUBROUTINE init_lattice(lat)
          INTEGER:: i, spin(2),lat(latsize)
          REAL:: rand1
          spin(1) = 1
          spin(2) = -1
          DO i= 1, latsize
            CALL RANDOM_NUMBER(rand1)
            lat(i) = spin(int(rand1*(2))+1)
            END DO
          END SUBROUTINE

! equilibrate cluster
        SUBROUTINE equilibrate(lat)
          INTEGER:: i, lat(latsize)
            DO i = 1, eq
              CALL flip_cluster(lat)
              END DO
            END SUBROUTINE

! build a wolff cluster and flip the cluster
        SUBROUTINE flip_cluster(lat)
          INTEGER:: lat(latsize), temp_lat(latsize)
          INTEGER :: queue(latsize)
          REAL:: rand1
          INTEGER:: head, tail
          INTEGER:: site, neighbor
          LOGICAL:: will_bond, parallel

          !Initialize head and tail
          head = 1
          tail = 2

          !Choose random site to start building cluster
          CALL RANDOM_NUMBER(rand1)
          queue(head) = INT(rand1*latsize) + 1

          temp_lat = lat

          !flip first site
          site = queue(head)
          lat(site) = -lat(site)

          !iterate through queue and add to cluster
          DO WHILE (head.NE.tail)

            site = queue(head)

            !Check North neighbor
            neighbor = north(site)
            will_bond = bond(neighbor,lat)
            parallel = ((lat(site).EQ.(lat(neighbor)*(-1))).AND.will_bond)

            IF ((lat(site).EQ.(lat(neighbor)*(-1))).AND.will_bond) THEN
              lat(neighbor) = lat(neighbor) * (-1)
              queue(tail) = neighbor
              tail = tail + 1
              END IF

            !Check east neighbor
            neighbor = east(queue(head))
            will_bond = bond(neighbor,lat)
            IF ((lat(site).EQ.(lat(neighbor)*(-1))).AND.will_bond) THEN
              lat(neighbor) = lat(neighbor) * (-1)
              queue(tail) = neighbor
              tail = tail + 1
              END IF

            !Check south neighbor
            neighbor = south(queue(head))
            will_bond = bond(neighbor,lat)
            IF ((lat(site).EQ.(lat(neighbor)*(-1))).AND.will_bond) THEN
              lat(neighbor) = lat(neighbor) * (-1)
              queue(tail) = neighbor
              tail = tail + 1
              END IF

            !Check west neighbor
            neighbor = west(queue(head))
            will_bond = bond(neighbor,lat)
            IF ((lat(site).EQ.(lat(neighbor)*(-1))).AND.will_bond) THEN
              lat(neighbor) = lat(neighbor) * (-1)
              queue(tail) = neighbor
              tail = tail + 1
              END IF

            !advance head
            head = head + 1
            END DO
          END SUBROUTINE

        !Return the index of the north neighbor
        PURE INTEGER FUNCTION north(site)
          INTEGER, INTENT(IN):: site
          north = MODULO((site- length), latsize ) + (1 - CEILING( &
            MODULO(site - length, latsize)/ FLOAT(latsize))) * latsize
          END FUNCTION north

        !Return the index of the south neighbor
        PURE INTEGER FUNCTION south(site)
          INTEGER, INTENT(IN):: site
          south = MODULO((site + length), latsize ) + (1 - CEILING( &
            MODULO(site + length, latsize)/ FLOAT(latsize))) * latsize
          END FUNCTION south

        !Return the index of the east neighbor
        PURE INTEGER FUNCTION east(site)
          INTEGER, INTENT(IN):: site
          east = site - MODULO(site - 1, length) + MODULO(site, length)
          END FUNCTION east

        !Return the index of the west neighbor
        PURE INTEGER FUNCTION west(site)
          INTEGER, INTENT(IN):: site
          west =  ((site - 1) / length) * length +  MODULO(site - 1, length) &
            + (1 - CEILING(MODULO(site - 1, length)/ FLOAT(length))) * length
          END FUNCTION west

        !Feed in the site a bond is being built to, and return true or false
        !depending on whether a bond is actually built
        LOGICAL FUNCTION bond(neigh, lat)
          INTEGER:: neigh, lat(latsize)
          REAL*8 :: mc, rand1
          INTEGER:: n,s,e,w, ein, ef, echan

          mc = 1 - EXP( - 2 *  (beta))
          !generate random number
          CALL RANDOM_NUMBER(rand1)
          !determine whether a bond is built
          bond = (rand1.LT.mc)
          END FUNCTION bond

        REAL*8 FUNCTION energy(lat)
          INTEGER:: lat(latsize), n, e, s, w, i

          energy = 0

          DO i = 1, latsize
            n = north(i)
            e = east(i)
            s = south(i)
            w = west(i)

            energy = energy + lat(i) * (lat(n) + lat(e) + lat(s) + lat(w))
            END DO

          END FUNCTION energy

      END PROGRAM
