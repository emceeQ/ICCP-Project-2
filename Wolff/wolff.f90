
PROGRAM wolff
      IMPLICIT NONE
      INTEGER, PARAMETER:: length = 5, latsize = length * length
      INTEGER, PARAMETER:: steps = 50
      REAL*8, PARAMETER:: tau = 10
      INTEGER:: lattice(latsize) = 1
      INTEGER:: i, j
      CHARACTER:: visual(latsize)

      !Initialize the seed
      CALL init_random_seed()

      !Initialize the lattice
      CALL init_lattice(lattice)

      DO i = 1, latsize
        if (lattice(i).EQ.1 ) THEN
          visual(i) = "*"
        else
          visual(i) = "O"
        END IF
        END DO

      WRITE(*, FMT = "(5(A,(1x)))")visual

      !Do the wolff algorithm and take measurements for the given number of
      !steps
      DO j = 1, steps

        print *,

        !Flip a cluster
        CALL flip_cluster(lattice)
        DO i = 1, latsize
        if (lattice(i).EQ.1 ) THEN
          visual(i) = "*"
        else
          visual(i) = "O"
        END IF
        END DO

        WRITE(*, FMT = "(5(A,(1x)))")visual
        END DO


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

! build a wolff cluster and flip the cluster
        SUBROUTINE flip_cluster(lat)
          INTEGER:: lat(latsize)
          INTEGER :: queue(latsize)
          REAL:: rand1
          INTEGER:: head, tail
          INTEGER:: site, neighbor
          LOGICAL:: will_bond

          !Initialize head and tail
          head = 1
          tail = 2

          !Choose random site to start building cluster
          CALL RANDOM_NUMBER(rand1)
          queue(head) = int(rand1*latsize) + 1

          !flip first site
          site = queue(head)
          lat(site) = -lat(site)

          !iterate through queue and add to cluster
          DO WHILE (head.NE.tail)

            site = queue(head)

            !Check North neighbor
            neighbor = north(site)
            will_bond = bond(neighbor,lat)
            IF (((lat(site).EQ.(lat(neighbor)*(-1)))).AND.will_bond) THEN
              lat(neighbor) = lat(neighbor) * (-1)
              queue(tail) = neighbor
              tail = tail + 1
              END IF

            !Check east neighbor
            neighbor = east(queue(head))
            will_bond = bond(neighbor,lat)
            IF (((lat(site).EQ.(lat(neighbor)*(-1)))).AND.will_bond) THEN
              lat(neighbor) = lat(neighbor) * (-1)
              queue(tail) = neighbor
              tail = tail + 1
              END IF

            !Check south neighbor
            neighbor = south(queue(head))
            will_bond = bond(neighbor,lat)
            IF (((lat(site).EQ.(lat(neighbor)*(-1)))).AND.will_bond) THEN
              lat(neighbor) = lat(neighbor) * (-1)
              queue(tail) = neighbor
              tail = tail + 1
              END IF

            !Check west neighbor
            neighbor = west(queue(head))
            will_bond = bond(neighbor,lat)
            IF (((lat(site).EQ.(lat(neighbor)*(-1)))).AND.will_bond) THEN
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
            MODULO(site + length, latsize)/ FLOAT(latsize))) * 25
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
            + (1 - CEILING(MODULO(site - 1, length)/ FLOAT(length))) * 5
          END FUNCTION west

        !Feed in the site a bond is being built to, and return true or false
        !depending on whether a bond is actually built
        LOGICAL FUNCTION bond(neigh, lat)
          INTEGER:: neigh, lat(latsize)
          REAL*8 :: mc, rand1
          INTEGER:: n,s,e,w, ein, ef, echan
          !calculate change in energy
          n = lat(north(neigh))
          s = lat(south(neigh))
          w = lat(west(neigh))
          e = lat(east(neigh))

          ein = lat(neigh) * (n + s + e + w)
          ef = -lat(neigh) * (n + s + e + w)
          echan = ef - ein
          mc = EXP(DFLOAT(echan) / tau)

          !generate random number
          CALL RANDOM_NUMBER(rand1)

          !determine whether a bond is built
          bond = (rand1.LT.0.5)

          END FUNCTION bond

      END PROGRAM