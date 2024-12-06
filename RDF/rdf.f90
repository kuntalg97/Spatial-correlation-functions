      PROGRAM three_body

      ! Finds the g(r) (RDF)

      ! Kuntal Ghosh
      ! February 2022
      
      IMPLICIT NONE

      INTEGER*8 :: i,j,k,ios,nlines,iframe,nframes,natoms,id
      INTEGER*8 :: ibin_r,nbin_r,nc_bin
      REAL*8, PARAMETER :: dr=0.02d0, dang=1.0d0,sigma=3.7d0,& 
                           r_fixed=3.0d0,pi=4.0d0*ATAN(1.0d0)
      REAL*8, ALLOCATABLE :: xyz(:,:,:),force_xyz(:,:,:),g(:)
      REAL*8 :: v(3),v1(3),v2(3),dmin,dmax,d, &
                rij,arb,box_length(3),norm_fact,rho,total,nc

      OPEN (1, FILE = '../traj.lammpstrj', STATUS = 'OLD')   
      OPEN (2, FILE = 'rdf.dat', STATUS = 'UNKNOWN')

!     Reading the file
      
      nlines = 0
      read_loop: DO
         READ (1,*,IOSTAT=ios)
         IF (ios/=0) EXIT read_loop
         nlines = nlines + 1
      END DO read_loop

      REWIND (1)

      READ (1,*)
      READ (1,*)
      READ (1,*)
      READ (1,*) natoms
      READ (1,*)
      DO i = 1,3
         READ (1,*) arb,box_length(i) !Reads the box length
      END DO      

      REWIND (1)

      rho = natoms/(box_length(1)*box_length(2)*box_length(3))
      nframes = nlines/(natoms + 9) ! Typically, skip the first 9 lines of each frame

      ALLOCATE (xyz(3,natoms,nframes),force_xyz(3,natoms,nframes))

!     Storing all the coordinates

      DO i = 1,nframes
         DO j = 1,9
            READ (1,*)
         END DO
         DO k = 1,natoms
            READ (1,*) id,xyz(1:3,k,i),force_xyz(1:3,k,i)
         END DO
      END DO

      CLOSE (1)

      dmin = 0.0d0
      dmax = 0.5d0*box_length(1)

      WRITE (*,*) "dmin, dmax", dmin,dmax

!     Computing g(r)

      nbin_r = INT(0.5d0*box_length(1)/dr) + 1      

      ALLOCATE (g(nbin_r))

      g(:) = 0.0d0

      DO iframe = 1,nframes
         DO i = 1,natoms-1
            DO j = i+1,natoms
       
               v(1:3) = xyz(1:3,j,iframe) - xyz(1:3,i,iframe)
               v(1:3) = v(1:3) - box_length(1:3)*NINT(v(1:3)/box_length(1:3))
               d = DSQRT(DOT_PRODUCT(v,v))                           

               IF (d<(0.5d0*box_length(1)) .AND. d>1.0d0) THEN
                  ibin_r = INT((d-dmin)/dr) + 1
               END IF                           

               IF (ibin_r>0 .AND. ibin_r<=nbin_r) THEN
                  g(ibin_r) = g(ibin_r) + 1.0d0
               END IF

            END DO
         END DO
      END DO
         
      d = dmin
      DO ibin_r = 1,nbin_r
         g(ibin_r) = g(ibin_r)/(4.0d0*DFLOAT(natoms)*dr*pi*(d**2)*DFLOAT(nframes)*rho)
         WRITE (2,*) d,g(ibin_r)
         d = d + dr
      END DO

      CLOSE (2)

!     Computing Nc (coordination number)
!     Trapezoidal rule for integration used here (can be used later)      

      total = 0.0d0
      d = dmin
      DO ibin_r = 3,nbin_r-1 !First bin is useless
         IF (d>sigma) EXIT
         total = total + 2.0d0*g(ibin_r)*(d**2)
         d = d + dr
      END DO
      nc_bin = INT((sigma-dmin)/dr) + 1
      nc = 0.5d0*dr*4.0d0*pi*rho*(total + g(2)*(dr**2) + g(nc_bin)*(sigma**2))
      WRITE (*,*) "Coordination number =", nc

      DEALLOCATE (xyz,force_xyz,g)

      END PROGRAM three_body

