      PROGRAM three_body

      ! Finds the ADF

      ! Kuntal Ghosh
      ! February 2022
      
      IMPLICIT NONE

      INTEGER*8 :: i,j,k,ios,nlines,iframe,nframes,natoms,id
      INTEGER*8 :: ibin_r,ibin_ang,nbin_r,nbin_ang,total
      REAL*8, PARAMETER :: dr=0.02d0, dang=1.0d0,sigma=10.0d0,& 
                           r_fixed=3.0d0,pi=4.0d0*ATAN(1.0d0)
      REAL*8, ALLOCATABLE :: xyz(:,:,:),force_xyz(:,:,:),p(:,:),a(:)
      REAL*8 :: v(3),v1(3),v2(3),dmin,dmax,d,angmin,angmax, &
                ang,rij,arb,box_length(3),norm_fact,rho,d1,d2

      OPEN (1, FILE = '../../cg.lammpstrj', STATUS = 'OLD')
      OPEN (2, FILE = 'adf.dat', STATUS = 'UNKNOWN')

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
      WRITE (*,*) box_length(1:3)

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


      angmin = 0.0d0
      angmax = 180.0d0
      dmin = 0.0d0
      dmax = 0.5d0*box_length(1)

      WRITE (*,*) "angmin, angmax, dmin, dmax", angmin,angmax,dmin,dmax

      nbin_ang = INT((angmax-angmin)/dang) + 1

      ALLOCATE (p(nbin_r,nbin_ang),a(nbin_ang))

!     Computing a(r) !Angular distribution function

      a(:) = 0.0d0

      DO iframe = 1,nframes
         DO i = 1,natoms
            DO j = 1,natoms
               IF (i/=j) THEN

                  v1(1:3) = xyz(1:3,j,iframe) - xyz(1:3,i,iframe)                  
                  v1(1:3) = v1(1:3) - box_length(1:3)*NINT(v1(1:3)/box_length(1:3))
                  d1 = DOT_PRODUCT(v1,v1)

                  DO k = j+1,natoms
                     IF (k/=i) THEN

                        v2(1:3) = xyz(1:3,k,iframe) - xyz(1:3,i,iframe)                        
                        v2(1:3) = v2(1:3) - box_length(1:3)*NINT(v2(1:3)/box_length(1:3))
                        d2 = DOT_PRODUCT(v2,v2)
                        
                        IF (d1<sigma .AND. d2<sigma) THEN
                           CALL angle (v1,v2,ang)

                           IF (ang>0.0d0 .AND. ang<180.0d0) THEN
                              ibin_ang = INT((ang-angmin)/dang) + 1
                           END IF

                           IF (ibin_ang>0 .AND. ibin_ang<=nbin_ang) THEN
                              a(ibin_ang) = a(ibin_ang) + 1.0d0
                              total = total + 1
                           END IF

                        END IF

                     END IF
                  END DO 
                  
               END IF
            END DO
         END DO
      END DO
         
      ang = angmin
      DO ibin_ang = 1,nbin_ang
         WRITE (2,*) (ang*pi/180.0d0),a(ibin_ang)/(DFLOAT(total)*dang)
         ang = ang + dang
      END DO

      CLOSE (3)


      DEALLOCATE (xyz,force_xyz,p,a)

      END PROGRAM three_body

      


      SUBROUTINE angle (v1,v2,a)

      IMPLICIT NONE

      REAL*8 :: cos_theta,v1(3),v2(3),d1,d2,a
      REAL*8, PARAMETER :: pi = 4.0d0*ATAN(1.0d0)

      cos_theta = DOT_PRODUCT(v1,v2)
      d1 = DSQRT(DOT_PRODUCT(v1,v1))
      d2 = DSQRT(DOT_PRODUCT(v2,v2))
      cos_theta = cos_theta/(d1*d2)
      a = DACOS(cos_theta)
      a = (a*180.0d0)/pi

      END SUBROUTINE angle   
