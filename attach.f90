!=============================================================================!
        module global
        
        real    :: Minf, Re, lambda, beta, tw, muw, t0, mu0
        real    :: a0ainf, c, m
        real    :: gamma=1.4, gamma1=0.4
        real    :: xtmin=0.0, xtmax=1.0
        real    :: etamin=0.0, etamax=20.0
        integer :: nx=100, ny=2000, nvar=5

        end module global
!=============================================================================!
        module constants

          real, parameter :: zero    = 0.0000000000000000000d+0
          real, parameter :: pt25    = 2.5000000000000000000d-1
          real, parameter :: pt33    = 3.3333333333333333333d-1
          real, parameter :: pt5     = 5.0000000000000000000d-1
          real, parameter :: pt66    = 6.6666666666666666666d-1
          real, parameter :: one     = 1.0000000000000000000d+0
          real, parameter :: onept25 = 1.2500000000000000000d+0   
          real, parameter :: onept33 = 1.3333333333333333333d+0
          real, parameter :: onept5  = 1.5000000000000000000d+0
          real, parameter :: two     = 2.0000000000000000000d+0
          real, parameter :: three   = 3.0000000000000000000d+0
          real, parameter :: four    = 4.0000000000000000000d+0
          real, parameter :: pi      = 3.1415926535897932385d+0
          real, parameter :: eps     = 1.0000000000000000000d-9

        end module constants
!=============================================================================!
        program fsc_solver 
!
!       Compute the compressible similarity solutions for an attachment-line 
!       boundary layer at Pr=1 using the Stewartson transformation.
!
!       i.e. -- Falkner-Skan-Cooke boundary layer
!
!       Re = U_\infty L / nu_\infty
!       Minf = Uinf / ainf
!       lambda = external sweep angle
!
!       S. Scott Collis
!
!       Written: 5-7-96
!
!       Revised: 8-14-98
!=============================================================================!
        use constants
        use global
        implicit none
        
        integer :: i, j, k
        real :: dxt, deta
        
        external calcdx, fsc
        real, external :: BSITG
        
        real, allocatable :: xt(:), x(:), eta(:), y(:), u(:,:), &
                             us(:), tt0(:), tte(:), uve(:), wve(:)
        real :: theta, dxdxt, aeainf, aea0, ue, we, ys, fpp, gp, t0tn0
        real :: b1, b2, a11, a12, a21, a22, det, bnorm, du1, du2
        real :: ve, delta1, delta2, Re_d1, deltaL
        integer :: iter, key(2)

!.... B-spline variables for integration

        integer :: korder
        real, allocatable :: knot(:), bs(:), int(:)
!=============================================================================!
        write(*,"(/,'Solve the compressible attachment-line eqn',/)")
        write(*,"('Enter M_inf, Re_delta1 = ',$)")
        read(*,*) Minf, Re_d1
        write(*,*)
        
        write(*,"('Enter Lambda (deg.) = ',$)")
        read(*,*) lambda
        write(*,*)
        lambda = pi * lambda / 180.0

        write(*,"('Enter Beta = ',$)")
        read(*,*) beta
        write(*,*)
        m = beta / ( two - beta )
        write(*,*) 'm = ',m
        
        write(*,"('Enter T0, Tw, muw = ',$)")
        read(*,*) t0, tw, muw
        write(*,*)

!.... viscosity varies linearly with temperature

        if (tw .ne. zero) then
          mu0 = muw * t0 / tw
        else
          mu0 = muw
          write(*,*) 'WARNING:  mu0 = muw since tw = 0'
        end if

        a0ainf = sqrt( one + pt5 * gamma1 * Minf**2 )

!.... This value of c makes the local sweep angle equal to lambda at xt = 1

        c = Minf * a0ainf * cos(lambda)
        
        allocate( xt(nx), x(nx) )
        
        x(1) = zero

!.... integrate to get x (not really needed for attachment-line case)

        call runge( x(1), 1, zero, xtmax, nx-1, calcdx, xt, x)
        
        do i = 1, nx
          if (m.eq.zero) then
            aeainf = a0ainf * sqrt( ( one + pt5 * gamma1 * Minf**2 * &
              cos(lambda)**2 ) / ( a0ainf**2 + pt5 * gamma1 * &
              c**2 ) )
          else
            aeainf = a0ainf * sqrt( ( one + pt5 * gamma1 * Minf**2 * &
              cos(lambda)**2 ) / ( a0ainf**2 + pt5 * gamma1 * &
              c**2 * xt(i)**(2*m) ) )
          end if
          aea0 = aeainf / a0ainf
          if (m.eq.zero) then
            theta = 180.0 / pi * atan2(Minf * sin(lambda), aea0 * c )
          else
            theta = 180.0 / pi * atan2(Minf * sin(lambda), aea0 * c * xt(i)**m)
          end if
          dxdxt = ( mu0 * tw ) / ( muw * t0 ) * aea0** &
                  ((one-three*gamma)/gamma1)
          if (m.eq.zero) then
            write(10,10) xt(i), x(i), aeainf, aea0, theta, &
              aea0 * c, Minf * sin(lambda), dxdxt
          else
            write(10,10) xt(i), x(i), aeainf, aea0, theta, &
              aea0 * c * xt(i)**m, Minf * sin(lambda), dxdxt
          end if
        end do

!.... compute some quantities at xt = 0

        aeainf = a0ainf * sqrt( ( one + pt5 * gamma1 * Minf**2 * &
                 cos(lambda)**2 ) / ( a0ainf**2 ) )
        aea0 = aeainf / a0ainf
        ue = aea0 * c * zero**m
        we = Minf * sin(lambda)
        theta = 180.0 / pi * atan2(we, ue)

!.... As defined, ys should be divided by Sqrt(Re)

        ys = one

!.... Output edge values

        write(*,*)
        write(*,*) 'Values at xt = 0'
        write(*,"('    x        AeAinf   AeA0          theta       ue         we         Me         ys')")
        write(*,10) x(1), aeainf, aea0, theta, ue, we, sqrt(ue**2+we**2)/aeainf, ys
        write(*,*)

!=============================================================================!
!.... integrate the compressible FSC equations
!=============================================================================!

        allocate( eta(ny), y(ny), u(nvar,ny), us(nvar), tt0(ny), tte(ny), &
                  uve(ny), wve(ny) )
        
!.... Fixed wall-temperature and non-slip boundary conditions

        t0tn0 = (one + pt5 * gamma1 * Minf**2 ) / &
                (one + pt5 * gamma1 * Minf**2 * (cos(lambda))**2 )
                
        write(*,*) 'T_0/T_N0 = ', t0tn0
        write(*,"('Enter f'''', g'' ==> ',$)")
        read(*,*) fpp, gp
        write(*,*)
        write(*,*)


        u(1,1) = fpp                    ! f'' (guess)
        u(2,1) = zero                   ! f'
        u(3,1) = zero                   ! f
        u(4,1) = gp                     ! g' (guess)
        u(5,1) = zero                   ! g
        key(1) = 1
        key(2) = 4

!.... Begin the Newton iteration

        iter = 0
 100    continue
        iter = iter + 1

!.... save the starting values

        us = u(:,1)
        
        call runge(u(1,1),nvar,etamin,etamax,ny-1,fsc,eta,u)

        write (*,20) iter, u(key(1),1), u(key(2),1), u(2,ny), u(5,ny)
 20     format(1p,i5,1x,4(e13.7,1x))

!.... Newton solve

        b1 = u(2,ny)
        b2 = u(5,ny)
        
!.... Perturb u(1,1)

        u(key(1),1) = (one + eps) * us(key(1))
        u(key(2),1) = us(key(2))
        
        call runge(u(1,1),nvar,etamin,etamax,ny-1,fsc,eta,u)

        a11 = (u(2,ny)-b1)/(eps*us(key(1)))
        a21 = (u(5,ny)-b2)/(eps*us(key(1)))

!.... Perturb u(4,1)

        u(key(1),1) = us(key(1))
        u(key(2),1) = (one + eps) * us(key(2))
        
        call runge(u(1,1),nvar,etamin,etamax,ny-1,fsc,eta,u)
        
        a12 = (u(2,ny)-b1)/(eps*us(key(2)))
        a22 = (u(5,ny)-b2)/(eps*us(key(2)))

!.... solve the 2x2 system for the correction

        det = a11*a22-a21*a12
        
        b1 = one - b1
        b2 = one - b2
        bnorm = sqrt(b1**2 + b2**2)
        
        du1 = (b1*a22 - b2*a12) / det
        du2 = (a11*b2 - a21*b1) / det
        
        u(key(1),1) = u(key(1),1) + du1
        u(key(2),1) = u(key(2),1) + du2
        
        if (bnorm .gt. 1.0e-7 .and. iter .le. 10) goto 100

!.... integrate one last time using the latest wall values

        iter = iter + 1
        call runge(u(1,1),nvar,etamin,etamax,ny-1,fsc,eta,u)

        write (*,20) iter, u(key(1),1), u(key(2),1), u(2,ny), u(5,ny)

!.... converged

        write (*,"(/,'Converged to ',1pe10.3)") bnorm

!.... compute the temperature profile at the attachment-line

        ve = sqrt( (Minf * sin(lambda))**2 )
        do j = 1, ny
          tt0(j) = one + (tw/t0 - one) * (one - u(5,j)) - &
                   (one - one/t0tn0) * u(5,j)**2
          tte(j) = tt0(j) / aea0**2
          uve(j) = 0.0
          wve(j) = Minf * sin(lambda) * u(5,j) /ve
!         write (20,10) eta(j), (u(k,j), k = 1, nvar), tte(j)
!         write (30,10) eta(j), &
!                       u(2,j) * cos(lambda) * cos(lambda) + &
!                       u(5,j) * sin(lambda) * sin(lambda),  &
!                      -u(2,j) * cos(lambda) * sin(lambda) + &
!                       u(5,j) * sin(lambda) * cos(lambda)
!         write (40,10) eta(j), one/tte(j), uve(j), zero, wve(j), tte(j)
        end do

!.... compute y at xt = 0, (I'm currently using B-splines to do this)
!.... Note that the Reynolds number has not been accounted for yet.

        korder = 5
        allocate ( knot(ny+korder), bs(ny) )
        call BSNAK( ny, eta, korder, knot )
        call BSINT( ny, eta, tte, korder, knot, bs )
        y(1) = zero
        do j = 2, ny
          y(j) = y(j-1) + BSITG(eta(j-1),eta(j),korder,knot,ny,bs)
        end do

!.... compute the displacement thickness, needs to be divided by sqrt(Re)

        allocate( int(ny) )
        do j = 1, ny
          int(j) = tte(j) - wve(j)
        end do
        call BSINT( ny, eta, int, korder, knot, bs )
        delta1 = BSITG(eta(1),eta(ny),korder,knot,ny,bs)

!.... compute the momentum thickness, needs to be divided by sqrt(Re)

        do j = 1, ny
          int(j) = wve(j) * ( one - wve(j) )
        end do
        call BSINT( ny, eta, int, korder, knot, bs )
        delta2 = BSITG(eta(1),eta(ny),korder,knot,ny,bs)
        write(*,*) 'I = theta/delta_l = ',delta2

!.... given Re_delta1 as input, compute Re = U^*_e l^* / nu^*_e

        Re = ( Re_d1 / delta1 )**2
        
        write(*,60) delta1, delta2, delta1/delta2, Re
 60     format('Delta_1 = ',1pe13.6,'  Delta_2 = ',1pe13.6, &
               '  H = ',1pe13.6,'  Re = ',1pe13.6)

!.... cprofile.dat is in chordwise coordinates
!.... sprofile.dat is in streamline coordinates

        write(*,"('Enter desired delta/L (0 = similarity) ==> ',$)") 
        read(*,*) deltaL
        write(*,*)

!.... Note that dUe/dx = one at the attachement line

        if (deltaL .ne. zero) then
          eta = y * deltaL
          write(*,70) delta1*deltaL, delta2*deltaL, delta1/delta2
70        format('Delta_1 = ',1pe13.6,'  Delta_2 = ',1pe13.6, &
                 '  H = ',1pe13.6)
        endif

!       if (deltaL .ne. zero) then
!         eta = y / delta1 * deltaL
!         write(*,70) delta1/delta1*deltaL, delta2/delta1*deltaL, &
!                     delta1/delta2
!       endif
        
        open(50,file='cprofile.dat')
        open(60,file='sprofile.dat')
        write(50,*) "# ", 1, ny, 5, eta(ny)
        write(60,*) "# ", 1, ny, 5, eta(ny)
        do j = 1, ny
          write (50,10) eta(j), one/tte(j), uve(j), zero, wve(j), tte(j)
          write (60,10) eta(j), one/tte(j), uve(j) * cos(lambda) + &
                        wve(j) * sin(lambda), zero, &
                        -uve(j) * sin(lambda) + wve(j) * cos(lambda), &
                        tte(j)
        end do
        close(50)
        close(60)
        
        stop
 10     format( 8(1pe10.4,1x) )
        end

!=============================================================================!
        subroutine fsc(eta, u, du)
!
!       The compressible Falkner-Skan-Cooke equations
!
!       Revised:  5-6-96
!=============================================================================!
        use global
        use constants
        implicit none

        real :: eta, u(nvar), du(nvar)
        real :: a1, a2 
!=============================================================================!
        a1 = (one + pt5 * gamma1 * Minf**2 ) / &
             (one + pt5 * gamma1 * Minf**2 * (cos(lambda))**2 )

        a2 = ( tw / t0 - one ) * a1
        a1 = a1 - one
        
        du(1) = beta * ( u(2)**2 - one - a1 * (one - u(5)**2) - &
                a2 * (one - u(5)) ) - u(3) * u(1)
 
        du(2) = u(1)

        du(3) = u(2)

        du(4) = -u(3) * u(4)

        du(5) = u(4)

        return
        end

!=============================================================================!
        subroutine calcdx( xt, x, dxdxt )
!
!       Compute x from the Stewartson's + similarity transformation
!
!=============================================================================!
        use global
        use constants
        implicit none
        
        real :: xt, x, dxdxt
        real :: aeainf, aea0
!=============================================================================!
        if (m .eq. zero ) then
          aeainf = a0ainf * sqrt( ( one + pt5 * gamma1 * Minf**2 * &
            cos(lambda)**2 ) / ( a0ainf**2 + pt5 * gamma1 * c**2  ) )
        else
          aeainf = a0ainf * sqrt( ( one + pt5 * gamma1 * Minf**2 * &
            cos(lambda)**2 ) / ( a0ainf**2 + pt5 * gamma1 * &
            c**2 * xt**(2*m) ) )
        end if
        aea0 = aeainf / a0ainf
        dxdxt = ( mu0 * tw ) / ( muw * t0 ) * aea0**((one-three*gamma)/gamma1)
        return
        end
        
