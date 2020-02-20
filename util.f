
C***********************************************************************
C>    \brief Advance one time step using fourth order (real) Runge-Kutta
C>    \param[in] neq number of equations
C>    \param[in] yo initial value
C>    \param[out] yf final value
C>    \param[in] to intial time
C>    \param[in] h time step
C>    \param[in] FUNC function to integrate
C***********************************************************************
      subroutine SRK4(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step using fourth order (real) Runge-Kutta
C
C***********************************************************************
      external    FUNC
      integer     neq
      real        to, h
      real        yo(neq), yf(neq)
      real        f(neq), k1(neq), k2(neq), k3(neq), k4(neq), q(neq)
      
      call FUNC(neq, yo, to, f)
      do j = 1 , neq
        k1(j) = h*f(j)
        q(j) = yo(j) + 0.5*k1(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k2(j) = h*f(j)
        q(j) = yo(j) + 0.5*k2(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k3(j) = h*f(j)
        q(j) = yo(j) + k3(j)
      end do
      call FUNC(neq, q, to+h, f)
      do j = 1 , neq
        k4(j) = h*f(j)
        yf(j) = yo(j)+k1(j)/6.+(k2(j)+k3(j))/3.+k4(j)/6.
      end do

      return
      end

C***********************************************************************
C>    \brief Advance one time step using fourth order (real) Runge-Kutta
C>    \param[in] neq number of equations
C>    \param[in] yo initial value
C>    \param[out] yf final value
C>    \param[in] to intial time
C>    \param[in] h time step
C>    \param[in] FUNC function to integrate
C***********************************************************************
      subroutine CRK4(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step using fourth order (complex) Runge-Kutta
C
c***********************************************************************
      external    FUNC
      integer     neq
      real        to, h
      complex     yo(neq), yf(neq)
      complex     f(neq), k1(neq), k2(neq), k3(neq), k4(neq), q(neq)
      
      call FUNC(neq, yo, to, f)
      do j = 1 , neq
        k1(j) = h*f(j)
        q(j) = yo(j) + 0.5*k1(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k2(j) = h*f(j)
        q(j) = yo(j) + 0.5*k2(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k3(j) = h*f(j)
        q(j) = yo(j) + k3(j)
      end do
      call FUNC(neq, q, to+h, f)
      do j = 1 , neq
        k4(j) = h*f(j)
        yf(j) = yo(j)+k1(j)/6.+(k2(j)+k3(j))/3.+k4(j)/6.
      end do

      return
      end

C***********************************************************************
      SUBROUTINE SPLINE (N,X,Y,FDP)      
C***********************************************************************
C
C.... Note: this routine is in the public domain and available
C     at https://web.stanford.edu/class/me200c/
C
C-----THIS SUBROUTINE COMPUTES THE SECOND DERIVATIVES NEEDED 
C-----IN CUBIC SPLINE INTERPOLATION.  THE INPUT DATA ARE:    
C-----N = NUMBER OF DATA POINTS          
C-----X = ARRAY CONTAINING THE VALUES OF THE INDEPENDENT VARIABLE      
C-----    (ASSUMED TO BE IN ASCENDING ORDER)       
C-----Y = ARRAY CONTAINING THE VALUES OF THE FUNCTION AT THE 
C-----    DATA POINTS GIVEN IN THE X ARRAY         
C-----THE OUTPUT IS THE ARRAY FDP WHICH CONTAINS THE SECOND  
C-----DERIVATIVES OF THE INTERPOLATING CUBIC SPLINE.         
      DIMENSION X(N),Y(N),A(N),B(N),C(N),R(N),FDP(N)  
C-----COMPUTE THE COEFFICIENTS AND THE RHS OF THE EQUATIONS. 
C-----THIS ROUTINE USES THE CANTILEVER CONDITION.  THE PARAMETER       
C-----ALAMDA (LAMBDA) IS SET TO 1. BUT THIS CAN BE USER-MODIFIED.      
C-----A,B,C ARE THE THREE DIAGONALS OF THE TRIDIAGONAL SYSTEM;         
C-----R IS THE RIGHT HAND SIDE.  THESE ARE NOW ASSEMBLED.    
      ALAMDA = 1.    
      NM2 = N - 2    
      NM1 = N - 1    
      C(1) = X(2) - X(1)       
      DO 1 I=2,NM1   
      C(I) = X(I+1) - X(I)     
      A(I) = C(I-1)  
      B(I) = 2.*(A(I) + C(I))  
      R(I) = 6.*((Y(I+1) - Y(I))/C(I) - (Y(I) - Y(I-1))/C(I-1))        
    1 CONTINUE       
      B(2) = B(2) + ALAMDA * C(1)        
      B(NM1) = B(NM1) + ALAMDA * C(NM1)  
C-----AT THIS POINT WE COULD CALL A TRIDIAGONAL SOLVER SUBROUTINE      
C-----BUT THE NOTATION IS CLUMSY SO WE WILL SOLVE DIRECTLY.  THE       
C-----NEXT SECTION SOLVES THE SYSTEM WE HAVE JUST SET UP.    
      DO 2 I=3,NM1   
      T = A(I)/B(I-1)          
      B(I) = B(I) - T * C(I-1) 
      R(I) = R(I) - T * R(I-1) 
    2 CONTINUE       
      FDP(NM1) = R(NM1)/B(NM1) 
      DO 3 I=2,NM2   
      NMI = N - I    
      FDP(NMI) = (R(NMI) - C(NMI)*FDP(NMI+1))/B(NMI)         
    3 CONTINUE       
      FDP(1) = ALAMDA * FDP(2) 
      FDP(N) = ALAMDA * FDP(NM1)         
C-----WE NOW HAVE THE DESIRED DERIVATIVES SO WE RETURN TO THE          
C-----MAIN PROGRAM.  
      RETURN         
      END  

C***********************************************************************
      SUBROUTINE SPEVAL (N,X,Y,FDP,XX,F) 
C***********************************************************************
C
C.... Note: this routine is in the public domain and available
C     at https://web.stanford.edu/class/me200c/
C
C-----THIS SUBROUTINE EVALUATES THE CUBIC SPLINE GIVEN       
C-----THE 2ND DERIVATIVE COMPUTED BY SUBROUTINE SPLINE.          
C-----THE INPUT PARAMETERS N,X,Y,FDP HAVE THE SAME 
C-----MEANING AS IN SPLINE.    
C-----XX = VALUE OF INDEPENDENT VARIABLE FOR WHICH 
C-----     AN INTERPOLATED VALUE IS REQUESTED      
C-----F =  THE INTERPOLATED RESULT       
      DIMENSION X(N),Y(N),FDP(N)      
C-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
#if USE_NR_HUNT
c
c     Search using bisection with a good guess
c
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call HUNT (X,N,XX,I)
      END IF
      IOLD = I
#elif 1
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call BISECT (X,N,XX,I)
      ENDiF
      IOLD = I
#else
c
c     This is really a slow way of searching
c
      NM1 = N - 1
      DO 1 I=1,NM1
      IF (XX.LE.X(I+1)) GO TO 10
    1 CONTINUE   
#endif
C-----NOW EVALUATE THE CUBIC   
   10 DXM = XX - X(I)          
      DXP = X(I+1) - XX        
      DEL = X(I+1) - X(I)      
      F = FDP(I)*DXP*(DXP*DXP/DEL - DEL)/6.        
     1   +FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6.     
     2   +Y(I)*DXP/DEL + Y(I+1)*DXM/DEL 
      RETURN        
      END 

C***********************************************************************
      subroutine BISECT(X,N,XX,I)
C***********************************************************************
      dimension X(N)
C***********************************************************************
      il = I-1
      ir = N-1
      do while (ir-il .gt. 1) 
        im = ISHFT(ir+il,-1) 
        if ( X(im+1) > xx ) then
          ir = im
        else
          il = im
        end if
      end do
      I = il+1
      return
      end

C***********************************************************************
      SUBROUTINE SPDER(N,X,Y,FDP,XX,F,FP,FPP)
C***********************************************************************
C
C.... Note: this routine is in the public domain and available
C     at https://web.stanford.edu/class/me200c/
C
C-----THIS SUBROUTINE EVALUATES THE CUBIC SPLINE GIVEN       
C-----THE 2ND DERIVATIVE COMPUTED BY SUBROUTINE SPLINE.          
C-----THE INPUT PARAMETERS N,X,Y,FDP HAVE THE SAME 
C-----MEANING AS IN SPLINE.    
C-----XX = VALUE OF INDEPENDENT VARIABLE FOR WHICH 
C-----     AN INTERPOLATED VALUE IS REQUESTED      
C-----F =  THE INTERPOLATED RESULT       
C-----FP = THE INTERPOLATED DERIVATIVE RESULT       
      INTEGER N
      DIMENSION X(N),Y(N),FDP(N)
      REAL XX, F, FP, FPP
C-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
#if USE_NR_HUNT
c
c     Search using bisection with a good guess
c
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call HUNT (X,N,XX,I)
      END IF
      IOLD = I
#elif 1
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call BISECT (X,N,XX,I)
      ENDiF
      IOLD = I
#else
c
c     This is really a slow way of searching
c
      NM1 = N - 1
      DO 1 I=1,NM1
      IF (XX.LE.X(I+1)) GO TO 10
    1 CONTINUE   
#endif
C-----NOW EVALUATE THE CUBIC   
   10 continue
C     write(*,*) I, X(I), XX, X(I+1)
      DXM = XX - X(I)
      DXP = X(I+1) - XX
      DEL = X(I+1) - X(I)
      F = FDP(I)*DXP*(DXP*DXP/DEL - DEL)/6.0
     1   +FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6.0
     2   +Y(I)*DXP/DEL + Y(I+1)*DXM/DEL
      FP= FDP(I)*(-3.0*DXP*DXP/DEL + DEL)/6.0
     1   +FDP(I+1)*(3.0*DXM*DXM/DEL - DEL)/6.0
     2   -Y(I)/DEL + Y(I+1)/DEL
      FPP=FDP(I)*DXP/DEL+FDP(I+1)*DXM/DEL
      RETURN     
      END

C***********************************************************************
      subroutine SLSRK14(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step with low-storage Runge Kutta 14
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h
      real     yo(neq), yf(neq)
      real     f(neq), yt(neq)
C***********************************************************************
      parameter (nsubstep=14)
      real A(0:14), B(0:13), c(0:13), w(0:14)
      data A / 0.0, 
     &         -0.7188012108672410, 
     &         -0.7785331173421570,
     &         -0.0053282796654044, 
     &         -0.8552979934029281, 
     &         -3.9564138245774565, 
     &         -1.5780575380587385,
     &         -2.0837094552574054, 
     &         -0.7483334182761610,
     &         -0.7032861106563359,  
     &          0.0013917096117681,
     &         -0.0932075369637460, 
     &         -0.9514200470875948,
     &         -7.1151571693922548, 
     &         0.0/
      data B / 0.0367762454319673,
     &         0.3136296607553959,
     &         0.1531848691869027,
     &         0.0030097086818182,
     &         0.3326293790646110,
     &         0.2440251405350864,
     &         0.3718879239592277,
     &         0.6204126221582444,
     &         0.1524043173028741,
     &         0.0760894927419266,
     &         0.0077604214040978,
     &         0.0024647284755382,
     &         0.0780348340049386,
     &         5.5059777270269628 /
      data c / 0.0,
     &         0.0367762454319673,
     &         0.1249685262725025,
     &         0.2446177702277698,
     &         0.2476149531070420,
     &         0.2969311120382472,
     &         0.3978149645802642,
     &         0.5270854589440328,
     &         0.6981269994175695,
     &         0.8190890835352128,
     &         0.8527059887098624,
     &         0.8604711817462826,
     &         0.8627060376969976,
     &         0.8734213127600976 /
      data w / -0.116683473041717417,
     &         0.213493962104674251,
     &         0.128620987881127052,
     &         4.610096100109887907,
     &         -5.386527768056724064,
     &         1.445540684241274576,
     &         -0.761388932107154526,
     &         0.543874700576422732,
     &         0.102277834602298279,
     &         0.07127466608688701188,
     &         -3.459648919807762457,
     &         37.20095449534884580,
     &         -39.09786206496502814,
     &         5.505977727026962754,
     &         0.0 /
      do j = 1, neq
        yf(j) = yo(j)
      end do
      do i = 0, nsubstep-1
        t = to + c(i)*h
        call FUNC(neq, yf, t, f)
        do j = 1, neq
          yt(j) = A(i)*yt(j) + h*f(j)
        end do
        do j = 1, neq
          yf(j) = yf(j) + B(i)*yt(j)
        end do
        if (i+1 .lt. nsubstep) then
          t = to + c(i+1)*h
        else
          t = to + h
        end if
      end do
      return
      end

C***********************************************************************
      subroutine CLSRK14(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step with low-storage Runge Kutta 14
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h
      complex  yo(neq), yf(neq)
      complex  f(neq), yt(neq)
C***********************************************************************
      parameter (nsubstep=14)
      real A(0:14), B(0:13), c(0:13), w(0:14)
      data A / 0.0, 
     &         -0.7188012108672410, 
     &         -0.7785331173421570,
     &         -0.0053282796654044, 
     &         -0.8552979934029281, 
     &         -3.9564138245774565, 
     &         -1.5780575380587385,
     &         -2.0837094552574054, 
     &         -0.7483334182761610,
     &         -0.7032861106563359,  
     &          0.0013917096117681,
     &         -0.0932075369637460, 
     &         -0.9514200470875948,
     &         -7.1151571693922548, 
     &         0.0/
      data B / 0.0367762454319673,
     &         0.3136296607553959,
     &         0.1531848691869027,
     &         0.0030097086818182,
     &         0.3326293790646110,
     &         0.2440251405350864,
     &         0.3718879239592277,
     &         0.6204126221582444,
     &         0.1524043173028741,
     &         0.0760894927419266,
     &         0.0077604214040978,
     &         0.0024647284755382,
     &         0.0780348340049386,
     &         5.5059777270269628 /
      data c / 0.0,
     &         0.0367762454319673,
     &         0.1249685262725025,
     &         0.2446177702277698,
     &         0.2476149531070420,
     &         0.2969311120382472,
     &         0.3978149645802642,
     &         0.5270854589440328,
     &         0.6981269994175695,
     &         0.8190890835352128,
     &         0.8527059887098624,
     &         0.8604711817462826,
     &         0.8627060376969976,
     &         0.8734213127600976 /
      data w / -0.116683473041717417,
     &         0.213493962104674251,
     &         0.128620987881127052,
     &         4.610096100109887907,
     &         -5.386527768056724064,
     &         1.445540684241274576,
     &         -0.761388932107154526,
     &         0.543874700576422732,
     &         0.102277834602298279,
     &         0.07127466608688701188,
     &         -3.459648919807762457,
     &         37.20095449534884580,
     &         -39.09786206496502814,
     &         5.505977727026962754,
     &         0.0 /
      do j = 1, neq
        yf(j) = yo(j)
      end do
      do i = 0, nsubstep-1
        t = to + c(i)*h
        call FUNC(neq, yf, t, f)
        do j = 1, neq
          yt(j) = A(i)*yt(j) + h*f(j)
        end do
        do j = 1, neq
          yf(j) = yf(j) + B(i)*yt(j)
        end do
        if (i+1 .lt. nsubstep) then
          t = to + c(i+1)*h
        else
          t = to + h
        end if
      end do
      return
      end

C***********************************************************************
      subroutine advance(FUNC, neq, t1, t2, nstep, t, x)
C***********************************************************************
C
C       Advance from t1 to t2
C
C***********************************************************************
        external FUNC
        real t(nstep+1), x(neq,nstep+1)
        real dt
C***********************************************************************
        dt = (t2 - t1)/nstep
        t(1) = t1
        do i = 1, nstep
          call srkck45(neq, x(1,i), x(1,i+1), t(i), dt, FUNC)
          t(i+1) = t(i) + dt
        end do
        return
      end

C***********************************************************************
      subroutine SRKCK45(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step Runge-Kutta Cash-Karp method
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h, t
      real     yo(neq), yf(neq), yt(neq)
      real     yk(neq,6), ye(neq)
C***********************************************************************
      real b(6,5)
      real a(6), c(6), d(6)
      data a / 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 /
      data b / 0.0, 0.2, 0.075, 0.3, -0.2037037037037037,
     &         0.029495804398148147,
     &         0.0, 0.0, 0.225, -0.9, 2.5, 0.341796875,
     &         0.0, 0.0, 0.0, 1.2, -2.5925925925925926,
     &         0.041594328703703706,
     &         0.0, 0.0, 0.0, 0.0, 1.2962962962962963,
     &         0.40034541377314814,
     &         0.0, 0.0, 0.0, 0.0, 0.0, 0.061767578125 /
      data c / 0.09788359788359788, 0.0, 0.4025764895330113,
     &         0.21043771043771045, 0.0, 0.2891022021456804 /
      data d / -0.004293774801587311, 0.0, 0.018668586093857853,
     &         -0.034155026830808066, -0.019321986607142856,
     &         0.03910220214568039 /
c
c     Test data
c
#ifdef OS_DEBUG
      do i = 1, 6
        do j = 1, 5
          write(*,*) i, j, b(i,j)
        end do
      end do
      stop
#endif
c
c     Stage 1 - 6
c
      do m = 1, 6
        t = to + a(m)*h
        do n = 1, neq
          yt(n) = yo(n)
        end do
        do k = 1, m-1
          do n = 1, neq
            yt(n) = yt(n) + b(m,k)*yk(n,k)
          end do
        end do
        call FUNC(neq, yt, t, yk(1,m))
        do n = 1, neq
          yk(n,m) = h * yk(n,m)
        end do
      end do
c
c     Final solution and error
c
      do n = 1, neq
        yf(n) = yo(n)
        ye(n) = 0.0
      end do
      do k = 1, 6
        do n = 1, neq
          yf(n) = yf(n) + c(k)*yk(n,k)
          ye(n) = ye(n) + d(k)*yk(n,k)
        end do
      end do

      return
      end

C***********************************************************************
      subroutine CRKCK45(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step Runge-Kutta Cash-Karp method
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h, t
      complex  yo(neq), yf(neq), yt(neq)
      complex  yk(neq,6), ye(neq)
C***********************************************************************
      real b(6,5)
      real a(6), c(6), d(6)
      data a / 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 /
      data b / 0.0, 0.2, 0.075, 0.3, -0.2037037037037037, 
     &         0.029495804398148147,
     &         0.0, 0.0, 0.225, -0.9, 2.5, 0.341796875,
     &         0.0, 0.0, 0.0, 1.2, -2.5925925925925926, 
     &         0.041594328703703706,
     &         0.0, 0.0, 0.0, 0.0, 1.2962962962962963,
     &         0.40034541377314814,
     &         0.0, 0.0, 0.0, 0.0, 0.0, 0.061767578125 /
      data c / 0.09788359788359788, 0.0, 0.4025764895330113,
     &         0.21043771043771045, 0.0, 0.2891022021456804 /
      data d / -0.004293774801587311, 0.0, 0.018668586093857853,
     &         -0.034155026830808066, -0.019321986607142856,
     &         0.03910220214568039 / 
c
c     Test data
c
#ifdef OS_DEBUG
      do i = 1, 6
        do j = 1, 5
          write(*,*) i, j, b(i,j)
        end do
      end do
      stop
#endif
c
c     Stage 1 - 6
c
      do m = 1, 6 
        t = to + a(m)*h
        do n = 1, neq
          yt(n) = yo(n)
        end do
        do k = 1, m-1
          do n = 1, neq
            yt(n) = yt(n) + b(m,k)*yk(n,k)
          end do
        end do 
        call FUNC(neq, yt, t, yk(1,m))
        do n = 1, neq
          yk(n,m) = h * yk(n,m)
        end do
      end do
c
c     Final solution and error
c
      do n = 1, neq
        yf(n) = yo(n)
        ye(n) = 0.0
      end do
      do k = 1, 6
        do n = 1, neq
          yf(n) = yf(n) + c(k)*yk(n,k)
          ye(n) = ye(n) + d(k)*yk(n,k)
        end do
      end do

      return
      end 
