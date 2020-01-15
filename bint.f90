program plot
      
      implicit none
      
      integer :: i,j,k,left,leftmk,mflag,n,npoint
      real :: dx, t(10), values(7), x, xl

      data n, k /7,3/, t /3*0.,2*1.,3.,4.,3*6./
      data values /7*0./, npoint /31/

      xl = t(k)
      dx = (t(n+1)-t(k))/float(npoint-1)

      do i = 1, npoint
	x = xl + float(i-1)*dx
	call interv( t, n, x, left, mflag )
	leftmk = left - k
	call bsplvb( t, k, 1, x, left, values(leftmk+1) )
	write (*,*) x, (values(j), j=3,7)
10	format(10(1pe12.6,1x))
	do j = 1, k 
	  values(leftmk+j) = 0.0
	end do
      end do

      stop

end program plot

subroutine interv( xt, lxt, x, left, mflag )

      integer left, lxt, mflag, ihi, ilo, istep, middle
      real x, xt(lxf)
      data ilo /1/
      ihi = ilo + 1

      if (ihi .ge. lxt) then

	if (x .ge. xt(lxt)) then
	  mflag = 1
	  left = lxt
	  return
	end if

	if (lxt .le. 1) then
	  mflag = -1
	  left = 1
	  return
	end if

	ilo = lxt - 1
	ihi = lxt

      else if (x .lt. xt(ihi)) then
	
	if (x .ge. xt(ilo)) then
	  mflag = 0
	  left = ilo
	  return
	end if
	  
	
subroutine bsint( n, tau, gtau, k, t, bcoef )
!
!.... This is the original calling sequence
!
! subroutine bint( tau, gtau, t, n, k, q, bcoef, iflag )

      implicit none

      integer :: iflag, k, n, i, ilp1mx, j, jj, km1, kpkm2, left, lenq, np1
      real :: bcoef(n), gtau(n), q(2*k-1*n), t(n+k), tau(n), taui

      np1 = n + 1
      km1 = k - 1
      kpkm2 = 2 * km1
      left = k

      lenq = n * ( k + km1 )
      do i = 1, lenq
	q(i) = 0.
      end do

      do i = 1, n
	taui = tau(i)
	ilp1mx = min(i+k,np1)
	left = max(left,I)
	if (taui .lt. t(left)) goto 998
15	if (taui .lt. t(left+1)) goto 16
	left = left + 1
	if (left .lt. ilp1mx) goto 15
	left = left - 1
	if (taui .gt. t(left+1)) goto 998

16	call bsplvb (t, k, 1, taui, left, bcoef )

	jj = i - left + 1 + (left-k)*(k+km1)
	do j = 1, k
	  jj = jj + kpkm2
	  q(jj) = bcoef(j)
	end do
      end do

      call banfac( q, k+km1, n, km1, km1, iflag )

      if (iflag .eq. 1) goto 999
      
      do i = 1, n
	bcoef(i) = gtau(i)
      end do

      call banslv( q, k+km1, n, km1, km1, bcoef )

      return

998   iflag = 2
999   write(*,699)
699   format('Linear system in bint is not invertible')
      
      return
      end

subroutine bsplvb( t, jhigh, index, x, left, biatx )

      implicit none
      
      integer :: index, jhigh, left, i, j=1, jp1
      real    :: biatx(jhigh), t(1), x, deltal(jhigh), deltar(jhigh), &
                 saved, term

      if (index .eq. 1) then
	j = 1
	biatx(1) = 1.
	do while (j .lt. jhigh)
	  jp1 = j + 1
	  deltar(j) = t(left+j) - x
	  deltal(j) = x - t(left+1-j)
	  saved = 0.
	  do i = 1, j
	    term = biatx(i) / (deltar(i) + deltal(jp1-i))
	    biatx(i) = saved + deltar(i)*term
	    saved = deltal(jp1-i)*term
	  end do
	  biatx(jp1) = saved
	  j = jp1
	end do
      end if

      return

end subroutine bsplvb
	
subroutine banfac( w, nroww, nrow, nbandl, nbandu, iflag )

      implicit none

      integer :: iflag, nbandl, nbandu, nrow, nroww
      integer :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
      real :: w(nroww,nrow), factor, pivot

      iflag = 1
      middle = nbandu - 1

      nrowm1 = nrow - 1
      
      if (nrowm1 .eq. 0) then
	iflag = 2
	return
      else if (nrowm1 .eq. 1) then
	if (w(middle,nrow) .ne. 0) return
	iflag = 2
	return
      end if

      if (nbandl .le. 0) then
	do i = 1, nrowm1
	  if (w(middle,i) .eq. 0) then
	    iflag = 2
	    return
	  end if
	end do
      end if

      if (nbandu .le. 0) then
	do i = 1, nrowm1
	  pivot = w(middle,i)
	  if (pivot .eq. 0.) then
	    iflag = 2
	    return
	  end if
	  jmax = min(nbandl, nrow - i)
	  do j = 1, jmax
	    w(middle+j,i) = w(middle+j,i) / pivot
	  end do
	end do
	if (w(middle,nrow) .ne. 0) then
	  return
	else
	  iflag = 2
	  return
	end if
      end if

      do i = 1, nrowm1
	pivot = w(middle,i)
	if (pivot .eq. 0) then
	  iflag = 2
	  return
	end if
	jmax = min(nbandl,nrow-i)
	do j = 1, jmax
	  w(middle+j,i) = w(middle+j,i) / pivot
	end do
	kmax = min(nbandu,nrow-i)
	do k = 1, kmax
	  ipk = i + k
	  midmk = middle - k
	  factor = w(midmk,ipk)
	  do j = 1, jmax
	    w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i) * factor
	  end do
	end do
      end do

      if (w(middle,nrow) .eq. 0) iflag = 2

      return
end subroutine banfac

subroutine banslv( w, nroww, nrow, nbandl, nbandu, b )

      implicit none
      
      integer :: nbandl, nbandu, nrow, nroww
      integer :: i,j,jmax,middle,nrowm1
      real :: w(nroww,nrow), b(nrow)

      middle = nbandu + 1
      if (nrow .eq. 1) goto 49
      nrowm1 = nrow - 1
      if (nbandl .eq. 0) goto 30

      do i = 1, nrowm1
	jmax = min(nbandl, nrow-i)
	do j = 1, jmax
	  b(i+j) = b(i+j) - b(i)*w(middle+j,i)
	end do
      end do

30    if (nbandu .gt. 0) goto 40

      do i = 1, nrow
	b(i) = b(i) / w(1,i)
      end do
      return

40    do i = nrow,2,-1
	b(i) = b(i)/w(middle,i)
	jmax = min(nbandu,i-1)
	do j = 1, jmax
	  b(i-j) = b(i-j) - b(i) * w(middle-j,i)
	end do
      end do

49    b(1) = b(1) / w(middle,1)

      return
      end
