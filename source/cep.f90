! This program solves the non-LTE problem for an atom with many levels in a plane-parallel
! isothermal atmosphere using an iterative algorithm based on the coupled escape probability developed 
! by Elitzur & Asensio Ramos
! Andres Asensio Ramos   (07/04/2005) -> inclusion of the coupled escape probability method (Elitzur & Asensio Ramos)
! Andres Asensio Ramos   (30/11/2005) -> CEP method with Newton with analytical derivatives
! Andres Asensio Ramos   (29/3/2006) -> Variable physical conditions from file

module constants
implicit none

! ---------------------------------------------------------
! PHYSICAL CONSTANTS AND RELATIONS AMONG THEM
! ---------------------------------------------------------
! PI : pi
! PE : electron charge
! PME : electron mass
! PK : Bolztmann constant
! PH : Planck constant
! PC : light speed
! CI :
! PI8C2 : 8 * pi / c^2
! PI4H : 4 * pi / h
! PHK : h / k
! ---------------------------------------------------------

	real(kind=8), parameter :: PI = 3.141592654d0, PE = 4.8032d-10, PME = 9.10938188d-28, UMA = 1.66053873d-24
	real(kind=8), parameter :: PK = 1.3806503d-16, PH = 6.62606876d-27, PC = 2.99792458d10, CI = 2.0706d-16
	real(kind=8), parameter :: PI8C2 = 8.d0 * PI / (PC**2), PI4H = 4.d0 * PI / PH, PHK = PH / PK
	real(kind=8), parameter :: TWOHC2 = 2.d0 * PH / PC**2

! ---------------------------------------------------------
! ADDITIONAL CONSTANTS
! ---------------------------------------------------------
! nordm : maximum order of the Ng acceleration
! n_ptos : maximum number of points
! n_iters : maximum number of iterations
! ---------------------------------------------------------

	integer, parameter :: nordm = 4
   integer, parameter :: n_ptos = 8001, n_iters = 4000	

end module constants

module variables
use constants

! ---------------------------------------------------------
! VARIABLES
! ---------------------------------------------------------
! nz : number of points in the atmosphere
! nfrq : number of frequency points in the profile
! n_freq_ud_doppler : number of points per Doppler width in the profile
! n_freq_perfil : number of Doppler widths for the line
! nang : number of angles for the angular integration
! taufac : aumento en tau en cada punto con respecto al anterior tau(i) = taufac * tau(i-1) (NOT USED)
! col_density : maximum optical depth
! sclht : escala de alturas (NOT USED)
! tempc : temperature of the isothermal atmosphere
! iter : iteration number
! output_file : output file
! ---------------------------------------------------------

	real(kind=8) :: col_density
	real(kind=8) :: tempc, radius
	integer :: nz, riter, iter
	integer :: parabolic, verbose
	integer :: which_scheme, read_previous_flag
	character(len=80) :: output_file, model_file, collision_type, input_file, filename_previous
	character(len=80) :: atmospheric_file

! ---------------------------------------------------------
! VARIABLES WHICH VARY WITH DEPTH 
! ---------------------------------------------------------
! z : geometricl depth 
! tau : optical depth for each transition
! chil : line opacity
! kappa : continuum opacity 
! chie : electron opacity (scattering)
! Sl : line source function
! Lstar : approximate Lamdba operator for a given transition
! Jbar : mean intensity for a given transition
! ---------------------------------------------------------

	real(kind=8), allocatable :: dz(:), tau(:,:), B(:), chil(:), kappa(:), chie(:), Sl(:), Lstar(:), Jbar(:)
	real(kind=8), allocatable :: flux(:)
	real(kind=8) :: beta_file(1600), alpha_file(1600), tau_file(1600), spline_data(1600)

! ---------------------------------------------------------
! VARIABLES WHICH VARY WITH FREQUENCY
! ---------------------------------------------------------
! perfil : line profile
! In : specific intensity
! ---------------------------------------------------------

   real(kind=8), allocatable :: x_e3(:), w_e3(:)

	
! ---------------------------------------------------------
! VARIABLES WHICH VARY WITH DEPTH AND NUMBER OF TRANSITIONS
! ---------------------------------------------------------
! Jbar_total : mean intensity for all transitions
! Lstar_total : approximate lamdba operator (diagonal) for all transitions
! ---------------------------------------------------------	
		
   real(kind=8), allocatable :: Jbar_total(:), dJbardn(:,:), lstar_total(:), flux_total(:)
	real(kind=8), allocatable :: dSldn(:,:), dJbar_totaldn(:,:), dtaudn(:,:), dchildn(:,:)

! ---------------------------------------------------------
! REST OF VARIABLES
! ---------------------------------------------------------
! precision : precision
! relat_error : maximum relative error
! ---------------------------------------------------------	
   
	real(kind=8) :: precision, relat_error

	
! ---------------------------------------------------------
! ATOMIC MODEL
! ---------------------------------------------------------
! nl : n. levels
! ni : n. ionic states
! nt : n. transitions
! nact : n. active transitions
! nli(nlm) : ion at which the level belongs
! itran(2,nt) : upper and lower levels of the transitions
! abun : abundance
! dlevel(1,nlm) : frequency of each level
! dlevel(2,nlm) : statistical weights
! dion(1,ni) : ionization frequency
! dtran(1,nt) : radiative cross section of the transitions = fij * pi * (e**2) / (me * c)
! dtran(2,nt) : frequency of the transition
! dtran(3,nt) : collisional rate
! dtran(4,nt) : Doppler width 
! vtherm : thermal velocity
! ---------------------------------------------------------

	integer :: nl, ni, nact, nt,nr 
	integer, allocatable :: nli(:), itran(:,:)
	real(kind=8), allocatable :: dlevel(:,:), dion(:,:), dtran(:,:)
	real(kind=8), allocatable :: coll(:,:), doppler(:,:)
	real(kind=8) :: vtherm, einstein_twolevel, collis_twolevel, temp_twolevel
	
! ---------------------------------------------------------
! ATMOSPHERE MODEL
! ---------------------------------------------------------	
! datph(1,nz) : temperature
! datph(2,nz) : electronic density
! ---------------------------------------------------------
	real(kind=8), allocatable :: datph(:,:)
	real(kind=8) :: hydrogen_density, vmicrot
	integer :: npoints
	
! ---------------------------------------------------------
! POPULATIONS
! ---------------------------------------------------------
! nh(nz) : hydrogen abundance
! popl(nz*nl) : population of each level in LTE
! pop(nz*nl) : population of each level
! abund : abundance
! ---------------------------------------------------------

	real(kind=8), allocatable :: nh(:), popl(:), pop(:)
	real(kind=8), allocatable :: abund(:)

	integer :: n_quadr_beta, optically_thin, start_mode, escape_prob_algorithm
	real(kind=8), allocatable :: tau_beta(:), beta_function(:), beta_spline(:), alphap_function(:), alphap_spline(:)
	
end module variables

module functions
use constants
use variables
implicit none

contains

! *********************************************************
! *********************************************************
! MATHEMATICAL ROUTINES
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! Read some lines from a file
! ---------------------------------------------------------
	subroutine lb(unit,nlines)
	integer :: unit, nlines, i
	character(len=128) :: temp
		do i = 1, nlines
			read(unit,FMT='(a128)') temp
		enddo
	end subroutine lb

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! linear interpolation of vector x(:) in y(:)
! ---------------------------------------------------------
   subroutine lin_interpol(xa,ya,x,y)
   real*8, INTENT(IN) :: xa(:), ya(:), x(:)
   real*8, INTENT(INOUT) :: y(:)
   integer :: i, n
   integer :: locat(1), loc

   	n = size(x)

   	do i = 1, n
			loc = 1
			do while (xa(loc) < x(i))
				loc = loc + 1
			enddo   		
   		if (loc == 1) then
				y(i) = ya(loc)
   		else
				y(i) = (ya(loc)-ya(loc-1))/(xa(loc)-xa(loc-1)) * (x(i)-xa(loc-1)) + ya(loc-1)
   		endif
   	enddo

   end subroutine lin_interpol	
	
! ---------------------------------------------------------
! LU decomposition of a matrix
!  INPUT:
!		- a is the matrix to decompose
!		
!  OUTPUT:
!		- a is the LU decomposition of a
!		- indx is a vector that records the row permutation effected by the partial pivoting
!		- d takes values +1/-1 depending on whether the number of row interchanges was odd or even
! ---------------------------------------------------------
	subroutine ludcmp(a,indx,d)
	integer, INTENT(INOUT) :: indx(:)
	real(kind=8), INTENT(INOUT) :: a(:,:), d
	real(kind=8), parameter :: TINY = 1.d-20
	integer :: i, imax, j, k, n
	real(kind=8) :: aamax, dum, sum, vv(size(a,1))
		d = 1.d0
		n = size(a,1)
		
		do i = 1, n
			aamax = 0.d0	
			aamax = maxval(dabs(a(i,:)))
			if (aamax == 0.d0) pause 'Singular matrix in LU decomposition'
			vv(i) = 1.d0 / aamax
		enddo
		
		do j = 1, n
			do i = 1, j-1
				sum = a(i,j)
				do k = 1, i-1
					sum = sum - a(i,k) * a(k,j)
				enddo
				a(i,j) = sum
			enddo
			aamax = 0.d0
			do i = j, n
				sum = a(i,j)
				do k = 1, j-1
					sum = sum - a(i,k) * a(k,j)
				enddo
				a(i,j) = sum
				dum = vv(i) * dabs(sum)
				if (dum >= aamax) then
					imax = i
					aamax = dum
				endif				
			enddo
			if (j /= imax) then
				do k = 1, n
					dum = a(imax,k)
					a(imax,k) = a(j,k)
					a(j,k) = dum
				enddo
				d = -d
				vv(imax) = vv(j)
			endif
			indx(j) = imax
			if (a(j,j) == 0.d0) a(j,j) = TINY
			if (j /= n) then
				dum = 1.d0 / a(j,j)
				do i = j+1, n
					a(i,j) = a(i,j) * dum
				enddo
			endif
		enddo
	
	end subroutine ludcmp

! ---------------------------------------------------------
! Solves the set of equations AX=b where A is the LU decomposition of a matrix
!  INPUT:
!		- a is the LU decomposition of the system matrix
!		- b is the right hand side vector of the system
! 		- indx is the vector returned by ludcmp
!  OUTPUT:
!		- b is the solution of the system
! ---------------------------------------------------------
	subroutine lubksb(a,indx,b)
	real(kind=8), INTENT(IN) :: a(:,:)
	real(kind=8), INTENT(INOUT) :: b(:)
	integer, INTENT(IN) :: indx(:)
	integer :: i, ii, n, j, ll
	real(kind=8) :: sum
		n = size(a,1)
		ii = 0
		do i = 1, n
			ll = indx(i)
			sum = b(ll)
			b(ll) = b(i)
			if (ii /= 0) then
				do j = ii, i-1
					sum = sum - a(i,j) * b(j)
				enddo
			else if (sum /= 0.d0) then
				ii = i
			endif
			b(i) = sum
		enddo
		do i = n, 1, -1
			sum = b(i)
			do j = i+1, n
				sum = sum - a(i,j) * b(j)
			enddo
			b(i) = sum / a(i,i)
		enddo
	end subroutine lubksb
	
! ---------------------------------------------------------
! Resuelve un sistema de ecuaciones
! ---------------------------------------------------------
!     SOLUTION OF A SYSTEM OF LINEAR EQUATIONS -- A*X = X
!
!     INPUT FOR LUSLV
!     A        -LEFT HAND SIDE
!     X        -RIGHT HAND SIDE
!     N        -ORDER OF THE SYSTEM
!     NR       -ACTUAL FIRST DIMENSION OF A
!
!     OUTPUT FOR LUSLV
!     A        -LU FACTORIZATION OF ORIGINAL A
!     X        -SOLUTION OF THE LINEAR EQUATIONS
!
!     ONCE LUSLV HAS BEEN CALLED, OTHER SYSTEMS WITH THE SAME LEFT HAND
!     SIDE BUT DIFFERENT RIGHT HAND SIDES MAY BE EFFICIENTLY SOLVED BY
!     USING RESLV.
!
!     INPUT FOR RESLV
!     A        -LU FACTORIZATION OF LEFT HAND SIDE, PRODUCED BY A PREVIOUS
!              -CALL TO LUSLV
!     X        -RIGHT HAND SIDE
!     N        -ORDER OF THE SYSTEM, MUST BE THE SAME AS WHEN LUSLV WAS CALLED
!     NR       -ACTUAL FIRST DIMENSION OF A
!
!     OUTPUT FOR RESLV
!     X        -SOLUTION OF THE LINEAR EQUATIONS
	SUBROUTINE LUSLV(A,X,N,NR)
	
      IMPLICIT real(kind=8)(A-H,O-Z)
		integer :: n, nr
      DIMENSION A(NR,*),X(*),D(100)
      INTEGER R,P(100), i, k, mr1, j, itemp, jp1
      COMMON /LUCOMM/ D,P
      DO 7 R = 1, N
         DO 1 K = 1, N
            D(K) = A(K,R)
    1    CONTINUE
         MR1 = R - 1
         IF(MR1.LT.1) GO TO 4
         DO 3 J = 1, MR1
            ITEMP = P(J)
            A(J,R) = D(ITEMP)
            D(ITEMP) = D(J)
            JP1 = J + 1
            DO 2 I = JP1, N
               D(I) = D(I) - A(I,J)*A(J,R)
    2       CONTINUE
    3    CONTINUE
    4    DMX = ABS(D(R))
         P(R) = R
         DO 5 I = R, N
            IF(DMX.GT.ABS(D(I))) GO TO 5
            DMX = ABS(D(I))
            P(R) = I
    5    CONTINUE
         ITEMP = P(R)
         A(R,R) = 1.0/D(ITEMP)
         D(ITEMP) = D(R)
         MR1 = R + 1
         IF(MR1.GT.N) GO TO 8
         DO 6 I = MR1, N
            A(I,R) = D(I)*A(R,R)
    6    CONTINUE
    7 CONTINUE
    8 CALL RESLV(A,X,N,NR)
      
      END subroutine luslv
		
		
      SUBROUTINE RESLV(A,X,N,NR)
		
      IMPLICIT real(kind=8)(A-H,O-Z)
		integer :: n, nr
      INTEGER P(100), i, ip1, j, k, kp1, itemp
      real(kind=8) D(100),A(NR,*),X(*)
      COMMON /LUCOMM/ D,P
      DO 9 I = 1, N
         D(I) = X(I)
    9 CONTINUE
      DO 11 I = 1, N
         ITEMP = P(I)
         X(I) = D(ITEMP)
         D(ITEMP) = D(I)
         IP1 = I + 1
         IF(IP1.GT.N) GO TO 12
         DO 10 J = IP1, N
            D(J) = D(J) - A(J,I)*X(I)
   10    CONTINUE
   11 CONTINUE
   12 K = N + 1
      DO 15 I = 1, N
         K = K - 1
         SUM = 0.0
         KP1 = K + 1
         IF(KP1.GT.N) GO TO 14
         DO 13 J = KP1, N
            SUM = SUM + A(K,J)*X(J)
   13    CONTINUE
   14    X(K) = A(K,K)*(X(K)-SUM)
   15 CONTINUE
      END subroutine reslv


      

! ---------------------------------------------------------
! Solve a linear system of equations
! ---------------------------------------------------------	
	subroutine llslv(S,X,N,NR)
	integer, INTENT(IN) :: n, nr
	real(kind=8), INTENT(INOUT) :: S(nr,nr), x(nr)
	real(kind=8) :: sum
	integer :: i, j, k
	
	do i = 1, n
		do j = 1, i
			sum = s(i,j)
			do k = 1, j-1
				sum = sum - S(i,k) * S(j,k)
			enddo !k
			if (j < i) then
				s(i,j) = sum / S(j,j)
			else
				S(i,j) = dsqrt(dabs(sum))
			endif
		enddo !j
	enddo !i
	
	entry llreslv(S,X,N,NR)
	
	do i = 1, n
		sum = x(i)
		do j = 1, i-1
			sum = sum - S(i,j)*x(j)
		enddo !j
		x(i) = sum / S(i,i)
	enddo !i
	do i = n, 1, -1
		sum = x(i)
		do j = n, i+1, -1
			sum = sum - S(j,i)*x(j)
		enddo !j
		x(i) = sum / S(i,i)
	enddo !i
	
	end subroutine llslv
	
! ---------------------------------------------------------
! Return the Voigt profile
! ---------------------------------------------------------
      function voigt(a,vv,j)
      implicit real(kind=8) (a-h,o-z)
      implicit integer (i-n)
      complex :: z
      real(kind=8) :: xdws(28),ydws(28)

      data a0,a1,a2,a3,a4,a5,a6,b0,b1,b2,b3,b4,b5,b6/&
      122.607931777104326,214.382388694706425,181.928533092181549,&
      93.155580458138441,30.180142196210589,5.912626209773153,&
      .564189583562615,122.60793177387535,352.730625110963558,&
      457.334478783897737,348.703917719495792,170.354001821091472,&
      53.992906912940207,10.479857114260399/

      data xdws/.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,&
      3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20./,ydws/&
      9.9335991E-02,1.9475104E-01,2.8263167E-01,3.5994348E-01,&
      4.2443639E-01,4.7476321E-01,5.1050407E-01,5.3210169E-01,&
      5.4072434E-01,5.3807950E-01,5.0727350E-01,4.5650724E-01,&
      3.9993989E-01,3.4677279E-01,3.0134040E-01,1.7827103E-01,&
      1.2934799E-01,1.0213407E-01,8.4542692E-02,7.2180972E-02,&
      6.3000202E-02,5.5905048E-02,5.0253846E-02,4.1812878E-02,&
      3.5806101E-02,3.1311397E-02,2.7820844E-02,2.5031367E-02/

      v=abs(vv)
      IF(A.NE.0) GOTO 1 
      IF(J.NE.0) GOTO 3 
      VOIGT=DEXP(-V*V)   
      RETURN
   3  IF(V.GT.XDWS(1)) GOTO 4   
      D=V*(1.-.66666667*V*V)
      GOTO 8
   4  IF(V.GT.XDWS(28)) GOTO 5  
      K=27  
      DO 7 I=2,27   
      IF(XDWS(I).LT.V) GOTO 7   
      K=I   
      GOTO 6
   7  CONTINUE  
   6  KK=K-1
      KKK=K+1   
      D1=V-XDWS(KK) 
      D2=V-XDWS(K)  
      D3=V-XDWS(KKK)
      D12=XDWS(KK)-XDWS(K)  
      D13=XDWS(KK)-XDWS(KKK)
      D23=XDWS(K)-XDWS(KKK) 
      D=YDWS(KK)*D2*D3/(D12*D13)-YDWS(K)*D1*D3/(D12*D23)+YDWS(KKK)*&
        D1*D2/(D13*D23)
      GOTO 8
   5  Y=.5/V
      D=Y*(1.+Y/V)  
   8  VOIGT=5.641895836E-1*D
   9  IF(VV.LT.0.) VOIGT=-VOIGT 
      RETURN
   1  Z=CMPLX(A,-V) 
      Z=((((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0)/&
      (((((((Z+B6)*Z+B5)*Z+B4)*Z+B3)*Z+B2)*Z+B1)*Z+B0)
      IF(J.NE.0) GOTO 2 
      VOIGT=REAL(Z) 
      RETURN
   2  VOIGT=.5*AIMAG(Z) 
      GOTO 9
      END function voigt
      
! ---------------------------------------------------------
! Return the exponential integral En(x)
! ---------------------------------------------------------      
      FUNCTION expint(n,x)
      INTEGER n,MAXIT
      REAL(kind=8) expint,x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)
      INTEGER i,ii,nm1
      REAL(kind=8) a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
        pause 'bad arguments in expint'
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0.)then
        expint=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.EPS)then
            expint=h*exp(-x)
            return
          endif
11      continue
        pause 'continued fraction failed in expint'
      else
        if(nm1.ne.0)then
          expint=1./nm1
        else
          expint=-log(x)-EULER
        endif
        fact=1.
        do 13 i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-EULER
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          if(abs(del).lt.abs(expint)*EPS) return
13      continue
        pause 'series failed in expint'
      endif
      return
      END function expint


! ---------------------------------------------------------
! Given x(:) and y(:) which tabulate a function and the derivative at the boundary points
! this function returns the second derivative of the spline at each point
! ---------------------------------------------------------
		subroutine splin1(x,y,yp1,ypn,y2)
		real(kind=8), INTENT(IN) :: x(:), y(:), yp1, ypn
		real(kind=8), INTENT(INOUT) :: y2(size(x))
		integer :: n, i, k
		real(kind=8) :: p, qn, sig, un, u(size(x))

			n = size(x)
			
			if (yp1 > .99d30) then
				y2(1) = 0.d0
				u(1) = 0.d0
			else
				y2(1) = -0.5d0
				u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			endif

			do i = 2, n-1
				sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))				
				p = sig * y2(i-1)+2.d0
				y2(i) = (sig-1.d0)/p
				u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/&
					(x(i+1)-x(i-1))-sig*u(i-1))/p
			enddo
			if (ypn > .99d30) then
				qn = 0.d0
				un = 0.d0
			else
				qn = 0.5d0
				un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
			endif
			
			y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

			do k = n-1, 1, -1
				y2(k) = y2(k)*y2(k+1)+u(k)
			enddo

		end subroutine splin1

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! splines of vector x(:) in y(:)
! ---------------------------------------------------------
		subroutine spline(xa,ya,y2a,x,y)
		real(kind=8), INTENT(INOUT) :: y(:)
		real(kind=8), INTENT(IN) :: xa(:), ya(:), x(:)
		real(kind=8) :: y2a(:)
		integer :: n_x, n, i, k, khi, klo
		real(kind=8) :: a, b, h, extrap
			
			n = size(xa)
			n_x = size(x)
			!call splin1(xa,ya,1.d30,1.d30,y2a)						

			do i = 1, n_x					

! Downward extrapolation 
				if (x(i) < xa(1)) then
!					y(i) = ya(1)
					y(i) = ya(1) + (ya(1)-ya(2))/(xa(1)-xa(2)) * (xa(1) - x(i))
				else 

! Upward extrapolation
				if (x(i) > xa(n)) then
!					y(i) = ya(n)
					y(i) = ya(n) + (ya(n)-ya(n-1)) / (xa(n)-xa(n-1)) * (x(i) - xa(n))
				else
! In range
						klo = 1
						khi = n
1						if(khi-klo > 1) then
							k = (khi+klo)/2
							if (xa(k) > x(i)) then
								khi = k
							else
								klo = k
							endif					
							go to 1
						endif

						h = xa(khi)-xa(klo)

						if (h == 0.d0) pause 'bad xa input in spline'
						a = (xa(khi)-x(i))/h
						b = (x(i)-xa(klo))/h

						y(i) = a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0		
					endif
				endif
			enddo

		end subroutine spline
		
!-----------------------------------------------------------------
! Returns the weights (w) and the abscissas (x) for a Gaussian integration using the 
! Gauss-Legendre formula, using n points
!-----------------------------------------------------------------
	subroutine gauleg(x1,x2,x,w,n)
	integer, INTENT(IN) :: n
	real(kind=8), INTENT(IN) :: x1,x2
	real(kind=8), INTENT(INOUT) :: x(n),w(n)
	real(kind=8), parameter :: eps = 3.d-14
	integer :: i,j,m
	real(kind=8) :: p1,p2,p3,pp,xl,xm,z,z1
      
	m=(n+1)/2
   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   do i=1,m
   	z=dcos(3.141592654d0*(i-.25d0)/(n+.5d0))
1  	continue
   	p1=1.d0
   	p2=0.d0
   	do j=1,n
   		p3=p2
      	p2=p1
      	p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
		enddo
   	pp=n*(z*p1-p2)/(z*z-1.d0)
   	z1=z
   	z=z1-p1/pp
  		if(abs(z-z1).gt.EPS)goto 1
   	x(i)=xm-xl*z
   	x(n+1-i)=xm+xl*z
   	w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
   	w(n+1-i)=w(i)
	enddo
	
	end subroutine gauleg


! *********************************************************
! *********************************************************
! MATHEMATICAL ROUTINES
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! Inicialization of all the variables
! ---------------------------------------------------------
   subroutine init

   real(kind=8) :: div, deca, deltatau
   integer :: i

		open(unit=15,file='config.dat',action='read',status='old')
		call lb(15,5)
		read(15,*) input_file
		close(15)
		
		
		open(unit=15,file=input_file,action='read',status='old')		
		
		call lb(15,4)
		read(15,*) verbose
		
		call lb(15,2)
		read(15,*) output_file
		if (verbose == 1) then
			print *, 'Output file : ', output_file
		endif
		      
      call lb(15,2)
      read(15,*) atmospheric_file
      if (verbose == 1) then
			print *, 'Creating atmosphere...'
		endif
		call grid
		
		call lb(15,2)
      read(15,*) model_file
		if (verbose == 1) then
			print *, 'File with the atomic/molecular model : ', model_file
		endif

		call lb(15,2)
      read(15,*) vmicrot
		if (verbose == 1) then
			print *, 'Microturbulent velocity [km/s] : ', vmicrot
		endif
                  
		if (verbose == 1) then
			print *, 'Reading atomic model...'		
		endif
		call atom

		call lb(15,2)
      read(15,*) escape_prob_algorithm
		if (verbose == 1) then
			print *, 'Escape probability (0-> KROLIK & McKEE, 1-> EXACT) : ', escape_prob_algorithm
		endif
						
		
		call lb(15,2)
		read(15,*) read_previous_flag
		if (verbose == 1) then
			print *, 'Read previous calculation : ', read_previous_flag
		endif
		
		call lb(15,2)
		read(15,*) filename_previous
		if (verbose == 1 .and. read_previous_flag == 1) then
			print *, 'File with the previous calculation : ', filename_previous
		endif
		
		call lb(15,2)
		read(15,*) start_mode
		if (verbose == 1) then
			print *, 'Starting solution : ', start_mode
		endif
		
		call lb(15,2)
		read(15,*) which_scheme
		if (verbose == 1) then
			print *, 'Numerical scheme : ', which_scheme
		endif
		
		call lb(15,2)
		read(15,*) precision
		if (verbose == 1) then
			print *, 'Precision of the solution : ', precision
		endif
		
		close(15)
								
! Allocate memory for all the variables depending on the size of the atmosphere
		allocate(B(nz))
		allocate(chil(nz))
		allocate(kappa(nz))
		allocate(chie(nz))
		allocate(Sl(nz))
		allocate(Lstar(nz))
		allocate(Jbar(nz))
		allocate(flux(nz))
		
! Calculate the weigths for the angular and frequency integration
      call ang_freq_weights
		
		
		allocate(Jbar_total(nz*nact))
		allocate(dJbar_totaldn(nz*nact,nz*nl))
		allocate(dJbardn(nz,nz*nl))
		allocate(lstar_total(nz*nact))
		allocate(flux_total(nz*nact))		
		allocate(pop(nl*nz))
		allocate(popl(nl*nz))
		allocate(tau(nact,0:nz))
		
		allocate(dSldn(nz,nl*nz))
		allocate(dchildn(nz,nl*nz))
		allocate(dtaudn(0:nz,nl*nz))
		
		
		kappa = 0.d0
		chie = 0.d0

! Initialize in LTE
		if (verbose == 1) then
			print *, 'Initializing in LTE...'
		endif
		call poplte
		pop(1:nl*nz) = popl(1:nl*nz)		
				
		if (verbose == 1) then
			print *, 'End of initialization'		
		endif
					
   end subroutine init         

      
! *********************************************************
! *********************************************************
! ROUTINES FOR THE ATOMIC MODEL
! *********************************************************
! *********************************************************

!----------------------------------------------------------
! Define the atomic/molecular model
!----------------------------------------------------------
	subroutine atom
	integer :: i, index, up, low, n_collis, j
	real(kind=8) :: masa, col(1), tem(1)
	real(kind=8), allocatable :: values(:), temperatures(:), collis_spline(:)
	
!----------------------------------------------------------
!     DEFINE ATOMIC DATA
! 
! OUTPUT: 
!         INTEGER:
!         NL=totalnumber of levels , NI=Number of ions
!         NT=Number of transitions
!         NACT=Number of active radiative transitions
!         NR=Number of active+passive radiative transitions
!         NLI(NL)=ion to which each level belongs
!         ITRAN(2,nt)=Up/down levels in each transition
!         REAL*8:
!         ABUND= Abundance 
!         DLEVEL(1,nl)=Frecuency of every level
!       
!         DLEVEL(2,nl)=g=Pesos estadisticos
!
!         DION(1,ni)=XI=Ionization frequency
!         DTRAN(1,nt)=Radiative cross section=fij*pi*(e**2)/me/c
!         DTRAN(2,nt)=frecuency of the transition
!         DTRAN(3,nt)=Collisional cross section
!         DTRAN(4,nt)=Radiation Temperature OR Doppler Width
!----------------------------------------------------------

!
! 
! Number of ionization states is always one
		ni=1
		
		open(unit=12,file=model_file,action='read',status='old')
		
		call lb(12,5)
		read(12,*) masa
		
		call lb(12,2)
		read(12,*) masa
		print *, 'Mass : ', masa
						
		call lb(12,2)
		read(12,*) nl
		print *, 'Number of levels : ', nl
		
		call lb(12,2)
		read(12,*) nt
		print *, 'Number of transitions : ', nt
		
		call lb(12,2)
		read(12,*) nr
		print *, 'Number of radiative transitions : ', nr
		
		call lb(12,2)
		read(12,*) nact
		print *, 'Number of active transitions : ', nact

! Allocate memory
		allocate(nli(nl))
		allocate(itran(2,nt))
		allocate(dlevel(2,nl))
		allocate(dion(1,ni))
		allocate(dtran(4,nt))
		allocate(coll(nz,nt))
		allocate(doppler(nz,nt))
		
		dion(1,1)=0.d0
		
! Ion to which the levels belong
      do i=1,nl
      	nli(i)=1
		enddo
		
		call lb(12,2)
		do i = 1, nl
			read(12,*) index, dlevel(1,i), dlevel(2,i)
			dlevel(1,i) = dlevel(1,i) * PC
		enddo
		
		call lb(12,2)
		read(12,*) collision_type
		
		call lb(12,2)
		read(12,*) n_collis
								
		if (trim(adjustl(collision_type)) == 'INTERPOL') then
			print *, 'Using interpolated collisions'
			allocate(values(n_collis))
			allocate(temperatures(n_collis))
			allocate(collis_spline(n_collis))
			call lb(12,2)
			read(12,*) (temperatures(i),i=1,n_collis)
			call lb(12,2)
		else
			print *, 'Using fit to collisional rates'
			allocate(values(2))
			call lb(12,5)
		endif
		
		print *, 'Einstein'
! Now the radiative transitions and transform Einstein coefficient:
! c^2/(8*pi*nu^2)*(gu/gl)*Aul = h*nu/(4*pi)*Blu
		do i = 1, nr
			read(12,*) index, itran(1,i), itran(2,i), dtran(2,i), dtran(1,i)
			einstein_twolevel = dtran(1,i)
			up = itran(1,i)
			low = itran(2,i)
			dtran(1,i) = dtran(1,i) / (pi8c2*dtran(2,i)**2)*(dlevel(2,up)/dlevel(2,low))
		enddo
		print *, 'Doppler'		
! Doppler Line Widths
      do i = 1,nr
      	do j = 1, nz
      		vtherm = sqrt(2.d0*PK*datph(1,i)/(masa*UMA))
				if (vtherm/1.d5 < vmicrot) then
					vtherm = vmicrot * 1.d5					
				endif
							
      		doppler(j,i) = dtran(2,i) * vtherm / PC      		
      	enddo
		enddo
		
		call lb(12,2)

! And now collisional rates. The first nr transitions have to be equal to the radiative transitions
		if (trim(adjustl(collision_type)) == 'FIT') then
			do i = 1, nt				
				read(12,*) index, itran(1,i), itran(2,i), (values(j),j=1,2)				
				do j = 1, nz
					coll(j,i) = nh(j) * values(1) * datph(1,j)**values(2)
				enddo
				collis_twolevel = coll(1,1)
			enddo
		endif
			
		if (trim(adjustl(collision_type)) == 'INTERPOL') then
			do i = 1, nt
				read(12,*) index, itran(1,i), itran(2,i), (values(j),j=1,n_collis)
				do j = 1, nz
					tem(1) = datph(1,j)
					call lin_interpol(temperatures,values,tem,col)
					coll(j,i) = nh(j) * col(1)
				enddo				
			enddo
		endif
				
		if (collision_type == 'INTERPOL') then
			deallocate(values)
			deallocate(temperatures)
			deallocate(collis_spline)
		else
			deallocate(values)
		endif

		end subroutine atom

! ---------------------------------------------------------
! Initialize the integration weights in angle and frequency
! ---------------------------------------------------------
	subroutine ang_freq_weights
	real(kind=8) :: div
	real(kind=8), allocatable :: x_quadr(:), w_quadr(:)
	integer :: i
				      
      n_quadr_beta = 80
      allocate(x_e3(n_quadr_beta))
		allocate(w_e3(n_quadr_beta))
      call gauleg(-7.d0,7.d0,x_e3,w_e3,n_quadr_beta)
            		
	end subroutine ang_freq_weights

		
!-----------------------------------------------------------------
! Returns the number of points and sets the z axis
! tau(i) = tau(i-1) * taufac
!-----------------------------------------------------------------		
	subroutine grid
	real(kind=8) :: deltaz, ztemp, tautemp, tau0, t1, t2, t3, t4, t11
	integer :: i
		
! First calculate the number of zones in the atmosphere
		open(unit=12,file=atmospheric_file,status='old',action='read')
		call lb(12,3)
		read(12,*) npoints
		
		nz = min(npoints,n_ptos)
		allocate(dz(nz))
		allocate(datph(2,nz))
		allocate(nh(nz))
		allocate(abund(nz))
		
		print *, 'Number of zones : ', nz		
		do i = 1, nz
			read(12,*) t1, t2, t3, t4
			dz(i) = t1
			datph(1,i) = t2
			datph(2,i) = 1.d0
			nh(i) = t3
			abund(i) = t4 / t3
		enddo
		
		close(12)
						
	end subroutine grid		
				
! *********************************************************
! *********************************************************
! GENERAL ROUTINES
! *********************************************************
! *********************************************************
	
! ---------------------------------------------------------
! Put the LTE populations
! ---------------------------------------------------------	
	subroutine poplte
	real(kind=8) :: u(ni), fi(ni), fl(nl), sum ,kte, tp, fac, tot
	integer :: ip, i, ipl
	
		do ip = 1, nz
			ipl = nl * (ip-1)
			tp = datph(1, ip)
			kte = PHK / tp
			
			u = 0.d0

! Partition function
			do i = 1, nl
				fl(i) = dlevel(2,i) * dexp(-dlevel(1,i) * kte)  !g*exp(-h*nu/(k*T))
				u(nli(i)) = u(nli(i)) + fl(i)
			enddo
			
! If we are in a multi-ionic atom
			if (ni > 1) then
				
				do i = 1, nl
					fl(i) = fl(i) / u(nli(i))
				enddo
				
				fac = datph(2,ip) * ci / (tp * dsqrt(tp))
				
				do i = 1, ni-1
					fi(i) = fac * u(i) / u(i+1) * dexp(dion(1,i+1)*kte)
				enddo
				
				fi(ni) = 1.d0
				tot = fi(ni)
				
				do i = ni - 1, 1, -1
					fi(i) = fi(i) * fi(i+1)
					tot = tot + fi(i)
				enddo
				tot = abund(ip) *  nh(ip) / tot
				do i = 1, nl
					popl(i+ipl) = fl(i) * fi(nli(i)) * tot
				enddo
			else
				tot = abund(ip) * nh(ip) / u(1)
				do i = 1, nl
					popl(i + ipl) = fl(i) * tot
				enddo
			endif						
			
		enddo 
	end subroutine poplte

! ---------------------------------------------------------
! Read data from a previous file
! ---------------------------------------------------------
	subroutine read_previous(filename_previous)
	character(len=80) :: filename_previous
	integer :: nl_old, nz_old, ip, ipl, i
	real(kind=8), allocatable :: z_old(:), z_new(:), pop_old(:), pop1(:), pop2(:)
	real(kind=8) :: temp, r_old, deltaz_new, deltaz_old
		open(unit=20,file=filename_previous,status='old',action='read')
		
		call lb(20,6)
		read(20,*) nl_old
		call lb(20,1)
		read(20,*) nz_old
		print *, 'Previous/new number of levels : ', nl_old, nl
		print *, 'Previous/new number of points : ', nz_old, nz
		
		call lb(20,7)
		read(20,*) r_old
		
		print *, 'Previous/new radius : ', r_old, radius
				
		deltaz_old = r_old / nz_old
		deltaz_new = radius / nz
		
		print *, 'Previous/new dz : ', deltaz_old, deltaz_new
		
		allocate(z_old(nz_old))
		allocate(z_new(nz))
		allocate(pop_old(nz_old*nl_old))
		allocate(pop1(nz_old))
		allocate(pop2(nz))
		
		call lb(20,3)
		call lb(20,nz_old)
		
		z_old(1) = 0.d0
		z_new(1) = 0.d0
		
		do ip = 2, nz_old
			z_old(ip) = z_old(ip-1) + deltaz_old
		enddo
		
		do ip = 2, nz
			z_new(ip) = z_new(ip-1) + deltaz_new
		enddo
		
		call lb(20,3)

! Read the old populations		
		do ip = 1, nz_old
			ipl = nl_old * (ip-1)
			read(20,*) temp, (pop_old(i+ipl), i = 1, nl_old)
		enddo

! Now do the interpolation into the new z axis		
		do i = 1, nl_old
			do ip = 1, nz_old
				ipl = nl_old*(ip-1)
				pop1(ip) = pop_old(i+ipl)
			enddo
			call lin_interpol(z_old,pop1,z_new,pop2)
			do ip = 1, nz
				ipl = nl*(ip-1)
				pop(i+ipl) = pop2(ip)
			enddo
		enddo
		
		close(20)
		
		deallocate(z_old)
		deallocate(z_new)
		deallocate(pop_old)
		deallocate(pop1)
		deallocate(pop2)
		
	end subroutine read_previous
	
!-----------------------------------------------------------------	
! Write the results
!-----------------------------------------------------------------
	subroutine write_results
	integer :: i, ip, ipl, ipt, it, up, low
	real(kind=8) :: sig, glu, acon, agnus, chilp, chilm, snlte, slte, sb(nz)
	real(kind=8) :: source(nz), epsil_equiv, suma, total, tau_estimation
	real(kind=8), allocatable :: freq_axis(:), out_flux(:,:)
	
		open (UNIT=3,FILE=trim(adjustl(output_file))//'.out',STATUS='replace',ACTION='write')	
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                GENERAL DATA                         '
		write(3,*) '-----------------------------------------------------'
		write(3,*) 'N. active transitions'
		write(3,*) nact
		write(3,*) 'N. levels'
		write(3,*) nl
		write(3,*) 'Equivalent epsilon (assuming T is constant)'
		temp_twolevel = datph(1,1)
		epsil_equiv = collis_twolevel / einstein_twolevel * &
			(1.d0 - dexp(-PHK * dtran(2,1) / temp_twolevel) )		
		write(3,*) epsil_equiv / (1.d0+epsil_equiv)
		write(3,*) 'Estimated total optical depth'
		
		tau_estimation = einstein_twolevel * PC**3 /&
			(8.d0*PI**(1.5)*dtran(2,1)**3*vmicrot*1.d5)
		tau_estimation = tau_estimation * &
			(dlevel(2,2) / dlevel(2,1) * pop(1) - pop(2)) * sum(dz)
		
		write(3,*) tau_estimation
		write(3,*) 'N. zones'
		write(3,*) nz
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '        DeltaZ GRID and TAU for 1st transition       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, nz
			write(3,FMT='(2E12.5)') dz(ip), tau(1,ip)
		enddo
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '             FINAL POPULATIONS                       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, nz
			ipl = nl * (ip-1)
			write(3,FMT='(I4, (1P5E12.5))') ip, (pop(i+ipl), i = 1, nl)
		enddo
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '             TOTAL POPULATIONS                       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, nz
			ipl = nl * (ip-1)
			suma = 0.d0
			do i = 1, nl
				suma = suma + pop(i+ipl)
			enddo
			write(3,FMT='(I4, (1P5E12.5))') ip, suma
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '        FINAL DEPARTURE COEFFICIENTS                 '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, nz
			ipl = nl * (ip-1)
			write(3,FMT='(I4, (1P5E12.5))') ip, (pop(i+ipl)/popl(i+ipl), i = 1, nl)
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                    FINAL LSTAR                      '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, nz
			ipt = nact * (ip-1)
			write(3,FMT='(I4, 1X, 5(1X,1PE15.8))') ip, (lstar_total(i+ipt), i = 1, nact)
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                    FINAL JBAR                       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, nz
			ipt = nact * (ip-1)
			write(3,FMT='(I4, 1X, 5(1X,1PE15.8))') ip, (Jbar_total(i+ipt), i = 1, nact)
		enddo
					
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                   FINAL EMISSIVITY                  '
		write(3,*) '-----------------------------------------------------'
				
		ipt = nact * (nz-1)
		write(3,FMT='(5(1X,1PE15.8))') (flux_total(i+ipt), i = 1, nact)
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '      FINAL RELATIVE EMISSIVITY PERCENTAGE           '
		write(3,*) '-----------------------------------------------------'
		
		ipt = nact * (nz-1)
		total = 0.d0
		do i = 1, nact
			total = total + flux_total(i+ipt)
		enddo		
		write(3,FMT='(5(1X,F6.2))') (flux_total(i+ipt) / total * 100.d0 , i = 1, nact)
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                  SOURCE FUNCTIONS                   '
		write(3,*) '-----------------------------------------------------'
		
		do it = 1, nact			
			up	= itran(1,it)
			low = itran(2,it)
			sig = dtran(1,it) !/ dtran(4,it)
			glu = dlevel(2,low) / dlevel(2,up)
			acon = TWOHC2 * dtran(2,it)**3   !2*h*nu^3/c^2
			agnus = acon * glu
			
			chilm = 0.d0
			do i = 1, nz
				ipl = nl * (i-1)
				snlte = agnus * pop(up+ipl) / (pop(low+ipl) - glu*pop(up+ipl))
				slte = agnus * popl(up+ipl) / (popl(low+ipl) - glu*popl(up+ipl))
				sb(i) = snlte / slte
				source(i) = snlte
			enddo
			
			write(3,*) '*****************************************************'
			write(3,FMT='(A,I4)') ' Transition : ', it
			write(3,FMT='(A,I4,2X,I4)') '  Levels : ', up, low
			write(3,FMT='(A,1PE10.4)') '  Planck function : ', slte
			write(3,FMT='(A,1PE10.4)') '  Wavelength [micron] : ', PC / dtran(2,it) * 1.d4
			write(3,*) '     tau_center         S/B          Tex(K)'
			do i = 1, nz
				write(3,FMT='(4X,1PE10.4,5X,1PE10.4,5X,1PE10.4)') tau(it,i)/sqrt(PI), sb(i), &
					(PH*dtran(2,it))/PK / log(1.d0+acon/source(i))
				
			enddo
			
			write(3,*)
		enddo
		
		close(3)
		
		open (UNIT=3,FILE=trim(adjustl(output_file))//'.flux',STATUS='replace',ACTION='write')	
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                    FINAL FLUX                       '
		write(3,*) '-----------------------------------------------------'
		
		allocate(freq_axis(100))
		allocate(out_flux(nr,100))
		out_flux = 0.d0
		do i = 1, 100
			freq_axis(i) = -6.d0 + (i-1.d0) / 99.d0 * 12.d0
		enddo
		
		call calcflux_cep(pop, freq_axis, out_flux)
		do ip = 1, 100
			write(3,FMT='(F10.6,1X, 5(1X,1PE15.8))') freq_axis(ip), (out_flux(i,ip),i=1,nr)
		enddo
		deallocate(freq_axis)
		deallocate(out_flux)
		close(3)

	end subroutine write_results


! *********************************************************
! *********************************************************
! RATE EQUATIONS ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Put the collisional rates in the matrix which is passed as a parameter for the point ip in the atmosphere
! C(i,j) :  downward collisiona rate (i>j)
! In the subroutine which writes the rate equations, this matrix is transposed : C(low,up) goes in A(up,low)
!-----------------------------------------------------------------
	subroutine collis(ip,C)
	integer, INTENT(IN) :: ip
	real(kind=8), INTENT(INOUT) :: C(nl, nl)
	integer :: it, up, low

		do it = 1, nt
			up = itran(1,it)
			low = itran(2,it)			

			C(up,low) = coll(ip,it)
		enddo
		
	end subroutine collis


!-----------------------------------------------------------------
! Solves the rate equations for the point ip
! It is based on the preconditioning scheme proposed by Rybicki & Hummer (1992)
!-----------------------------------------------------------------
	subroutine rate_eq(ip, dxmax)
	integer, INTENT(IN) :: ip
	real(kind=8), INTENT(INOUT) :: dxmax
	integer :: low, up, ipl, ipt, i, j, il, it
	real(kind=8) :: A(nl, nl), B(nl), cul, acon, blu, glu, djbar, poptot, tol
	integer, allocatable :: indx(:)

! Offsets for the storing
		ipl = nl*(ip-1)
		ipt = nr*(ip-1)
		
		A = 0.d0
		
! Upward collisional rates A(low,up)
		call collis(ip,A)		

! Calculate the downward rates using detailed balance and interchange positions
		do i = 1, nl
			do j = 1, i-1
				cul = A(i,j)
! Rate into de j-th level				
				A(j,i) = cul
! Rate into de i-th level				
				A(i,j) = popl(i+ipl) / popl(j+ipl) * cul
			enddo
		enddo
		
! Active radiative rates (preconditioning)
		do it = 1, nact
		
			up = itran(1,it)
			low = itran(2,it)
			acon = TWOHC2*dtran(2,it)**3  !acon = (2*h*nu**3)/c**2
			blu = PI4H * dtran(1,it) / dtran(2,it)
			glu = dlevel(2,low) / dlevel(2,up)
			
			djbar = Jbar_total(it+ipt) - Lstar_total(it+ipt) * acon * pop(up+ipl) * glu /&
					 ( pop(low+ipl) - glu * pop(up+ipl) )
			
			A(low,up) = A(low,up) + glu * blu * (acon * (1.d0 - Lstar_total(it+ipt)) + djbar )
			A(up,low) = A(up,low) + blu * djbar
			
		enddo		
		
! Background radiative rates
		do it = nact + 1, nr
		
			up = itran(1,it)
			low = itran(2,it)
			acon = TWOHC2 * dtran(2,it)**3  !acon = (2*h*nu**3)/c**2
			blu = PI4H * dtran(1,it) / dtran(2,it)
			glu = dlevel(2,low) / dlevel(2,up)
			
			A(low,up) = A(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
			A(up,low) = A(up,low) + blu * Jbar_total(it+ipt)
			
		enddo		
		
! The system is homogeneous with null determinant. We have to substitute one of the equations by a closure
! relation. The sum over one column gives the total output rate from level i
		do i = 1, nl

			do j = 1, i-1
				A(i,i) = A(i,i) - A(j,i)
			enddo
			do j = i+1, nl
				A(i,i) = A(i,i) - A(j,i)
			enddo
			
		enddo		
		
! Conservation equation. We substitute the last equation
		B = 0.d0
		A(nl, 1:nl) = 1.d0
		
		poptot = abund(ip) * nh(ip)
		B(nl) = poptot
		
! Solve the linear system
		
		allocate(indx(nl))
		tol = 1.d-10
		call ludcmp(A,indx,tol)
		call lubksb(A,indx,B)
		deallocate(indx)
		

! Maximum relative change
		do il = 1, nl
			dxmax = max(dxmax, dabs( B(il) - pop(il+ipl) ) / dabs(B(il)) )
			pop(il+ipl) = B(il)
		enddo
	
	end subroutine rate_eq
		
!-----------------------------------------------------------------
! Do the population correction for all the points in the atmosphere
!-----------------------------------------------------------------
	subroutine correct_populations
	real(kind=8) :: relat_error_p
	integer :: ip
	
		relat_error_p = relat_error
		relat_error = 0.d0
		
		do ip = 1, nz
			call rate_eq(ip,relat_error)			
		enddo	
		
	end subroutine correct_populations
	   
!---------------------------------------------------------
!---------------------------------------------------------
! CEP routines
!---------------------------------------------------------
!---------------------------------------------------------
	
!-----------------------------------------------------------------
! Pre-calculate the optical depths and beta function
!-----------------------------------------------------------------	
	subroutine precalculate_beta
	integer :: i
		allocate(tau_beta(1000))
		allocate(beta_function(1000))
		allocate(beta_spline(1000))
		allocate(alphap_function(1000))
		allocate(alphap_spline(1000))
		
		print *, 'Calculating grid for beta and alpha function in the log range:', log10(1.d-30), &
			log10(1.d20)
		do i = 1, 1000
			tau_beta(i) = (i-1.d0) / 999.d0 * (log10(1.d20)-log10(1.d-30)) + log10(1.d-30)
		enddo
		
		do i = 1, 1000			
			beta_function(i) = beta(10.d0**tau_beta(i))
			alphap_function(i) = alphap(10.d0**tau_beta(i))
		enddo
					
		call splin1(tau_beta,beta_function,1.d30,1.d30,beta_spline)
		
	end subroutine precalculate_beta
   
!---------------------------------------------------------
! This function calculates the beta function
!---------------------------------------------------------
	function beta(tau_in)
   real(kind=8) :: beta, tau_in, salida, tau, paso, coef
   real(kind=8) :: d, b, q, dbdx, x(1), y(1)
   integer :: i, k
   
! KROLIK & McKEE APPROXIMATION   
   	coef = 1.d0/(6.d0*sqrt(3.d0))
   	tau = tau_in / dsqrt(PI)
   	if (tau < 1.d-4) then 
			salida = 1.d0
		else
			if (tau >= 3.41d0) then
   			salida = 1.d0 / (dsqrt(PI)*tau) * ( dsqrt(log(tau)) + 0.25d0 / dsqrt(log(tau)) + 0.14d0)
   		else
				salida = 1.d0 - 0.415d0*tau + 0.355d0*tau*log(tau)
   			dbdx = 0.355 -0.415 + 0.355*log(tau)
   			k = 1
   			d = tau * coef
   			b = 2.d0 * d
   			q = 1.d0
          	do while (q > 1.d-3)          
            	salida = salida - d * tau
            	dbdx = dbdx - b
             	k = k + 1             	
             	d = -d * tau * sqrt((k+1.d0)/(k+2.d0))*(k-1.d0)/(k*(k+2.d0))
             	b = (k+1.d0) * d
             	q = abs(b/dbdx)
            enddo	    			
	   	endif
		endif
   	beta = salida
   	
! Calculation of the beta function by integrating the E3 exponential integral function
   	salida = 0.d0
   	if (tau_in > 1.d-5) then
   		do i = 1, n_quadr_beta
   			salida = salida + w_e3(i)*(0.5d0-expint(3,tau_in*exp(-x_e3(i)**2)/sqrt(PI))) / tau_in				
   		enddo
   	endif
   	
   	if (tau_in <= 1.d-5) then
			salida = 0.99999999d0 !1.d0
   	endif
   	beta = salida
   	   	
   end function beta
	
!---------------------------------------------------------
! This function calculates the beta function
!---------------------------------------------------------
	function alphap(tau_in)
   real(kind=8) :: alphap, tau_in, salida, tau, paso, coef, profile
   real(kind=8) :: d, b, q, dbdx, x(1), y(1)
   integer :: i, k
      	
! Calculation of the derivative of the alpha function by integrating the E2 exponential integral function
   	salida = 0.d0   	
   	do i = 1, n_quadr_beta
			profile = exp(-x_e3(i)**2)/sqrt(PI)
   		salida = salida + w_e3(i)*profile*expint(2,tau_in*profile)
   	enddo
   
   	alphap = salida
   	   	
   end function alphap
   
!---------------------------------------------------------
! This function calculates the beta function by interpolating with a spline routine
!---------------------------------------------------------
	function beta2(tau_in)
   real(kind=8) :: beta2, tau_in, salida, tau, paso, coef
   real(kind=8) :: d, b, q, dbdx, x(1), y(1)
   integer :: i, k
      	
! KROLIK & McKEE APPROXIMATION   
   	if (escape_prob_algorithm == 0) then

   		coef = 1.d0/(6.d0*sqrt(3.d0))
   		tau = tau_in / dsqrt(PI)
			
			if (tau < 0.d0) then
				if (tau < -60.d0) tau = -60.d0
				beta2 = (1.d0-exp(-tau)) / tau
				return
			endif
			
   		if (tau < 1.d-4) then 
				salida = 1.d0
			else
				if (tau >= 3.41d0) then
   				salida = 1.d0 / (dsqrt(PI)*tau) * ( dsqrt(log(tau)) + 0.25d0 / dsqrt(log(tau)) + 0.14d0)
   			else
					salida = 1.d0 - 0.415d0*tau + 0.355d0*tau*log(tau)
   				dbdx = 0.355 -0.415 + 0.355*log(tau)
   				k = 1
   				d = tau * coef
   				b = 2.d0 * d
   				q = 1.d0
          		do while (q > 1.d-3)          
            		salida = salida - d * tau
            		dbdx = dbdx - b
             		k = k + 1             	
             		d = -d * tau * sqrt((k+1.d0)/(k+2.d0))*(k-1.d0)/(k*(k+2.d0))
             		b = (k+1.d0) * d
             		q = abs(b/dbdx)
            	enddo	    			
	   		endif
			endif
   		beta2 = salida

! Interpolation on a table calculated with the exact expression
		else
		
			if (tau_in == 0.d0) then
				beta2 = 1.d0
				return
   		endif   	

			if (tau_in < 0.d0) then
				print *, 'NEGATIVE tau in BETA ', tau_in
				if (tau_in < -30.d0) tau_in = -30.d0
				beta2 = (1.d0-exp(-tau_in)) / tau_in
				tau_in = -tau_in
	!			return
			endif


   		x(1) = log10(tau_in)
   		call spline(tau_beta,beta_function,beta_spline,x,y)
   		salida = y(1)   	   	
   		beta2 = salida
			
		endif
   	   	   	
   end function beta2 
	
!---------------------------------------------------------
! This function calculates the derivative of the alpha function by interpolating with a spline routine
!---------------------------------------------------------
	function alphap2(tau_in)
   real(kind=8) :: alphap2, tau_in, salida, tau, paso, coef
   real(kind=8) :: d, b, q, dbdx, x(1), y(1)
   integer :: i, k
      	
   	
! KROLIK & McKEE APPROXIMATION   
   	if (escape_prob_algorithm == 0) then

   		coef = 1.d0/(6.d0*sqrt(3.d0))
   		tau = tau_in / dsqrt(PI)
			
			if (tau < 0.d0) then
				if (tau < -60.d0) tau = -60.d0
				salida = (1.d0-exp(-tau)) / tau
				dbdx = (-1.d0+(1.d0+tau)*exp(-tau)) / tau**2
				alphap2 = salida + tau * dbdx
				return
			endif
			
   		if (tau < 1.d-4) then 
				salida = 1.d0
			else
				if (tau >= 3.41d0) then
   				salida = 1.d0 / (dsqrt(PI)*tau) * ( dsqrt(log(tau)) + 0.25d0 / dsqrt(log(tau)) + 0.14d0)
   			else
					salida = 1.d0 - 0.415d0*tau + 0.355d0*tau*log(tau)
   				dbdx = 0.355 -0.415 + 0.355*log(tau)
   				k = 1
   				d = tau * coef
   				b = 2.d0 * d
   				q = 1.d0
          		do while (q > 1.d-3)          
            		salida = salida - d * tau
            		dbdx = dbdx - b
             		k = k + 1             	
             		d = -d * tau * sqrt((k+1.d0)/(k+2.d0))*(k-1.d0)/(k*(k+2.d0))
             		b = (k+1.d0) * d
             		q = abs(b/dbdx)
            	enddo	    			
	   		endif
			endif
   		alphap2 = salida + tau * dbdx
			
		else
			
			if (tau_in == 0.d0) then
				alphap2 = 1.d0
				return
			endif

			if (tau_in < 0.d0) then
				print *, 'NEGATIVE tau in ALPHA ', tau_in
				if (tau_in < -30.d0) tau_in = -30.d0
				alphap2 = (1.d0-exp(-tau_in)) / tau_in + tau_in * (-1.d0+(1.d0+tau_in)*exp(-tau_in)) / tau_in**2
				tau_in = -tau_in	
	!			return
			endif

   		x(1) = log10(tau_in)
   		call spline(tau_beta,alphap_function,alphap_spline,x,y)
   		salida = y(1)   	   	
   		alphap2 = salida
			
		endif
   	   	   	
   end function alphap2	              	
   
!---------------------------------------------------------
! This subroutine calculates the flux for every transition using the coupled escape probability method
!---------------------------------------------------------
   subroutine calcflux_cep(x, freq_axis, flux_out)
	real(kind=8) :: x(:), freq_axis(:), flux_out(:,:)
	integer :: up, low, it, ip, ipl, ipt, ip1, npmip1, i
	real(kind=8) :: sig2, sig, glu, acon, chim, chip, tau0, deltaz, chilp
	real(kind=8), allocatable :: line_profile(:)
	
		allocate(line_profile(size(freq_axis)))

! Calculate line absorption coefficient, line source function and optical depths
		do it = 1, nr			
			up = itran(1,it)
			low = itran(2,it)
			sig2 = dtran(1,it) !/ dtran(4,it) !sig=(h*nu*Blu)/(4*pi*Dnudoppler)
			glu = dlevel(2,low) / dlevel(2,up)
			acon = TWOHC2 * dtran(2,it)**3  !acon = (2*h*nu**3)/c**2

			tau0 = 0.d0
			chip = 0.d0
			chim = 0.d0
			ip1 = 1
			tau(it,0) = 0.d0
			do ip = 1, nz
				sig = sig2 / doppler(ip,it)
				ipl = nl * (ip-1)

				chil(ip) = sig * (x(low+ipl) - glu * x(up+ipl))
				Sl(ip) = acon * glu * x(up+ipl) / (x(low+ipl) - glu * x(up+ipl))

         	if (ip == 1) then
					tau(it,ip) = chil(ip) * dz(ip)					
         	else
					tau(it,ip) = tau(it,ip-1) + chil(ip) * dz(ip)
         	endif

			enddo

! Calculate the flux at each wavelength (Eq. (77) in the notes)
			line_profile = exp(-freq_axis**2) / sqrt(PI)
			do i = 1, size(freq_axis)
				do ip = 1, nz
					flux_out(it,i) = flux_out(it,i) + 2.d0 * PI * Sl(ip) * &
						(expint(3,(tau(it,nz)-tau(it,ip))*line_profile(i)) - &
						expint(3,(tau(it,ip)-tau(it,0))*line_profile(i)) + &
						expint(3,(tau(it,ip-1)-tau(it,0))*line_profile(i)) - &
						expint(3,(tau(it,nz)-tau(it,ip-1))*line_profile(i)))
				enddo
			enddo
		enddo
		
		deallocate(line_profile)
  
   end subroutine calcflux_cep
	
!---------------------------------------------------------
! This subroutine calculates Jbar for every transition using the coupled escape probability method
!---------------------------------------------------------
   subroutine calcJbar_Lstar_cep(x)
	real(kind=8) :: x(:)
	integer :: up, low, it, ip, ipl, ipt, ip1, npmip1, i
	real(kind=8) :: sig2, sig, glu, acon, chim, chip, tau0, deltaz, chilp
	
! Calculate line absorption coefficient, line source function and optical depths
! Calculate also their derivatives with respect to the level populations
		do it = 1, nr
			dSldn = 0.d0
			dchildn = 0.d0
			up = itran(1,it)
			low = itran(2,it)
			sig2 = dtran(1,it) !/ dtran(4,it) !sig=(h*nu*Blu)/(4*pi*Dnudoppler)
			glu = dlevel(2,low) / dlevel(2,up)
			acon = TWOHC2 * dtran(2,it)**3  !acon = (2*h*nu**3)/c**2

			tau(it,0) = 0.d0
			dtaudn(0,:) = 0.d0
			do ip = 1, nz
				ipl = nl * (ip-1)
	
				sig = sig2 / doppler(ip,it)
				chil(ip) = sig * (x(low+ipl) - glu * x(up+ipl))
				Sl(ip) = acon * glu * x(up+ipl) / (x(low+ipl) - glu * x(up+ipl))
				
				dSldn(ip,up+ipl) = acon*glu / (x(low+ipl) - glu * x(up+ipl))**2 * x(low+ipl)
				dSldn(ip,low+ipl) = -acon*glu / (x(low+ipl) - glu * x(up+ipl))**2 * x(up+ipl)
				
				dchildn(ip,up+ipl) = -sig * glu
				dchildn(ip,low+ipl) = sig
				
				B(ip) = acon * glu * popl(up+ipl) / (popl(low+ipl) - glu * popl(up+ipl))
																	
         	if (ip == 1) then
					tau(it,ip) = chil(ip) * dz(ip)
					dtaudn(ip,:) = dchildn(ip,:) * dz(ip)
         	else
					tau(it,ip) = tau(it,ip-1) + chil(ip) * dz(ip)
					dtaudn(ip,:) = dtaudn(ip-1,:) + dchildn(ip,:) * dz(ip)
         	endif
			enddo
			
			call rt1d_cep(nz,it)
			
! Fill up the blocks of the vector Lstar_total and Jbar_total for all the transitions from the value of
! Jbar and Lstar coming from rt1d_cep
			do ip = 1, nz
				ipt = nr*(ip-1)
				if (optically_thin == 1) then
					Lstar_total(it+ipt) = 0.d0
					Jbar_total(it+ipt) = 0.d0
					dJbar_totaldn(it+ipt,:) = 0.d0
					flux_total(it+ipt) = 0.d0
				else
					Lstar_total(it+ipt) = lstar(ip)
					Jbar_total(it+ipt) = Jbar(ip)
					dJbar_totaldn(it+ipt,:) = dJbardn(ip,:)
					flux_total(it+ipt) = flux(ip)
				endif
				
			enddo
			
		enddo
  
   end subroutine calcJbar_Lstar_cep
   
!---------------------------------------------------------
! This subroutine solves the radiative transfer equation using the 
! coupled escape probability scheme developed by Elitzur & Asensio Ramos
!---------------------------------------------------------
   subroutine rt1d_cep(np,it)
	integer :: np, it
   integer :: i, j
   real(kind=8) :: tau_local, escape, der_alpha, der_beta, error, delta_tau(4), alpha_tau(4), tau_aqui(nact,nz+1)
	real(kind=8) :: alphap_tau(4), sign(4)

   Jbar = 0.d0
	Lstar = 0.d0
	flux = 0.d0
	dJbardn = 0.d0
		
	do i = 1, np
		tau_local = tau(it,i) - tau(it,i-1)
	
		escape = beta2(tau_local)
		
		Jbar(i) = Sl(i) * (1.d0 - escape)
		
		der_alpha = alphap2(tau_local)
		der_beta = (der_alpha - escape) / tau_local
		
		dJbardn(i,:) = (1.d0 - escape) * dSldn(i,:)
		dJbardn(i,:) = dJbardn(i,:) - Sl(i) * der_beta * (dtaudn(i,:)-dtaudn(i-1,:))
						
		Lstar(i) = (1.d0 - escape)
		do j = 1, np
			if (j /= i) then
				alpha_tau = 0.d0

! Alpha function and its derivative				
				delta_tau(1) = dabs(tau(it,i) - tau(it,j))
				escape = beta2(delta_tau(1))
				alpha_tau(1) = delta_tau(1) * escape
				if (tau(it,i) > tau(it,j)) then
					alphap_tau(1) = alphap2(delta_tau(1))
				else
					alphap_tau(1) = -alphap2(delta_tau(1))
				endif
												
				delta_tau(2) = dabs(tau(it,i-1) - tau(it,j))
				escape = beta2(delta_tau(2))				
				alpha_tau(2) = delta_tau(2) * escape 
				if (tau(it,i-1) > tau(it,j)) then
					alphap_tau(2) = alphap2(delta_tau(2))
				else
					alphap_tau(2) = -alphap2(delta_tau(2))
				endif
												
				delta_tau(3) = dabs(tau(it,i) - tau(it,j-1))				
				escape = beta2(delta_tau(3))
				alpha_tau(3) = delta_tau(3) * escape 
				if (tau(it,i) > tau(it,j-1)) then
					alphap_tau(3) = alphap2(delta_tau(3))
				else
					alphap_tau(3) = -alphap2(delta_tau(3))
				endif
								
				delta_tau(4) = dabs(tau(it,i-1) - tau(it,j-1))
				escape = beta2(delta_tau(4))
				alpha_tau(4) = delta_tau(4) * escape 
				if (tau(it,i-1) > tau(it,j-1)) then
					alphap_tau(4) = alphap2(delta_tau(4))
				else
					alphap_tau(4) = -alphap2(delta_tau(4))
				endif

! Mean intensity and its derivative with respect to the level populations
				Jbar(i) = Jbar(i) + 0.5d0 * Sl(j) / tau_local * &
					(alpha_tau(1) - alpha_tau(2) - alpha_tau(3) + alpha_tau(4))
					
! We separate the derivative in three terms
! Derivative of the source function
				dJbardn(i,:) = dJbardn(i,:) + 0.5d0 * dSldn(j,:) / tau_local * &
					(alpha_tau(1) - alpha_tau(2) - alpha_tau(3) + alpha_tau(4))

! Derivative of the alpha functions
				dJbardn(i,:) = dJbardn(i,:) + 0.5d0 * Sl(j) / tau_local * &
					(alphap_tau(1) * (dtaudn(i,:)-dtaudn(j,:)) - &
					 alphap_tau(2) * (dtaudn(i-1,:)-dtaudn(j,:)) - &
					 alphap_tau(3) * (dtaudn(i,:)-dtaudn(j-1,:)) + &
					 alphap_tau(4) * (dtaudn(i-1,:)-dtaudn(j-1,:)))
					 
! Derivative of the 1/tau
				dJbardn(i,:) = dJbardn(i,:) - 0.5d0 * Sl(j) / tau_local**2 * (dtaudn(i,:)-dtaudn(i-1,:)) * &
					(alpha_tau(1) - alpha_tau(2) - alpha_tau(3) + alpha_tau(4))
									
			endif
		enddo

! Total flux of the line (following Eq. (45) in the notes)
		flux(np) = flux(np) + (1.d0-Jbar(i)/Sl(i)) * Sl(i) * tau_local
	enddo
	
   end subroutine rt1d_cep
   
end module functions

! *********************************************************
! *********************************************************
! STATISTICAL EQUILIBRIUM EQUATIONS
! *********************************************************
! *********************************************************

module escape_probability
use functions
use constants
implicit none
	real(kind=8) :: A21, C21, T, nu0, g1, g2, total_tau, x(2)
contains

!-------------------------------------------------------------------
! Returns the value of the equations of the nonlinear set to be solved
!-------------------------------------------------------------------
	subroutine funcv(x,n,fvec)
	integer :: n
	real(kind=8) :: x(n), fvec(n)
	integer :: ip
	integer :: low, up, ipl, ipt, i, j, il, it
	real(kind=8) :: A(nl, nl), B(nl), populat(nl), producto(nl), cul, acon, blu, glu, djbar, poptot

		call calcJbar_Lstar_cep(x)
		
		do ip = 1, nz

! Offsets for storing
			ipl = nl*(ip-1)
			ipt = nr*(ip-1)
		
			A = 0.d0
		
! Upward collisional rates A(low,up)
			call collis(ip,A)		

! Calculate the downward rates and interchange the positions
			do i = 1, nl
				do j = 1, i-1
					cul = A(i,j)
! Rate into de j-th level				
					A(j,i) = cul
! Rate into de i-th level				
					A(i,j) = popl(i+ipl) / popl(j+ipl) * cul
				enddo
			enddo
		
! Active radiative rates (preconditioning)
			do it = 1, nact
		
				up = itran(1,it)
				low = itran(2,it)
				acon = TWOHC2 * dtran(2,it)**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				djbar = Jbar_total(it+ipt)
			
				A(low,up) = A(low,up) + glu * blu * (acon + djbar )
				A(up,low) = A(up,low) + blu * djbar
			
			enddo		
		
! Background radiative rates
			do it = nact + 1, nr
		
				up = itran(1,it)
				low = itran(2,it)
				acon = TWOHC2 * dtran(2,it)**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				A(low,up) = A(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
				A(up,low) = A(up,low) + blu * Jbar_total(it+ipt)
			
			enddo		
		
! The system is homogeneous with null determinant. We have to substitute one of the equations by a closure
! relation. The sum over one column gives the total output rate from level i
			do i = 1, nl

				do j = 1, i-1
					A(i,i) = A(i,i) - A(j,i)
				enddo
				do j = i+1, nl
					A(i,i) = A(i,i) - A(j,i)
				enddo
				
			enddo		
		
! Conservation equation
			B = 0.d0
			A(nl, 1:nl) = 1.d0
		
			poptot = abund(ip) * nh(ip)
			B(nl) = poptot
			
			do il = 1, nl
				populat(il) = pop(il+ipl)
			enddo
			
! Return the rate equations			
			producto = matmul(A,populat)-B
						
			do il = 1, nl				
				fvec(il+ipl) = producto(il)
			enddo		
		enddo			
		
	end subroutine funcv
	
!-------------------------------------------------------------------
! Returns the value of the equations of the nonlinear set to be solved
!-------------------------------------------------------------------
	subroutine funcv_analytic(x,n,fvec,fjac)
	integer :: n
	real(kind=8) :: x(n), fvec(n), fjac(n,n)
	integer :: ip
	integer :: low, up, ipl, ipt, i, j, il, it, il2, iz, ind
	real(kind=8) :: A(nl, nl), A2(nl,nl), B(nl), populat(nl), producto(nl), cul, acon, blu, glu, djbar, poptot
	real(kind=8), allocatable :: jac(:,:,:)

		allocate(jac(nl,nl,n))
		
		call calcJbar_Lstar_cep(x)
				
		do ip = 1, nz

! Offsets for storing
			ipl = nl*(ip-1)
			ipt = nr*(ip-1)
		
			A = 0.d0
			jac = 0.d0
		
! Upward collisional rates A(low,up)
			call collis(ip,A)		

! Calculate the downward rates and interchange the positions
			do i = 1, nl
				do j = 1, i-1
					cul = A(i,j)
! Rate into de j-th level				
					A(j,i) = cul
! Rate into de i-th level				
					A(i,j) = popl(i+ipl) / popl(j+ipl) * cul
				enddo
			enddo
					
! Active radiative rates (preconditioning)
			do it = 1, nact
		
				up = itran(1,it)
				low = itran(2,it)
				acon = TWOHC2 * dtran(2,it)**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				djbar = Jbar_total(it+ipt)
			
				A(low,up) = A(low,up) + glu * blu * (acon + djbar )				
				A(up,low) = A(up,low) + blu * djbar
				
				jac(low,up,:) = jac(low,up,:) + glu * blu * dJbar_totaldn(it+ipt,:)
				jac(up,low,:) = jac(up,low,:) + blu * dJbar_totaldn(it+ipt,:)
			
			enddo		
		
! Background radiative rates
			do it = nact + 1, nr
		
				up = itran(1,it)
				low = itran(2,it)
				acon = TWOHC2 * dtran(2,it)**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				A(low,up) = A(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
				A(up,low) = A(up,low) + blu * Jbar_total(it+ipt)
				
				jac(low,up,:) = jac(low,up,:) + glu * blu * dJbar_totaldn(it+ipt,:)
				jac(up,low,:) = jac(up,low,:) + blu * dJbar_totaldn(it+ipt,:)
			
			enddo		
		
! The system is homogeneous with null determinant. We have to substitute one of the equations by a closure
! relation. The sum over one column gives the total output rate from level i
			do i = 1, nl

				do j = 1, i-1
					A(i,i) = A(i,i) - A(j,i)
					jac(i,i,:) = jac(i,i,:) - jac(j,i,:)
				enddo
				do j = i+1, nl
					A(i,i) = A(i,i) - A(j,i)
					jac(i,i,:) = jac(i,i,:) - jac(j,i,:)
				enddo
				
			enddo		
			
			A2 = A
					
! Conservation equation
			B = 0.d0
			A(nl, 1:nl) = 1.d0
			jac(nl, 1:nl, :) = 0.d0
		
			poptot = abund(ip) * nh(ip)
			B(nl) = poptot
			
			do il = 1, nl
				populat(il) = x(il+ipl)
			enddo
			
! Return the rate equations			
			producto = matmul(A,populat)-B
						
			do il = 1, nl				
				fvec(il+ipl) = producto(il)								
			enddo

! Analytical Jacobian. Since the equations are F_j = sum(A_jk * n_k, k) - b_j, the Jacobian J_ji is 
! J_ji = dF_j/dn_i = sum(dA_jk/dn_i * n_k,k) + sum(A_jk * dn_k/dn_i,k) = sum(dA_jk/dn_i * n_k,k) + A_ji
			do il = 1, nl
				do il2 = 1, nl
					fjac(il+ipl,il2+ipl) = A(il,il2)
				enddo
				do il2 = 1, nl*nz					
					fjac(il+ipl,il2) = fjac(il+ipl,il2) + sum(jac(il,:,il2)*populat)
				enddo
			enddo
						
		enddo
				
		deallocate(jac)
		
	end subroutine funcv_analytic	
	
!-------------------------------------------------------------------
! Calculates the Jacobian using forward-differences
!-------------------------------------------------------------------            
	subroutine fdjac(n,x,fvec,df)
	integer :: n,np
	real(kind=8) :: df(n,n),fvec(n),x(n),EPS
	PARAMETER (EPS=1.d-6)
	integer :: i,j
	real(kind=8) :: h,temp,f(n)
	
		do j=1,n
			temp=x(j)
			h=EPS*dabs(temp)
			if(h.eq.0.d0)h=EPS
			x(j)=temp+h
			h=x(j)-temp
			call funcv(x,n,f)
			x(j)=temp
			do i=1,n 
				df(i,j)=(f(i)-fvec(i))/h
			enddo 
		enddo 
        
	end subroutine fdjac

!-------------------------------------------------------------------
! Returns the equations and the Jacobian of the nonlinear set to be solved at point x
!-------------------------------------------------------------------
	subroutine usrfun(x,n,fvec,fjac,derivatives)
	integer :: n, i, j
	logical :: derivatives
	real(kind=8) :: x(n), fvec(n), fjac(:,:)
	real(kind=8), allocatable :: fjac_num(:,:)
	
		allocate(fjac_num(n,n))
		
		fjac = 0.d0
! Analytical derivatives	
		if (derivatives) then
			call funcv_analytic(x,n,fvec,fjac)			
		else
! Numerical derivatives
			call funcv(x,n,fvec)
			call fdjac(n,x,fvec,fjac)
		endif
		
		deallocate(fjac_num)
	end subroutine usrfun
	
!-------------------------------------------------------------------
! Solves a system of nonlinear equations using the Newton method
!-------------------------------------------------------------------
	subroutine mnewt(ntrial,x,n,tolx,tolf,derivatives)
	integer :: n,ntrial
	real(kind=8) :: tolf,tolx,x(n)
	integer :: i,k
	logical :: derivatives
	logical :: calculate_jacobian
	real(kind=8) :: d,errf,errx,fvec(n),p(n), tol, amax
	real(kind=8), allocatable :: fjac(:,:)
	integer, allocatable :: indx(:)
	
		allocate(fjac(n,n))
		errx = 1.d0
		do k=1,ntrial
					
			call usrfun(x,n,fvec,fjac,derivatives)   !User subroutine supplies function values at x in fvec
			
! Scale each line of the Jacobian by its maximum to decrease the dynamical range
			do i = 1, n
				amax = maxval(dabs(fjac(:,i)))
				if (amax < 1.d-20) amax = 1.d-20
				fjac(:,i) = fjac(:,i) / amax
			enddo
                
			do i=1,n  !Check function convergence.
				errf=errf+dabs(fvec(i))
			enddo 
			if (errf <= tolf) then
				print *, 'Function absolute value is smaller than tolf : ', errf
				!return
			endif
			p = -fvec
			
			allocate(indx(n))
			tol = 1.d-10
			call ludcmp(fjac,indx,tol)
			call lubksb(fjac,indx,p)
			deallocate(indx)
		                			                
			errx=0.d0  ! Check root convergence.
			x = x + p
                
			!do i=1,n   !Update solution.
			!	errx=errx+dabs(p(i))                    
			!enddo
			
			errx = maxval(abs(p/x))
			
			if(errx <= tolx) then				
				print *, 'Iteration ',k, '   rel_error: ',errx
				return
			endif
			
			print *, 'Iteration ',k, '   rel_error: ',errx
		enddo 		
		
		deallocate(fjac)
        
	end subroutine mnewt

end module escape_probability


! *********************************************************
! *********************************************************
! *********************************************************
! MAIN PROGRAM
! *********************************************************
! *********************************************************
! *********************************************************

program jacobi
use constants
use variables
use functions
use escape_probability

   
! Initialization
   call init

! Precalculation of the beta and alpha functions
   call precalculate_beta
			
! If we start from optically thin, the system is linear. Solve it once	
	if (start_mode == 1 .and. read_previous_flag /= 1) then
		print *, 'Calculating optically thin solution...'
		optically_thin = 1
		call calcJbar_Lstar_cep(pop)

      call correct_populations
		optically_thin = 0
		print *, 'Done...'
	endif

! Use the Newton method with analytical derivatives
	if (which_scheme == 2) then		
		call mnewt(40,pop,nz*nl,precision,0.d0,.TRUE.)
	endif	

! Use the Newton method with numerical derivatives
	if (which_scheme == 3) then		
		call mnewt(40,pop,nz*nl,precision,0.d0,.FALSE.)	
	endif		

! Use the CEP ALI method
	if (which_scheme == 1) then
		relat_error = 1.d0
		iter = 1

		!precision = 1.d-8
		print *, 'Starting main loop...'
! Iterate until reaching the precision we want or we go above the maximum number of iterations
		do while (relat_error > precision .and. iter < n_iters)

! Calculate Jbar_total and Lstar_total solving the RT equation		
			call calcJbar_Lstar_cep(pop)		

! Do the population correction
      	call correct_populations

			iter = iter + 1

			print *, 'Iteration ',iter-1, '   rel_error: ',relat_error

		enddo
	endif
			
	call write_results
	
end program jacobi
