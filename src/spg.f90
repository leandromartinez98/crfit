subroutine callspg(n,x,f,outputunit,evalf,evalg,optpars)

  implicit none
  integer :: optpars(10)
  integer :: n, fcnt, inform, iprint, iter, maxfc, maxit, spginfo, outputunit
  double precision :: x(n), f, gpsupn, epsopt
  external :: evalf, evalg

  maxfc = optpars(1)

 ! epsopt  (i)   tolerance for the convergence criterion

  epsopt = 1.d-4

 ! maxit   (i)   maximum number of iterations

  maxit = maxfc

 ! iprint  (i)   controls output level (0 = no print)

  if ( outputunit == 0 ) then
    iprint = 0
  else if ( outputunit < 0 ) then
    iprint = 2
  else
    iprint = 1
  end if

  call spg(n,x,epsopt,maxit,maxfc,iprint,f,gpsupn,iter,fcnt,&
           spginfo,inform,outputunit,evalf,evalg) 
 
  inform = 0

  return
end subroutine callspg

! There are no explicit restrictions in this problem, so just
! return (there is nothing to project)

subroutine proj(n,x,flag)
  integer :: n, flag
  double precision :: x(n)
  flag = 0
  return
end subroutine

!     =================================================================
!     File: spg.f
!     =================================================================

!     =================================================================
!     Module: Spectral Projected Gradient Method
!     =================================================================

!     Last update of any of the component of this module:

!     December 20, 2007.

!     Users are encouraged to download periodically updated versions of
!     this code at the TANGO Project web page:
!
!     www.ime.usp.br/~egbirgin/tango/

!     *****************************************************************
!     *****************************************************************

      subroutine spg(n,x,epsopt,maxit,maxfc,iprint,f,gpsupn,iter,fcnt,&
      spginfo,inform,outputunit,evalf,evalg)

      implicit none

!     SCALAR ARGUMENTS
      double precision gpsupn,epsopt,f
      integer fcnt,inform,iprint,iter,maxfc,maxit,n,spginfo,outputunit

!     ARRAY ARGUMENTS
      double precision x(n)

!     Subroutine SPG implements the Spectral Projected Gradient Method 
!     (Version 2: "Feasible continuous projected path") to find a 
!     local minimizers of a given function with convex constraints, 
!     described in
!
!     E.G. Birgin, J.M. Martinez and M. Raydan, "Nonmonotone spectral
!     projected gradient methods for convex sets", SIAM Journal on
!     Optimization 10, pp. 1196-1211, 2000.
!
!     The user must supply the external subroutines evalf, evalg and 
!     proj to evaluate the objective function and its gradient and to 
!     project an arbitrary point onto the feasible region.
!
!     This version 20 DEC 2007 by E.G.Birgin, J.M.Martinez and M.Raydan.

!     Other parameters (i means input, o means output):
!
!     n       (i)   number of variables
!     x       (i/o) initial guess on input, solution on output
!     epsopt  (i)   tolerance for the convergence criterion
!     maxit   (i)   maximum number of iterations
!     maxfc   (i)   maximum number of functional evaluations
!     iprint  (i)   controls output level (0 = no print)
!     f       (o)   functional value at the solution
!     gpsupn  (o)   sup-norm of the projected gradient at the solution
!     iter    (o)   number of iterations
!     fcnt    (o)   number of functional evaluations
!     spginfo (o)   indicates the reason for stopping
!     inform  (o)   indicates an error in an user supplied subroutine

!     spginfo:
!
!     0: Small continuous-projected-gradient norm
!     1: Maximum number of iterations reached
!     2: Maximum number of functional evaluations reached
!
!     spginfo remains unset if inform is not equal to zero on output

!     inform:
!
!       0: ok
!     -90: error in the user supplied evalf subroutine
!     -91: error in the user supplied evalg subroutine
!     -92: error in the user supplied proj  subroutine

!     PARAMETERS
      integer m,nmax
      double precision lmax,lmin
      parameter ( m    =     100 )
      parameter ( nmax =  100000 )
      parameter ( lmin = 1.0d-30 )
      parameter ( lmax = 1.0d+30 )

!     LOCAL SCALARS
      integer i,lsinfo
      double precision fbest,fnew,lambda,sts,sty

!     LOCAL ARRAYS
      double precision g(nmax),gnew(nmax),gp(nmax),s(nmax),y(nmax),&
              d(nmax),xbest(nmax),xnew(nmax),lastfv(0:m-1)

!     EXTERNAL SUBROUTINES
      external ls,sevalf,sevalg,sproj,evalf,evalg

!     INTRINSIC FUNCTIONS
      intrinsic abs,max,min,mod

!     ==================================================================
!     Initialization
!     ==================================================================

!     Print problem information

      if ( iprint .eq. 1 ) then
          write(outputunit,fmt=1000)
          write(outputunit,fmt=1010) n
      else if ( iprint .eq. 2 ) then
          write(* ,fmt=1000)
          write(* ,fmt=1010) n
      end if

!     Set some initial values:

!     error tracker
      inform = 0

!     for counting number of iterations as well as functional evaluations
      iter = 0
      fcnt = 0

!     for the non-monotone line search
      do i = 0,m - 1
          lastfv(i) = - 1.0d+99
      end do

!     Project initial guess

      call sproj(n,x,inform)

!     Compute function and gradient at the initial point

      call sevalf(n,x,f,evalf)

      fcnt = fcnt + 1

      call sevalg(n,x,g,evalg)

!     Store functional value for the non-monotone line search

      lastfv(0) = f

!     Compute continuous-project-gradient and its sup-norm

      do i = 1,n
          gp(i) = x(i) - g(i)
      end do

      call sproj(n,gp,inform)
      if (inform .ne. 0) return

      gpsupn = 0.0d0
      do i = 1,n
          gp(i) = gp(i) - x(i)
          gpsupn = max( gpsupn, abs( gp(i) ) )
      end do

!     Initial steplength

      if ( gpsupn .ne. 0.0d0) then
          lambda =  min( lmax, max( lmin, 1.0d0 / gpsupn ) )
      else
          lambda = 0.0d0
      end if

!     Initiate best solution and functional value

      fbest = f

      do i = 1,n
          xbest(i) = x(i)
      end do

!     ==================================================================
!     Main loop
!     ==================================================================

 100  continue

!     Print iteration information

      if ( iprint .eq. 1 ) then
          if ( mod(iter,10) .eq. 0 ) then
              write(outputunit,fmt=1020)
          end if
          write(outputunit,fmt=1030) iter,f,gpsupn
      else if ( iprint .eq. 2 ) then
          if ( mod(iter,10) .eq. 0 ) then
              write(* ,fmt=1020)
          end if
          write(* ,fmt=1030) iter,f,gpsupn
      end if

!      open(20,file='spg-tabline.out')
!      write(20,fmt=1040) n,iter,fcnt,f,gpsupn
!      close(20)

!     ==================================================================
!     Test stopping criteria
!     ==================================================================

!     Test whether the continuous-projected-gradient sup-norm
!     is small enough to declare convergence

      if ( gpsupn .le. epsopt ) then
          spginfo = 0

          if ( iprint .eq. 1 ) then
              write(outputunit,1100)
          else if ( iprint .eq. 2 ) then
              write(*, 1100)
          end if

          go to 200
      end if

!     Test whether the number of iterations is exhausted

      if (iter .ge. maxit) then
          spginfo = 1

          if ( iprint .eq. 1 ) then
              write(outputunit,1110)
          else if ( iprint .eq. 2 ) then
              write(*, 1110)
          end if

          go to 200
      end if

!     Test whether the number of functional evaluations

      if (fcnt .ge. maxfc) then
          spginfo = 2

          if ( iprint .eq. 1 ) then
              write(outputunit,1120)
          else if ( iprint .eq. 2 ) then
              write(*, 1120)
          end if

          go to 200
      end if

!     ==================================================================
!     Iteration
!     ==================================================================

      iter = iter + 1

!     Compute search direction

      do i = 1,n
          d(i) = x(i) - lambda * g(i)
      end do

      call sproj(n,d,inform)
      if (inform .ne. 0) return

      do i = 1,n
          d(i) = d(i) - x(i)
      end do

!     Perform safeguarded quadratic interpolation along the spectral 
!     continuous projected gradient

      call ls(n,x,f,g,d,m,lastfv,maxfc,fcnt,fnew,xnew,lsinfo,inform,evalf)
      if ( inform .ne. 0 ) return

      if ( lsinfo .eq. 2 ) then
          spginfo = 2

          if ( iprint .eq. 1 ) then
              write(outputunit,1120)
          else if ( iprint .eq. 2 ) then
              write(*, 1120)
          end if

          go to 200
      end if

!     Set new functional value and save it for the non-monotone line 
!     search

      f = fnew
      lastfv(mod(iter,m)) = f

!     Gradient at the new iterate

      call sevalg(n,xnew,gnew,evalg)

!     Compute s = xnew - x and y = gnew - g, <s,s>, <s,y>, the 
!     continuous-projected-gradient and its sup-norm

      sts = 0.0d0
      sty = 0.0d0
      do i = 1,n
          s(i)  = xnew(i) - x(i)
          y(i)  = gnew(i) - g(i)
          sts   = sts + s(i) ** 2
          sty   = sty + s(i) * y(i)
          x(i)  = xnew(i)
          g(i)  = gnew(i)
          gp(i) = x(i) - g(i)
      end do

      call sproj(n,gp,inform)
      if ( inform .ne. 0 ) return

      gpsupn = 0.0d0
      do i = 1,n
          gp(i) = gp(i) - x(i)
          gpsupn = max( gpsupn, abs( gp(i) ) )
      end do

!     Spectral steplength

      if ( sty .le. 0.0d0 ) then
          lambda = lmax

      else
          lambda = max( lmin, min( sts / sty, lmax ) )
      end if

!     Best solution and functional value

      if ( f .lt. fbest ) then
          fbest = f

          do i = 1,n
              xbest(i) = x(i)
          end do
      end if

!     ==================================================================
!     Iterate
!     ==================================================================

      go to 100

!     ==================================================================
!     End of main loop
!     ==================================================================

 200  continue

!     ==================================================================
!     Write statistics
!     ==================================================================

      if ( iprint .eq. 1 ) then
          write(outputunit,fmt=2000) iter,fcnt,f,gpsupn
      else if ( iprint .eq. 2 ) then
          write(* ,fmt=2000) iter,fcnt,f,gpsupn
      end if

!     ==================================================================
!     Finish returning the best point
!     ==================================================================

      f = fbest

      do i = 1,n
          x(i) = xbest(i)
      end do

!     ==================================================================
!     NON-EXECUTABLE STATEMENTS
!     ==================================================================

 1000 format(/,1X,78('='),&
             /,1X,'This is the SPECTRAL PROJECTED GRADIENT (SPG) for ',&
                  'for convex-constrained',/,1X,'optimization. If you ',&
                  'use this code, please, cite:',/,&
             /,1X,'E. G. Birgin, J. M. Martinez and M. Raydan, ',&
                  'Nonmonotone spectral projected',/,1X,'gradient ',&
                  'methods on convex sets, SIAM Journal on ',&
                  'Optimization 10, pp.',/,1X,'1196-1211, 2000, and',/,&
             /,1X,'E. G. Birgin, J. M. Martinez and M. Raydan, ',&
                  'Algorithm 813: SPG - software',/,1X,'for ',&
                  'convex-constrained optimization, ACM Transactions ',&
                  'on Mathematical',/,1X,'Software 27, pp. 340-349, ',&
                  '2001.',/,1X,78('='))

 1010 format(/,1X,'Entry to SPG.',/,1X,'Number of variables: ',I7)

 1020 format(/,4X,'ITER',10X,'F',8X,'GPSUPN')

 1030 format(  1X,I7,1X,1P,D16.8,1X,1P,D7.1)
!1040 format(  1X,I7,1X,I7,1X,I7,1X,1P,D16.8,1X,1P,D7.1,1X,'(Abnormal ',&
!               'termination. Probably killed by CPU time limit.)')

 1100 format(/,1X,'Flag of SPG: Solution was found.')
 1110 format(/,1X,'Flag of SPG: Maximum of iterations reached.')
 1120 format(/,1X,'Flag of SPG: Maximum of functional evaluations reached.')
!1130 format(/,1X,'Flag of SPG: Too small step in the line search.',&
!             /,1X,'Probably, an exaggerated small norm of the ',&
!                  'continuous projected gradient',&
!             /,1X,'is being required for declaring convergence.')

 2000 format(/,1X,'Number of iterations               : ',9X,I7,&
             /,1X,'Number of functional evaluations   : ',9X,I7,&
             /,1X,'Objective function value           : ',1P,D16.8,&
             /,1X,'Sup-norm of the projected gradient : ',9X,1P,D7.1)

      end

!     *****************************************************************
!     *****************************************************************

      subroutine ls(n,x,f,g,d,m,lastfv,maxfc,fcnt,fnew,xnew,lsinfo,&
      inform,evalf)

      implicit none

!     SCALAR ARGUMENTS
      integer inform,lsinfo,maxfc,fcnt,m,n
      double precision f,fnew

!     ARRAY ARGUMENTS
      double precision d(n),g(n),lastfv(0:m-1),x(n),xnew(n)

!     Nonmonotone line search with safeguarded quadratic interpolation
 
!     lsinfo:
!
!     0: Armijo-like criterion satisfied
!     2: Maximum number of functional evaluations reached

!     PARAMETERS
      double precision gamma
      parameter ( gamma     = 1.0d-04 )

!     LOCAL SCALARS
      integer i
      double precision alpha,atemp,fmax,gtd

!     EXTERNAL SUBROUTINES
      external sevalf, evalf

!     INTRINSIC FUNCTIONS
      intrinsic abs,max

!     Initiate

      inform = 0
      fmax = lastfv(0)
      do i = 1,m - 1
          fmax = max( fmax, lastfv(i) )
      end do

      gtd = 0.0d0
      do i = 1,n
         gtd = gtd + g(i) * d(i)
      end do

      alpha = 1.0d0

      do i = 1,n
          xnew(i) = x(i) + alpha * d(i)
      end do

      call sevalf(n,xnew,fnew,evalf)

      fcnt = fcnt + 1

!     Main loop

 100  continue

!     Test stopping criteria

      if ( fnew .le. fmax + gamma * alpha * gtd ) then
          lsinfo = 0
          return
      end if

      if (fcnt .ge. maxfc) then
          lsinfo = 2
          return
      end if 

!     Safeguarded quadratic interpolation

      if ( alpha .le. 0.1d0 ) then
          alpha = alpha / 2.0d0

      else
          atemp = ( - gtd * alpha ** 2 ) / ( 2.0d0 * ( fnew - f - alpha * gtd ) )

          if ( atemp .lt. 0.1d0 .or. atemp .gt. 0.9d0 * alpha ) then
              atemp = alpha / 2.0d0
          end if

          alpha = atemp
      end if

!     New trial

      do i = 1,n
          xnew(i) = x(i) + alpha * d(i)
      end do

      call sevalf(n,xnew,fnew,evalf)

      fcnt = fcnt + 1

!     Iterate

      go to 100

      end

!     *****************************************************************
!     *****************************************************************

      subroutine sevalf(n,x,f,evalf)

      implicit none

!     SCALAR ARGUMENTS
      integer n 
!      integer inform,n
      double precision f

!     ARRAY ARGUMENTS
      double precision x(n)

!     LOCAL SCALARS
!      integer flag

!     EXTERNAL SUBROUTINES
      external evalf,reperr

      call evalf(x,f)

!     This is true if f if Inf, - Inf or NaN
      if ( .not. f .gt. - 1.0d+99 .or. .not. f .lt. 1.0d+99 ) then
          f = 1.0d+99
      end if

      end

!     *****************************************************************
!     *****************************************************************

      subroutine sevalg(n,x,g,evalg)

      implicit none

!     SCALAR ARGUMENTS
      integer n
!      integer inform,n

!     ARRAY ARGUMENTS
      double precision g(n),x(n)

!     LOCAL SCALARS
!      integer flag

!     EXTERNAL SUBROUTINES
      external evalg,reperr

      call evalg(n,x,g)

      end

!     *****************************************************************
!     *****************************************************************

      subroutine sproj(n,x,inform)

      implicit none

!     SCALAR ARGUMENTS
      integer inform,n

!     ARRAY ARGUMENTS
      double precision x(n)

!     LOCAL SCALARS
      integer flag

!     EXTERNAL SUBROUTINES
      external proj,reperr

      call proj(n,x,flag)

      if ( flag .ne. 0 ) then
          inform = - 92
          call reperr(inform)
          return
      end if

      end

!     ******************************************************************
!     ******************************************************************

      subroutine reperr(inform)

      implicit none

!     SCALAR ARGUMENTS
      integer inform

      if ( inform .eq. -90 ) then
          write(* ,fmt=100) 'EVALF'

      else if ( inform .eq. -91 ) then
          write(* ,fmt=100) 'EVALG'

      else if ( inform .eq. -92 ) then
          write(* ,fmt=100) 'PROJ '
      end if

!     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'*** There was an error in the user supplied ',&
                  'subroutine ',A10,' ***',/)

      end
