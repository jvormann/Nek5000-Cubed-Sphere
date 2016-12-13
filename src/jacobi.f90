!###################################
!##Subroutines and functions########
!##for Jacobi Polynomials###########
!##Author: Jan Vormann#############
!###################################
MODULE Jacobi_Subroutines
IMPLICIT NONE
CONTAINS
!
!#####################################################
!##Functions to calculate Jacobi Polynomials##########
!#####################################################
!
FUNCTION jacobi_pol(x,m,a,b)
!Calculate the value of the Jacobi polynomials up to degree m+1 at x
!Using recursion formula (Karniadakis, Sherwin(2013), p 586)
REAL(KIND=8), DIMENSION(0:m) :: jacobi_pol 
REAL(KIND=8) :: x,a,b
REAL(KIND=8), DIMENSION(4) :: ac
INTEGER :: m,j
jacobi_pol(0)=1.0_8
jacobi_pol(1)=0.5_8*(a-b+(a+b+2.0_8)*x)
DO j=2,m
  ac=jacobi_rec_const(j-1,a,b)
  jacobi_pol(j)=(1.0_8/ac(1))*((ac(2)+ac(3)*x)*jacobi_pol(j-1)-ac(4)*jacobi_pol(j-2))
ENDDO 
RETURN
END FUNCTION

!Functions to evaluate the constants for Jacobi recursion
function jacobi_rec_const(n,a,b)
REAL (KIND=8), DIMENSION(4) :: jacobi_rec_const
INTEGER :: n
REAL(KIND=8) :: a,b,nd
nd=DBLE(n)
jacobi_rec_const(1)=2.0_8*(nd+1.0_8)*(nd+a+b+1.0_8)*(2.0_8*nd+a+b)
jacobi_rec_const(2)=(2.0_8*nd+a+b+1.0_8)*(a**2.0_8-b**2.0_8)
jacobi_rec_const(3)=(2.0_8*nd+a+b)*(2.0_8*nd+a+b+1.0_8)*(2.0_8*n+a+b+2.0_8)
jacobi_rec_const(4)=2.0_8*(nd+a)*(nd+b)*(2.0_8*nd+a+b+2.0_8)
RETURN
END FUNCTION
!
!#######################################################
!##Functions for the derivatives of Jacobi polynomials##
!#######################################################
!
FUNCTION ddx_jacobi_pol(x,m,a,b)
!Calculate the value of the derivative of the Jacobi polynomials up to degree m at x
!Using recursion formula (Karniadakis, Sherwin(2013), p 586)
REAL(KIND=8), DIMENSION(0:m) :: ddx_jacobi_pol, jac_pol
REAL(KIND=8) :: x,a,b
REAL(KIND=8), DIMENSION(3) :: bc
INTEGER :: m,j
ddx_jacobi_pol(0)=0.0_8
ddx_jacobi_pol(1)=0.5_8*(a+b+2.0_8)
jac_pol=jacobi_pol(x,m,a,b)
DO j=2,m
  bc=ddx_jacobi_rec_const(j,a,b,x)
  ddx_jacobi_pol(j)=(bc(2)/bc(1))*jac_pol(j)+(bc(3)/bc(1))*jac_pol(j-1)
ENDDO
RETURN
END FUNCTION
!
function ddx_jacobi_rec_const(n,a,b,x)
!Function to evaluate the constants for Jacobi recursion (derivatives)
REAL (KIND=8), DIMENSION(3) :: ddx_jacobi_rec_const
INTEGER :: n
REAL(KIND=8) :: a,b,nd,x
nd=DBLE(n)
ddx_jacobi_rec_const(1)=(2.0_8*nd+a+b)*(1.0_8-x**2.0_8)
ddx_jacobi_rec_const(2)=nd*(a-b-(2.0_8*nd+a+b)*x)
ddx_jacobi_rec_const(3)=2.0_8*(nd+a)*(nd+b)
RETURN
END FUNCTION
!
!#################################################
!##Root finding for Jacobi polynomials############
!#################################################
!
function jacobi_roots(m,eps,a,b)
!Find the roots of a Jacobi Polynomial, using Newton-Raphson-Iteration and polynomial deflation
!See (Karniadakis, Sherwin(2013), p 598)
!Polynomal deflation ensures that we find a new root at each iteration
REAL (KIND=8), DIMENSION(0:m-1) :: jacobi_roots, xr
REAL(KIND=8), DIMENSION(0:m) :: djp,jp
REAL (KIND=8) :: r,pi,delta,eps,sum1,a,b
INTEGER :: k,isum,im,icount,m
pi=4.0_8*atan(1.0_8)
DO k=0,m-1 !Loop over all existing roots
 r=-cos((2.0_8*REAL(k)+1.0_8)*pi/(2.0_8*m)) !The Chebyshev roots have an explicit form, so we use them as initial values
 IF (k .GT. 0) r=(r+jacobi_roots(k-1))/2.0_8 !A better approximation is given by the average of r and the last root found
 delta=eps+1.0_8
 DO WHILE (delta .GT. eps) !Iterate until the root is found to a given tolerance eps
  sum1=0.0_8
  icount=0
  DO isum=0,k-1
    sum1=sum1+(1.0_8/(r-jacobi_roots(isum)))
  ENDDO
  jp=jacobi_pol(r,m,a,b)
  djp=ddx_jacobi_pol(r,m,a,b)
  delta=-jp(m)/(djp(m)-jp(m)*sum1)
  r=r+delta
  icount=icount+1
  IF (icount .GT. 10000) THEN
    PRINT*,'Iteration stopped after 10000 steps.'
    EXIT
  ENDIF
  ENDDO
  jacobi_roots(k)=r
ENDDO
RETURN
END FUNCTION
!
!#######################################################################
!###Functions for Gauss-Lobatto-Legendre-Integration and Differentaion##
!#######################################################################exit
!
FUNCTION gll_quad_points(m)
!Returns the Gauss-Lobatto-Legendre Quadrature Points for order m
!Karniadakis, Sherwin (2013), p 61
REAL (KIND=8), DIMENSION(0:m-1) :: gll_quad_points
INTEGER :: m
gll_quad_points(0)=-1.0_8
gll_quad_points(1:m-2)=jacobi_roots(m-2,1E-8_8,1.0_8,1.0_8)
gll_quad_points(m-1)=1.0_8
RETURN
END FUNCTION
!
FUNCTION gll_weights(m)
!Returns the weight for the Gauss-Lobatto-Legendre-Integration of order m
!Karniadakis, Sherwin (2013), p 61
REAL (KIND=8), DIMENSION(0:m-1) :: gll_weights, qp, jp
INTEGER :: m,im
qp=gll_quad_points(m)
gll_weights(0)=2.0_8/DBLE(m*(m-1))
gll_weights(m-1)=2.0_8/DBLE(m*(m-1))
DO im=1,m-2
  jp(0:m-1)=jacobi_pol(qp(im),m-1,0.0_8,0.0_8)
  gll_weights(im)=2.0_8/(DBLE(m*(m-1))*(jp(m-1)**2.0_8))
ENDDO
RETURN
END FUNCTION
!
FUNCTION gll_diff_matrix(m)
!Returns the Differentation Matrix for use with Gauss-Lobatto-Legendre-Points
!Karniadakis, Sherwin (2013), p 68
REAL (KIND=8), DIMENSION (0:m-1, 0:m-1) :: gll_diff_matrix
REAL (KIND=8), DIMENSION(0:m-1) :: qp, jp_im,jp_jm
INTEGER :: m, im, jm
qp=gll_quad_points(m)
DO im=0,m-1
 DO jm=0,m-1
  IF ((im .NE. jm) .AND. (im .GE. 0) .AND. (jm .LE. m-1)) THEN
    jp_im(0:m-1)=jacobi_pol(qp(im),m-1,0.0_8,0.0_8)
    jp_jm(0:m-1)=jacobi_pol(qp(jm),m-1,0.0_8,0.0_8)
    gll_diff_matrix(im,jm)=jp_im(m-1)/(jp_jm(m-1)*(qp(im)-qp(jm)))
  ELSE IF ((im .GE. 1) .AND. (im==jm) .AND. (jm .LE. m-2)) THEN
    gll_diff_matrix(im,jm)=0
  ELSE IF ((im==0) .AND. (jm==0)) THEN
    gll_diff_matrix(im,jm)=-DBLE(m*(m-1))/4.0_8
  ELSE IF ((im==m-1) .AND. (jm==m-1)) THEN
    gll_diff_matrix(im,jm)=DBLE(m*(m-1))/4.0_8
  ENDIF
 ENDDO
ENDDO
RETURN
END FUNCTION
!
FUNCTION gll_integral(f,l,u,m)
!Returns the Integral over the function f from l to u, using m Quadraturepoints
!Uses Gauss Lobatto Legendre Quadrature
REAL (KIND=8) :: gll_integral,l,u,fak
REAL (KIND=8), DIMENSION(0:m-1) :: gll_w, gll_qp
REAL (KIND=8) :: f
INTEGER :: m, im
!fak maps the Intervall [l:u] to [-1:1]
fak=(u-l)/2.0_8
gll_w=gll_weights(m)
gll_qp=gll_quad_points(m)
gll_integral=0.0_8
DO im=0,m-1
 gll_integral=gll_integral+gll_w(im)*f(fak*gll_qp(im)+(l+u)/2.0_8)
ENDDO
gll_integral=fak*gll_integral
RETURN
END FUNCTION
!
FUNCTION gll_derivative(f,l,u,m)
!Returns the derivative of the function f at the m GLL Collocation points
!The Points are mapped to the intervall [l,u]
!Uses Gauss Lobatto Legendre Quadrature
REAL (KIND=8), DIMENSION (0:m-1) :: gll_derivative
REAL (KIND=8), DIMENSION (0:m-1) :: gll_qp, f_val !Quadraturepoints, function values at qp
REAL (KIND=8), DIMENSION (0:m-1,0:m-1) :: gll_dm !Differentiationmatrix
REAL (KIND=8) :: f, fak, l, u
INTEGER :: m,im,jm
fak=(u-l)/2.0_8
gll_qp=gll_quad_points(m)
gll_dm=gll_diff_matrix(m)
DO jm=0,m-1
 f_val(jm)=f(fak*gll_qp(jm)+(l+u)/2.0_8)
ENDDO
gll_derivative=0.0_8
DO im=0,m-1
  DO jm=0,m-1
    gll_derivative(im)=gll_derivative(im)+gll_dm(im,jm)*f_val(jm)
  ENDDO
ENDDO
!Multiply with Jacobian dxsi/dx=2/(u-l)
!Chain rule: du/dx=(du/dxsi)(dxsi/dx)
gll_derivative=(2.0_8/(u-l))*gll_derivative
RETURN
END FUNCTION
END MODULE Jacobi_Subroutines
