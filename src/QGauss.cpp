//class for 1D quadrature rule (QGAUSS)

#include "QGauss.h"

 QGauss:: QGauss(const unsigned int &nQuad)
   : quadrature_points(nQuad)
   , weights(nQuad)
 {
	 if (nQuad == 0) return;
 
    std::vector<long double> points = jacobi_polynomial_roots<long double>(nQuad, 0, 0);
 
   for (unsigned int i = 0; i < (points.size() + 1) / 2; ++i)
     {
       this->quadrature_points[i]        = points[i];
       this->quadrature_points[nQuad - i - 1]= 1. - points[i];
 
       // derivative of Jacobi polynomial
       const long double pp = 0.5 * (nQuad + 1) * jacobi_polynomial_value(nQuad - 1, 1, 1, points[i]);
       const long double x      = -1. + 2. * points[i];
       const double      w      = 1. / ((1. - x * x) * pp * pp);
       this->weights[i]         = w;
       this->weights[nQuad - i - 1] = w;
     }
 }
 
 
 template <typename Number>
   Number
   QGauss::jacobi_polynomial_value(const unsigned int degree,
                           const int          alpha,
                           const int          beta,
                           const Number       x)
   {
    if (alpha < 0 || beta < 0){
		
		std::cout << "Negative alpha/beta coefficients not supported \n";
		throw std::exception();	
	} 
	
     // the Jacobi polynomial is evaluated using a recursion formula.
     Number p0, p1;
 
     // The recursion formula is defined for the interval [-1, 1], so rescale
     // to that interval here
     const Number xeval = Number(-1) + 2. * x;
 
     // initial values P_0(x), P_1(x):
     p0 = 1.0;
     if (degree == 0)
       return p0;
     p1 = ((alpha + beta + 2) * xeval + (alpha - beta)) / 2;
     if (degree == 1)
       return p1;
 
     for (unsigned int i = 1; i < degree; ++i)
       {
         const Number v  = 2 * i + (alpha + beta);
         const Number a1 = 2 * (i + 1) * (i + (alpha + beta + 1)) * v;
         const Number a2 = (v + 1) * (alpha * alpha - beta * beta);
         const Number a3 = v * (v + 1) * (v + 2);
         const Number a4 = 2 * (i + alpha) * (i + beta) * (v + 2);
 
         const Number pn = ((a2 + a3 * xeval) * p1 - a4 * p0) / a1;
         p0              = p1;
         p1              = pn;
       }
     return p1;
   }
 
 
 
   template <typename Number>
   std::vector<Number>
   QGauss::jacobi_polynomial_roots(const unsigned int degree,
                           const int          alpha,
                           const int          beta)
   {
     std::vector<Number> x(degree, 0.5);
 
     // compute zeros with a Newton algorithm.
 
     // Set tolerance. For long double we might not always get the additional
     // precision in a run time environment (e.g. with valgrind), so we must
     // limit the tolerance to double. Since we do a Newton iteration, doing
     // one more iteration after the residual has indicated convergence will be
     // enough for all number types due to the quadratic convergence of
     // Newton's method
 
     const Number tolerance =
       4 * std::max(static_cast<Number>(std::numeric_limits<double>::epsilon()),
                    std::numeric_limits<Number>::epsilon());
 
     // The following implementation follows closely the one given in the
     // appendix of the book by Karniadakis and Sherwin: Spectral/hp element
     // methods for computational fluid dynamics (Oxford University Press,
     // 2005)
 
     // If symmetric, we only need to compute the half of points
     const unsigned int n_points = (alpha == beta ? degree / 2 : degree);
     for (unsigned int k = 0; k < n_points; ++k)
       {
         // we take the zeros of the Chebyshev polynomial (alpha=beta=-0.5) as
         // initial values, corrected by the initial value
         Number r = 0.5 - 0.5 * std::cos(static_cast<Number>(2 * k + 1) /
                                         (2 * degree) * std::acos(-1)); //PI = acos(-1)
         if (k > 0)
           r = (r + x[k - 1]) / 2;
 
         unsigned int converged = invalid_unsigned_int;
         for (unsigned int it = 1; it < 1000; ++it)
           {
             Number s = 0.;
             for (unsigned int i = 0; i < k; ++i)
               s += 1. / (r - x[i]);
 
             // derivative of P_n^{alpha,beta}, rescaled to [0, 1]
             const Number J_x =
               (alpha + beta + degree + 1) *
               jacobi_polynomial_value(degree - 1, alpha + 1, beta + 1, r);
 
             // value of P_n^{alpha,beta}
             const Number f = jacobi_polynomial_value(degree, alpha, beta, r);
             const Number delta = f / (f * s - J_x);
             r += delta;
             if (converged == invalid_unsigned_int &&
                 std::abs(delta) < tolerance)
               converged = it;
 
             // do one more iteration to ensure accuracy also for tighter
             // types than double (e.g. long double)
             if (it == converged + 1)
               break;
           }
 
        if(converged == invalid_unsigned_int){
			std::cout << "Newton iteration for zero of Jacobi polynomial did not converge. \n";
			throw std::exception();
		} 
 
         x[k] = r;
       }
 
     // in case we assumed symmetry, fill up the missing values
     for (unsigned int k = n_points; k < degree; ++k)
       x[k] = 1.0 - x[degree - k - 1];
 
     return x;
   }