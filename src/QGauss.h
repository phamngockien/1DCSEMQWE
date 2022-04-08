//class for 1D quadrature rule (QGAUSS)
//*************************************************************************
/*This class is written based on Deal.II library
Bangerth, W., Hartmann, R., & Kanschat, G., 2007. deal. II—a general-purpose object-oriented finite element library. ACM Transactions on Mathematical Software (TOMS), 33(4), 24, doi: 10.1145/1268776.1268779.
You can use the existing class when coupling the code in Deal.II*/
//*************************************************************************
	
#ifndef QGAUSS_H
#define QGAUSS_H

#include <limits>
#include <cmath>
#include <exception>
#include <iostream>
#include <vector>

class QGauss
{
	 public:
	 QGauss(const unsigned int &nQuad);
	 
	 //in these functions quadrature points and weights are returned by their addresses without copying
	 const std::vector<long double> &get_points() const;
     const std::vector<double> &get_weights() const;
	 
	 private:
	 std::vector<long double> quadrature_points;
	 std::vector<double> weights;
	 
	 static const unsigned int invalid_unsigned_int = static_cast<unsigned int>(-1);
	 
	 template <typename Number>
     Number
     jacobi_polynomial_value(const unsigned int degree,
											  const int          alpha,
											  const int          beta,
											  const Number       x);
 
   template <typename Number>
   std::vector<Number>
   jacobi_polynomial_roots(const unsigned int degree,
											const int          alpha,
											const int          beta);
	
};

 inline const std::vector<long double> &QGauss::get_points() const
 {
   return quadrature_points;
 }
 
 inline const std::vector<double> &QGauss::get_weights() const
 {
   return weights;
 }
#endif