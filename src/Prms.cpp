//CLASS FOR STORING ALL NECESSARY PARAMETERS

#include "Prms.h"

PRMS:: PRMS(){}


//=======================================================================
//public functions
//The precomputed part for all receivers and all types of source
void PRMS::pre_compute()
{
	PI = std::acos(-1);
	mu0 = 4. * PI * 1e-7;
	
	//We need xInt (Interval global x) to compute the ouputs
	getBesselZeros();
    
	//pre-compute the Bessel function kernel and the quadrature points
	getBesselxW();

}


//=======================================================================
//private functions
//---------------------------------------------------------------------------------------------------
//function returns xInt   - breakpoints for dividing up the global integral
// where the Bessel function of order nu are zeros
//Thus x is the solution of J_nu(x) = 0.
//We solve it by Newton iteration method
void PRMS::getBesselZeros() {

    const int nu = 1; //compute the Zeros of J1, which is the maxima of J0

    std::vector<double> x(nIntervals);

    //initial guess using asymtotic zeros
    for (unsigned int j = 0; j < nIntervals; ++j) {
        x[j] = (j + 1) * PI + nu * PI / 2 - PI / 4;
    }

    //for evaluate values of Bessel functions
    std::vector<double> y0(nIntervals), y1(nIntervals); //y0 = J_nu(x); y1 = J_nu+1 (x)
    std::vector<double> step(nIntervals);

    //for checking convergence
    double temp_norm;

    //Newton iterations
    //usually stop at 5
    for (unsigned int k = 0; k < 10; ++k) {

        temp_norm = 0.;

        for (unsigned int j = 0; j < nIntervals; ++j) {
            //y0[j] = std::cyl_bessel_j(nu, x[j]);
            //y1[j] = std::cyl_bessel_j(nu + 1, x[j]);
			y0[j] =  boost::math::cyl_bessel_j(nu, x[j]);
            y1[j] =  boost::math::cyl_bessel_j(nu + 1, x[j]);

            //the step size 
            step[j] = (-1.) * y0[j] / (nu * y0[j] / x[j] - y1[j]); //=-J_nu/J'_nu

            //update the new solution to x
            x[j] += step[j];

            //check convergence
            temp_norm += std::pow(step[j], 2);
        }

        //if the solution converged then break the look
        temp_norm = sqrt(temp_norm);
        if (temp_norm < 1e-12) {
            break;
        }
    }

    //++++++evaluating breakpoints++++++++
    xInt.resize(nIntervals);
    //get points at which J1(x) = 0;
    xInt = x;
    //add x = 0 (1e-20) at the begin of xInt
    //note: now xInt.size() = nIntervals + 1
    xInt.insert(xInt.begin(), 1e-20);

    //++++++++++++++++++++++++++++++++++++++++++++++++++++
}
//---- end getBesselZeros
//---------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
//function returns: quadrature intervals (real coords) and Bessel*weight at each quadrature point
//Output:
//Bx     - global vector of all quadrature points between all breakpoints: Bx = lambda * rho
//J0xW   - vector of weight(Bx) * J0(Bx) with J0 is Bessel function of the first kind of order 0
//J1xW   - vector of weight(Bx) * J1(Bx) with J1 is Bessel function of the first kind of order 1
// where rho is the horizontal distance between Tx (transmitter) and Rx_list (receiver)
// lambda = sqrt(k_x^2 + k_y^2) with k_x and k_y are the spatial wave number in x- and y-directions
void PRMS::getBesselxW() {

    //reserve memory for Bx, J0xW, J1xW
    Bx.resize(nIntervals * nQuad);
    J0xW.resize(nIntervals * nQuad);
    J1xW.resize(nIntervals * nQuad);

    //utilize Gauss quadrature rule from Deal.II lib, which defined in unit [0,1]
    QGauss quadrature_formula(nQuad);

    //get quadrature weights
    std::vector<double> weights;
    weights = quadrature_formula.get_weights();

    //quadrature points in unit coords for each interval
    std::vector<long double> x_hat;
    x_hat = quadrature_formula.get_points();

    //temporary global x
    double temp_x;

    for (unsigned int l = 0; l < nIntervals; ++l) {

        for (unsigned int k = 0; k < nQuad; ++k) {
            temp_x = (xInt[l + 1] - xInt[l]) * x_hat[k] + xInt[l];
            Bx[k + l * nQuad] = temp_x;
            //J0xW[k + l * nQuad] = std::cyl_bessel_j(0, temp_x) * weights[k];
            //J1xW[k + l * nQuad] = std::cyl_bessel_j(1, temp_x) * weights[k];
			J0xW[k + l * nQuad] = boost::math::cyl_bessel_j(0, temp_x) * weights[k];
            J1xW[k + l * nQuad] = boost::math::cyl_bessel_j(1, temp_x) * weights[k];
			
        }
    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}
//-- end getBesselxW
//------------------------------------------------------------------------------------------------
