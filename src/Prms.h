//CLASS FOR STORING ALL NECESSARY PARAMETERS

#ifndef PRMS_H
#define PRMS_H

#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "QGauss.h"

 #include <boost/math/special_functions/bessel.hpp>  //for bessel function

class PRMS
{
	public:
	
	PRMS();
	
	unsigned int dipole_type; //type of dipole 0 for electric, 1 for magnetic
	
	unsigned int nQuad;  //number of quadrature points
	
	//absolute and relative tolerances for QWE convergence criteria
	double rtol;
	double atol;
	
	unsigned int nfqs; //number of frequencies
	std::vector<double> frequencies; //their values
	
	unsigned int nLayers; //Number of layers inside the model
	std::vector<double> sig; //arrays of conductivities of the layers (S/m)
    std::vector<double> z; //arrays of depths of the boundaries inside the model
	
	unsigned int nTx, nRx; //number of transmitters and receivers
	std::vector<std::vector<double>> Tx_list, Rx_list; //their locations
	
	//Tx parameters
	std::vector<double> moment_list, azimuthTx_list, dipTx_list;
	
	 //for computing the Integral
    const unsigned int nIntervals = 40; //number of intervals 40

    //The following vectors are precomputed to save time
    std::vector<double> xInt;//get points at which J1(x) = 0;
    std::vector<double> Bx;
    std::vector<double> J0xW;
    std::vector<double> J1xW;
	
	// the thicknesses of the layers
	std::vector<double> h;
	
	double PI;
	double mu0; 
	
	// pre-computation 
	// xInt, Bx, J0xW, J1xW for integral evaluation
	void  pre_compute();
	
	private:
	
	//function returns: quadrature intervals (real coords) and Bessel*weight at each quadrature point
    void getBesselxW();

    //function returns nIntervals discrete points in a vector x = lambda *rho
    void getBesselZeros();
};
#endif
