#ifndef HOMOBACKGROUND_H
#define HOMOBACKGROUND_H


#include <vector>
#include <complex>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>


//compute the back ground field in homogeneous media
// (sources in unbounded medium)
class HomoBackground
{
public:
	
 HomoBackground(const unsigned int &dipole_type,
							  const double &frequency,
							const double &moment,
							const std::vector<double> &Tx,
							const std::vector<std::vector<double>> &Rx_list,
							const double &sig);
							
virtual ~HomoBackground();
	
private:

	double PI;
	double mu0;
	double omega;
	
	// electric dipole source   +x-oriented
	std::vector<double> e_values_ED(const double &moment,
															const std::vector<double> &Tx,
															const std::vector<double> &Rx,
															const double &sig) const;
												
	std::vector<double> h_values_ED(const double &moment,
															const std::vector<double> &Tx,
															const std::vector<double> &Rx,
															const double &sig) const;
																									
	
	//--------------------------------------------------------------------------------------------
	// magetic dipole source    +z-oriented
	std::vector<double> e_values_MD(const double &moment,
															const std::vector<double> &Tx,
															const std::vector<double> &Rx,
															const double &sig) const;
	
	std::vector<double> h_values_MD(const double &moment,
															const std::vector<double> &Tx,
															const std::vector<double> &Rx,
															const double &sig) const;
};

#endif