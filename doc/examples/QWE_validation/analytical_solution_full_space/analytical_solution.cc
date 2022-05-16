
#include "HomoBackground.h"

int main()
{
	//x-oriented ED : 0, z-oriented MD: 1
	unsigned int dipole_type = 0;
	
	//frequency (Hz)
	double frequency = 1.;
	double moment = 1.;
	double sig = 1.;
	
	//TX position
	std::vector<double> Tx(3);
	Tx[0] = 0.; Tx[1] = 0.; Tx[2] = 0.; //at the origin (0,0,0)
	
	//list of Rx positions
	unsigned int n_points = 31;
	std::vector<std::vector<double>> Rx_list(n_points, std::vector<double>(3));
	for (unsigned int k = 0; k < n_points; ++k){
		Rx_list[k][0] = 300.;
		Rx_list[k][1] = 0.;
		Rx_list[k][2] = 0. + k*100.;
	}
	
	HomoBackground analytical_solution(dipole_type,
																frequency,
																moment,
																Tx,
																Rx_list,
																sig);
	
	return 0;
}


