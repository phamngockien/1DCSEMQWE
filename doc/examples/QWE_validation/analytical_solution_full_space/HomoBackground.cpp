#include "HomoBackground.h"

//this class computes the analytical solution for either x-oriented ED or z-oriented MD in a homogeneous full space
//number of transmitter is 1
//receivers' positions are given by the input

HomoBackground::HomoBackground(const unsigned int &dipole_type,
																const double &frequency,
																const double &moment,
																const std::vector<double> &Tx,
																const std::vector<std::vector<double>> &Rx_list,
																const double &sig)
{
	PI = std::acos(-1);
	mu0 = 4. * PI * 1e-7;
	omega = 2.* PI * frequency;
	
	const unsigned int n_points = Rx_list.size();
	std::vector<std::vector<double>> e_values_list(n_points, std::vector<double>(6));
	std::vector<std::vector<double>> h_values_list(n_points, std::vector<double>(6));
	
	if (dipole_type == 0) { // x-oriented ED
		
		//call the function for compute analytical solution for x-oriented ED
		for (unsigned int i = 0; i < n_points; ++i) {
			e_values_list[i] = e_values_ED(moment, Tx, Rx_list[i],sig);
			h_values_list[i] = h_values_ED(moment, Tx, Rx_list[i],sig);
		}
		
	}else if (dipole_type == 1) { // z-oriented MD
		
		//call the function for compute analytical solution for z-oriented MD
		for (unsigned int i = 0; i < n_points; ++i) {
			e_values_list[i] = e_values_MD(moment, Tx, Rx_list[i],sig);
			h_values_list[i] = h_values_MD(moment, Tx, Rx_list[i],sig);
		}
	}else{
			std::cout<<"the code now only work for x-oriented ED (dipole_type 0) or z-oriented MD (dipole_type 1) : \n";
			throw std::exception();
	}
	
	//save the results in a csv file
	std::ofstream file;
	file.open("analytical_solution.csv", std::ofstream::trunc);
	
	if (file){ //test file openning is good
		 file.precision(32);
		  //Header line
		file << "source_type,TX_x,TX_y,TX_z,moment,frequency,sig,";
		file << "RX_x,RX_y,RX_z,";
		file << "Ex_re" <<"," << "Ex_im" << "," << "Ey_re" <<"," << "Ey_im" <<"," << "Ez_re" <<"," << "Ez_im" << "," 
			 << "Hx_re" <<"," << "Hx_im" <<","  << "Hy_re" <<","  << "Hy_im" <<"," << "Hz_re" <<","  << "Hz_im" <<"," ;
		file << "\n";
		

		for (unsigned int j = 0; j <  n_points; ++j) {
			file << dipole_type <<"," << Tx[0] <<"," << Tx[1] <<"," << Tx[2] <<","
				 << moment <<","  << frequency <<","  << sig <<","
				 << Rx_list[j][0] <<"," << Rx_list[j][1] <<"," << Rx_list[j][2] <<"," 
				 << e_values_list[j][0] << "," << e_values_list[j][1] <<","  // Ex -re, imag
				 <<  e_values_list[j][2] <<"," << e_values_list[j][3] << ","   // Ey  -re, imag
				 <<  e_values_list[j][4] << "," << e_values_list[j][5] <<","  // Ez  -re, imag
				<< h_values_list[j][0] << "," << h_values_list[j][1] <<","  // Hx -re, imag
				 <<  h_values_list[j][2] <<"," << h_values_list[j][3] << ","   // Hy  -re, imag
				 <<  h_values_list[j][4] << "," << h_values_list[j][5] ;  // Hz  -re, imag
			file << std::endl;
		}
		 
		 file.close();
	}else {
		std::cout<<"Could not open the output file, please check the HomoBackground constructor !";
		std::cout<< "\n\n";
		throw std::exception();	
	}		
		
}


HomoBackground::~HomoBackground() = default;

//==========================================================================================
// electric dipole source  x-oriented
std::vector<double> HomoBackground::e_values_ED(const double &moment,
																									const std::vector<double> &Tx,
																									const std::vector<double> &Rx,
																									const double &sig) const
{
	 std::complex<double> i(0, 1);

    const std::complex<double> k = sqrt(- i * omega * mu0  * sig);

    double x = (Rx[0] - Tx[0]);
    double y = (Rx[1] - Tx[1]);
    double z = (Rx[2] - Tx[2]);
	double r = sqrt(x*x + y*y + z*z);
	//prevent nan values
    if (r < 1e-6) {r = 1e-6;}

    std::vector<double> values(6);

   
    std::complex<double> E_x = moment * std::exp(-1. * i * k * r)
													* ((x * x / r / r) * (-k * k * r * r + 3. + 3. * i * k * r) +
														(k * k * r * r - i * k * r - 1.))
													/ (4 * PI * sig * r * r * r);
    
	std::complex<double> E_y = moment * std::exp(-1. * i * k * r)
													* (x * y / r / r) * (-k * k * r * r + 3. + 3. * i * k * r)   
													/ (4 * PI * sig * r * r * r);
    
	std::complex<double> E_z = moment * std::exp(-1. * i * k * r)
													* (x * z / r / r) * (-k * k * r * r + 3. + 3. * i * k * r)   
													/ (4 * PI * sig * r * r * r);

    values[0] = E_x.real();
	values[1] = E_x.imag();
    values[2] = E_y.real();
	values[3] = E_y.imag();
    values[4] = E_z.real();
    values[5] = E_z.imag();
    
    return values; 
}


std::vector<double> HomoBackground::h_values_ED(const double &moment,
																									const std::vector<double> &Tx,
																									const std::vector<double> &Rx,
																									const double &sig) const
{
	 std::complex<double> i(0, 1);

    const std::complex<double> k = sqrt(- i * omega * mu0  * sig);

    double x = (Rx[0] - Tx[0]);
    double y = (Rx[1] - Tx[1]);
    double z = (Rx[2] - Tx[2]);
	double r = sqrt(x*x + y*y + z*z);
	//prevent nan values
    if (r < 1e-6) {r = 1e-6;}

    std::vector<double> values(6);

        std::complex<double> H_x = 0.;
        std::complex<double> H_y = moment * std::exp(-1. * i * k * r)
														* (1. + i * k * r) * (-z)
														/ (4 * PI * r * r * r);
        std::complex<double> H_z = moment * std::exp(-1. * i * k * r)
														* (1. + i * k * r) * (y)
														/ (4 * PI * r * r * r);

    values[0] = H_x.real();
	values[1] = H_x.imag();
    values[2] = H_y.real();
	values[3] = H_y.imag();
    values[4] = H_z.real();
    values[5] = H_z.imag();
    
    return values;
}

//==========================================================================================
// magetic dipole source    +z-oriented

std::vector<double> HomoBackground::e_values_MD(const double &moment,
																									const std::vector<double> &Tx,
																									const std::vector<double> &Rx,
																									const double &sig) const
{
	std::complex<double> i(0, 1);

    const std::complex<double> k = sqrt(- i * omega * mu0  * sig);

    double x = (Rx[0] - Tx[0]);
    double y = (Rx[1] - Tx[1]);
    double z = (Rx[2] - Tx[2]);
	double r = sqrt(x*x + y*y + z*z);
	//prevent nan values
    if (r < 1e-6) {r = 1e-6;}

    std::vector<double> values(6);

        std::complex<double> E_x = i * omega * mu0 * moment * (1. + i * k * r)
														* std::exp(-1. * i * k * r)
														* (y / r )
														/ (4 * PI * r * r);
        std::complex<double> E_y = i * omega * mu0 * moment * (1. + i * k * r)
														* std::exp(-1. * i * k * r)
														* (-x / r )
														/ (4 * PI * r * r);
        std::complex<double> E_z = 0;

		values[0] = E_x.real();
		values[1] = E_x.imag();
		values[2] = E_y.real();
		values[3] = E_y.imag();
		values[4] = E_z.real();
		values[5] = E_z.imag();
    

    return values;
}					


std::vector<double> HomoBackground::h_values_MD(const double &moment,
																									const std::vector<double> &Tx,
																									const std::vector<double> &Rx,
																									const double &sig) const
{		
	 std::complex<double> i(0, 1);

    const std::complex<double> k = sqrt(- i * omega * mu0  * sig);

    double x = (Rx[0] - Tx[0]);
    double y = (Rx[1] - Tx[1]);
    double z = (Rx[2] - Tx[2]);
	double r = sqrt(x*x + y*y + z*z);
	//prevent nan values
    if (r < 1e-6) {r = 1e-6;}

    std::vector<double> values(6);

    std::complex<double> H_x = moment * std::exp(-1. * i * k * r)
													* ((x*z/r/r)*(-k*k*r*r+3.+3.*i*k*r))
													/ (4 * PI * r * r * r);

    std::complex<double> H_y = moment * std::exp(-1. * i * k * r)
													* ((y*z/r/r)*(-k*k*r*r+3.+3.*i*k*r))
													/ (4 * PI * r * r * r);

    std::complex<double> H_z = moment * std::exp(-1. * i * k * r)
													* ((z*z/r/r)*(-k*k*r*r+3.+3.*i*k*r)+(k*k*r*r-1.-i*k*r))
													/ (4 * PI * r * r * r);


    values[0] = H_x.real();
	values[1] = H_x.imag();
    values[2] = H_y.real();
	values[3] = H_y.imag();
    values[4] = H_z.real();
    values[5] = H_z.imag();
    
    return values;
}												  
