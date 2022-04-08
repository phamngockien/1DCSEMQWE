//
// Class to output the results in a file
//

#include "output_1D_CSEM_solution.h"

output_1D_CSEM_solution::output_1D_CSEM_solution(PRMS &prms){
	//save the address of prms object 
	prms_ptr = &prms;

	//delete al previous csv files those store the old data
	system("rm -rf *.csv");
}

output_1D_CSEM_solution::~output_1D_CSEM_solution() = default;

//function to save the input parameters
void output_1D_CSEM_solution::save_model_parameters()  const
{
	std::ofstream file;
	
	file.open("model_parameters.csv", std::ofstream::app);
	
	if (file){ //test file openning is good
		 file.precision(32);
		
		//save the Layers' parameters
		//file << "Number_of_layer," << prms.nLayers << std::endl;
		file << "Top_depth(m),resistivity(ohm-m)" << std::endl;
		for (unsigned int k = 0; k< prms_ptr->sig.size(); ++k){
			file << prms_ptr->z[k] << "," << 1./prms_ptr->sig[k] << "\n";
		}
		file << "\n";
		
		file.close();
	}else {
		std::cout<<"Could not open the output file, please check the save_model_parameters function !";
		std::cout<< "\n\n";
		throw std::exception();	
	}		
}


//function to save the TX-RX configurations
void output_1D_CSEM_solution::save_TX_RX_configs() const
{
	std::ofstream file;
	file.open("TX_RX_configurations.csv", std::ofstream::app);
	
	if (file){ //test file openning is good
		 file.precision(32);
		 
		 //Header line
		file << "source_type,TX_x,TX_y,TX_z,moment,azimuth,dip,";
		file << "RX_x,RX_y,RX_z";
		file << "\n";
		
		//Save the data
		//loop over the list of transmitters (TX_list)
		for (unsigned int Tx_index = 0; Tx_index < prms_ptr->nTx ; ++Tx_index){
			
			//loop over the list of receivers (RX_list)
			for (unsigned int Rx_index = 0; Rx_index < prms_ptr->nRx ; ++Rx_index){
				//first column is the source_type
				switch (prms_ptr->dipole_type){
					case 0: file << "Electric_Dipole_Source," ;  break;
					case 1: file << "Magnetic_Dipole_Source,";  break;
					default: std::cout<<" The source type in the input file is wrong! \n\n"; throw std::exception();	
				}
				
				//then TX parameters	
				file << prms_ptr->Tx_list[Tx_index][0] << ","  << prms_ptr->Tx_list[Tx_index][1] << "," << prms_ptr->Tx_list[Tx_index][2] << ","  
				 <<prms_ptr->moment_list[Tx_index] <<"," <<prms_ptr->azimuthTx_list[Tx_index] <<","  << prms_ptr->dipTx_list[Tx_index] << ",";
				
				//finally the positions of the receivers
				file << prms_ptr->Rx_list[Rx_index][0] << ","  << prms_ptr->Rx_list[Rx_index][1] << "," << prms_ptr->Rx_list[Rx_index][2];
				file << "\n";
			}
		}			
	
	file.close();
	}else {
		std::cout<<"Could not open the output file, please check the save_TX_RX_configs function !";
		std::cout<< "\n\n";
		throw std::exception();	
	}		
}


//function to save data in blocks for a sperate frequency
void output_1D_CSEM_solution::save_block(const unsigned int & freq_index)  const
{
	std::ofstream file;
	std::string filename = "CSEM_1D_RESPONSES_at_frequency_index_" + std::to_string(freq_index+1) + ".csv";

	file.open(filename, std::ofstream::app);
	
	if (file){ //test file openning is good
		 file.precision(32);		

		//create labels to save the fields at each receiver
		file << "freq_Hz,"   
			 << "Ex_re" <<"," << "Ex_im" << "," << "Ey_re" <<"," << "Ey_im" <<"," << "Ez_re" <<"," << "Ez_im" << "," 
			 << "Hx_re" <<"," << "Hx_im" <<","  << "Hy_re" <<","  << "Hy_im" <<"," << "Hz_re" <<","  << "Hz_im" <<"," 
			 << "Bx_re" << "," << "Bx_im" <<","  << "By_re" <<"," << "By_im" << "," << "Bz_re" << "," 
			 << "Bz_im"; //to compare with Key (2009)
		file << std::endl;
		 
		 file.close();
	}else {
		std::cout<<"Could not open the output file, please check the save_block function!";
		std::cout<< "\n\n";
		throw std::exception();	
	}								
}

//function to save the EM fields at all receivers induced by the Tx_index source at the freq_index frequency
void output_1D_CSEM_solution::save(const unsigned int &freq_index,
															   const std::vector<std::vector<std::complex<double>>> &EM_fields_1D)   const
{

	std::ofstream file;
	std::string filename = "CSEM_1D_RESPONSES_at_frequency_index_" + std::to_string(freq_index+1) + ".csv";
	file.open(filename, std::ofstream::app);
	
	//observation frequency
	double freq = prms_ptr->frequencies[freq_index];
	
	if (file){ //test file openning is good
		 file.precision(32);
			
		//save the 1D CSEM solution
		for (unsigned int j = 0; j < prms_ptr->nRx; ++j) {
			file << freq <<","
				 << EM_fields_1D[j][0].real() << "," << EM_fields_1D[j][0].imag() <<","  // Ex
				 << EM_fields_1D[j][1].real() <<"," << EM_fields_1D[j][1].imag() << ","   // Ey
				 << EM_fields_1D[j][2].real() << "," << EM_fields_1D[j][2].imag() <<","  // Ez
				 << EM_fields_1D[j][3].real() << "," << EM_fields_1D[j][3].imag() << ","   // Hx
				 << EM_fields_1D[j][4].real() << "," << EM_fields_1D[j][4].imag() <<","   // Hy
				 << EM_fields_1D[j][5].real() << "," << EM_fields_1D[j][5].imag() << ","   // Hz
				 << EM_fields_1D[j][3].real()   * 4. * prms_ptr->PI * 1e-7 << "," 
				 << EM_fields_1D[j][3].imag() * 4. * prms_ptr->PI * 1e-7 << ","   // Bx
				 << EM_fields_1D[j][4].real()   * 4. * prms_ptr->PI * 1e-7 << "," 
				 << EM_fields_1D[j][4].imag() * 4. * prms_ptr->PI * 1e-7 << ","   // By
				 << EM_fields_1D[j][5].real()   * 4. * prms_ptr->PI * 1e-7 << "," 
				 << EM_fields_1D[j][5].imag() * 4. * prms_ptr->PI * 1e-7;        // Bz
			file << std::endl;
		}
			
		file.close();
	}else {
		std::cout<<"Could not open the output file, please check the save function!";
		std::cout<< "\n\n";
		throw std::exception();	
	}		
}
