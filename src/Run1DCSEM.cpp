//CLASS FOR READ INPUT FILE, RUN THE 1D CSEM MODELING AND WRITE RESUTLS IN OUTPUT FILE

#include "Run1DCSEM.h"

Run1DCSEM:: Run1DCSEM(){}
 
void Run1DCSEM::Run1DForwardProblem()
{
	// initiate a class for all necessary parameters
	PRMS prms;
	
	//Input: load all parameters from input file
	read_input_file(prms);
	print_input(prms);
	
	// pre-computation 
	// PI, mu
	// xInt, Bx, J0xW, J1xW for integral evaluation by QWE
	prms.pre_compute();
	
	//an object to save the results
	output_1D_CSEM_solution out(prms);
	
	//save the model's parameters
	out.save_model_parameters();
	
	//TODO save the transmitter-receiver configurations
	out.save_TX_RX_configs();

	//Timer
	auto start = std::chrono::steady_clock::now();
	
	//Compute EM Fields and Output
	//loop over the frequencies
	for (unsigned int m = 0; m < prms.nfqs; ++m){
		
		//to store the EM fields in seperate blocks w.r.t frequency
		out.save_block(m);
		
		//loop over transmitters
		for (unsigned int k = 0; k < prms.nTx; ++k){	
		
			//get the solution at	the m_th frequency for the k_th transmitter
			std::vector<std::vector<std::complex<double>>> EMfields;
			GET_CSEM1D_FD_QWE get_EM_Fields(prms, k, m);
			
			//solve for this Tx at this frequency
			get_EM_Fields.solve_1D_CSEM(EMfields);	
			
			//save the responses to output file
			out.save(m, EMfields);	
		}
	}

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Time for getting the EM field is: " << elapsed_seconds.count() << "seconds \n";	
}

//=======================================================================
//private functions

void Run1DCSEM::read_input_file(PRMS &prms)
{
	std::ifstream file("1D_CSEM_INPUT.txt");  //open a file to perform read operation using file object
   if (file) {  //same as file.good()
		std::string line;
		while(getline(file,line)){ //read data from file object and put it into string.
			bool comment_line;
			std::string parameter_name, parameter_value;
			read_parameters_line(line, comment_line, parameter_name, parameter_value);
			
			if (!comment_line){
				if (parameter_name == "# dipole type"){  //read the parameter and its value
					prms.dipole_type = std::stoi(parameter_value); 
				    
					//test the dipole type if >1 or <0 then it is not implemented
					if (prms.dipole_type > 1 || prms.dipole_type < 0) 
					{
						std::cout<<"dipole type must be in the range [0,1] ! Please change in the input file";
						std::cout<< "\n\n";
						throw std::exception();	
					}
				}
						
				if (parameter_name == "# nquad"){
					prms.nQuad = std::stoi(parameter_value); 
				}
				
				if (parameter_name == "# rtol"){
					prms.rtol = std::stod(parameter_value); 
				}
				
				if (parameter_name == "# atol"){
					prms.atol = std::stod(parameter_value); 
				}
			
				if (parameter_name == "# transmitters"){
					prms.nTx = std::stoi(parameter_value); 
					read_Tx_parameters(file, line, prms);
				}

				if (parameter_name == "# frequencies"){
					prms.nfqs = std::stoi(parameter_value); 
					read_freq_list(file, line, prms);
				}
			
				if (parameter_name == "# layers"){
					prms.nLayers = std::stoi(parameter_value); 
					read_layers_parmeters(file, line,prms);
				}
					
				if (parameter_name == "# receivers"){
					prms.nRx = std::stoi(parameter_value); 
					read_Rx_list(file, line, prms);
				}			
			}	
		}
		file.close(); //close the file object.
	} else {
		std::cout<<"Could not open the input file!";
		std::cout<< "\n\n";
		throw std::exception();	
	}		
}


void Run1DCSEM::print_input(PRMS &prms) const
{
	std::cout << "===========================================================================================\n";
	std::cout << "                                    1D CSEM QWE																			\n";
	std::cout << "\n";
	std::cout << "   A program for computing the EM fields from an arbitrarily type of dipole source						\n";
	std::cout << "      in an N-layered model using vector potential formulation 						\n";
	std::cout << "\n";
	std::cout << "\n";
	std::cout << " 	         This is a free program, but WITHOUT ANY WARRANTY									\n";
	std::cout << " \n";
	std::cout << "                                 READING INPUT FILE																		\n";
	std::cout << "============================================================================================\n";
	
	//uncomment the following codes when you want to check the input parameters 
	/*
	std::cout << "\n";
	std::cout<<"The type of dipole: " <<prms.dipole_type << std::endl; 
	std::cout << "(Types: 0, 1 are electric and magnetic dipole, respectively). \n\n";
	std::cout<<"Number of points to compute integral by quadrature rule: " << prms.nQuad << std::endl;
	std::cout << "--------------------------------------------------------------------------------------------\n";
	
	std::cout << "\n";
	std::cout<<"Number of Transmitters: "<< prms.nTx <<",   including: \n";
	print_Tx_parmeters(prms);
	std::cout << "--------------------------------------------------------------------------------------------\n";
	
	std::cout << "\n";
	std::cout<<"Number of frequencies: "<< prms.nfqs <<",   including: \n";
	print_vec(prms.frequencies);
	std::cout << "--------------------------------------------------------------------------------------------\n";
	
	std::cout << "\n";
	std::cout<<"Number of layers: "<< prms.nLayers <<",   with parameters: \n";
	std::cout<<"Top depth \t conductivity \n";
	for (unsigned  int k = 0; k < prms.sig.size(); ++k){
		std::cout<< prms.z[k] << " \t  \t " << prms.sig[k] << std::endl;
	}
	std::cout<<prms.z.back() << " (base of the lower-most interface)" <<std::endl;
	std::cout << "--------------------------------------------------------------------------------------------\n";
		
	std::cout << "\n";
	std::cout<<"Number of Receivers: "<< prms.nRx <<",  including points: \n";
	print_point_vec(prms.Rx_list);
	std::cout << "--------------------------------------------------------------------------------------------\n";
	*/
	
	std::cout << "\n";
	std::cout<<"Computing the EM fields .... " << std::endl; 
}

void Run1DCSEM::read_parameters_line(std::string &line,
																	 bool &comment_line,
																	 std::string &parameter_name,
																	 std::string &parameter_value) 
{
	//replace all tabs by spaces
	 for_each (line.begin(), line.end(), replace_tab_by_space);
		
	//erase all spaces at the begining of the line
	ltrim(line);					
	
	comment_line = false;
     // if the line is blank, consider it a comment_line.
	 if (line.empty()) comment_line = true;
	 // If the first char is a comment_line char (! or %), then the whole line is a comment_line.
	 if (line[0] == '!' || line[0] == '%') comment_line = true;
	
	if (!comment_line){
		char delim; 
		
		//seperate code and comment_line parts in the line
		std::vector<std::string> code_comment;
		//for "!" signal
		if (line.find('!') != std::string::npos){
			delim = '!';
			tokenize(line, delim, code_comment);
			line = code_comment[0]; //code part is the first part
		}
		//for "%" signal
		if (line.find('%') != std::string::npos){
			delim = '%';
			tokenize(line, delim, code_comment);
			line = code_comment[0]; //code part is the first part
		}
			
		//seperate parameter name and its value
		//std::size_t found = line.find_first_of("+-0123456789."); //This is the value of the parameter or model data
		std::vector<std::string> name_value;
		
		//NOTE the data of the parameter is read in different way in read_input_file functions
		if (line.find(':') == std::string::npos /*&& found != 0*/) { //The line does not contain parameter name and value
				std::cout << "Reading input file faills! Please check from this part: \n";
				std::cout << "There must be a ':' seperating the parameter name and its value at line: \n"<< line << "\n\n";
				throw std::exception();
		} else if (line.find(':') != std::string::npos){ //The line contains parameter name and value
				delim = ':';
				tokenize(line, delim, name_value);
				rtrim(name_value[0]); //right trim the name
				trim(name_value[1]); //trim from both ends the value
				parameter_name = name_value[0];
				parameter_value = name_value[1];
				//lower case the parameter_name
				std::transform(parameter_name.begin(), parameter_name.end(), parameter_name.begin(), [](unsigned char c){ return std::tolower(c); });
		}
	}  	
}

void Run1DCSEM::read_Tx_parameters(std::ifstream &file, //input
																	std::string &line, //input
																	PRMS &prms)
{
	std::vector<double> temp_point(3); //for x, y, z positions
	double moment, azimuth, dip;
	for (unsigned int k = 0; k < prms.nTx; ++k){
		if (file >> temp_point[0] >> temp_point[1] >> temp_point[2] >> moment >> azimuth >> dip){
			
			//test the dipTx in the range [-90, 90] degree, 0 value is for the (Oxy) plane
			if (dip > 90. || dip < -90.) 
			{
				std::cout<<"the dipTx of the dipole source must be in the range [-90,90] ! \n";
				std::cout<< "Please change in the input file from this part \n";
				std::cout<< line <<"\n\n";
				throw std::exception();	
			}
			
			//test the azimuthTx in the range [0, 360] degree, 0 value is for the Ox direction
			if (azimuth > 360. || azimuth < 0.) 
			{
				std::cout<<"the azimuthTx of the dipole source must be in the range [0, 360] ! \n";
				std::cout<< "Please change in the input file from this part \n";
				std::cout<< line <<"\n\n";
				throw std::exception();	
			}
			
			prms.Tx_list.push_back(temp_point);
			prms.moment_list.push_back(moment);
			prms.azimuthTx_list.push_back(azimuth);
			prms.dipTx_list.push_back(dip);
		}else{
			std::cout<<"Reading Input fails! Please check from this part: \n";
			std::cout<< line <<"\n\n";
			throw std::exception();
		}	
	}									
}

void Run1DCSEM::read_Rx_list(std::ifstream &file,
														 std::string &line,
														 PRMS &prms)
{
	std::vector<double> temp_point(3);
	for (unsigned int k = 0; k < prms.nRx; ++k){
		if (file >> temp_point[0] >> temp_point[1] >> temp_point[2]){
			prms.Rx_list.push_back(temp_point);
		}else{
			std::cout<<"Reading Input fails! Please check from this part: \n";
			std::cout<< line <<"\n\n";
			throw std::exception();
		}	
	}
}

void Run1DCSEM::read_freq_list(std::ifstream &file, //input
														std::string &line, //input
														PRMS &prms)
{
	double temp;
	for (unsigned int k = 0; k < prms.nfqs; ++k){ //ignore parameter_name and load data
		if(file >> temp){		//skip header and load data lines
			prms.frequencies.push_back(temp);
		}else{
			std::cout<<"Reading Input fails! Please check from this part: \n";
			std::cout<< line <<"\n\n";
			throw std::exception();
		}
	}
}

void Run1DCSEM::read_layers_parmeters(std::ifstream &file, //input
																		std::string &line,  //input
																		PRMS &prms)
{
	double temp_z, temp_rho;
	for (unsigned int k = 0; k < prms.nLayers; ++k){
		if(file >> temp_z >> temp_rho){
			prms.z.push_back(temp_z);
			prms.sig.push_back(1./temp_rho);
		}else{
			std::cout<<"Reading layers' parameters fails! Please check from this part: \n";
			std::cout<< line <<"\n\n";
			throw std::exception();
		}
	}
	prms.z[0] = -1. * std::numeric_limits<double>::max(); //the uppermost interface is infinity
	prms.z.push_back(std::numeric_limits<double>::max()); //the lowermost interface is infinity	
	
	 //find thickness of each layer
	//set the thickness of the top layer to infinity
	prms.h.push_back(std::numeric_limits<double>::max());
	//compute the thicknesses of inner layers
	for (unsigned int j = 2; j < prms.z.size(); ++j) {
		prms.h.push_back(prms.z[j] - prms.z[j-1]);
		if (prms.h.back() <=0){
			std::cout << "There is an error in the input of depth in the INPUT file !!! \n\n";
			throw std::exception();	
		}
	}
	//set the thickness of the bottom layer to infinity
	prms.h.back() = std::numeric_limits<double>::max();
	
	if (prms.h.size() != prms.sig.size()){
		std::cout << "There is an error in the computation of thicknesses !!! \n\n";
		throw std::exception();	
	}
}

//split a string by delemiter
void Run1DCSEM::tokenize(const std::string &str, 
												const char &delim,
												std::vector<std::string> &out)
{
	size_t start;
	size_t end = 0;

	while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
	{
		end = str.find(delim, start);
		out.push_back(str.substr(start, end - start));
	}
}


//Testing print
void Run1DCSEM::print_point_vec(const auto &vec) const
{
	for (auto &v : vec){
		std::for_each(v.begin(), v.end(), [](double i){std::cout << i << "\t";});
		std::cout <<std::endl;
	}
}

void Run1DCSEM::print_vec(const auto &vec) const
{
	std::for_each(vec.begin(), vec.end(), [](auto i){std::cout << i << "\n";});
}

void Run1DCSEM::print_Tx_parmeters(PRMS &prms) const
{
	for (unsigned int k = 0; k < prms.Tx_list.size(); ++k){
		std::cout << "Transmitter " << k+1 << " located at: \t" 
					 << prms.Tx_list[k][0] << "\t" << prms.Tx_list[k][1] << "\t"<< prms.Tx_list[k][2] << "\n";
		std::cout << "that has the moment: " << prms.moment_list[k] 
					  <<", \t azimuth:  " << prms.azimuthTx_list[k]
					  <<", \t dip:  " << prms.dipTx_list[k] <<"\n\n";
	}
}