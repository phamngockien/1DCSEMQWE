//CLASS FOR READ INPUT FILE, RUN THE 1D CSEM MODELING AND WRITE RESUTLS IN OUTPUT FILE

#ifndef RUN1DCSEM_H
#define RUN1DCSEM_H

#include <fstream>
#include <string>
#include <algorithm> 	//for find_if
#include <vector>
#include <exception>
#include <complex>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <limits> //for infinity

#include "Prms.h"
#include "GET_CSEM1D_FD_QWE.h"
#include "output_1D_CSEM_solution.h"

class Run1DCSEM
{
	public:
	
	Run1DCSEM();
	
	void Run1DForwardProblem();
	
	private:
	
	void read_input_file(PRMS &prms);
	
	void print_input(PRMS &prms) const;
	
	void read_parameters_line(std::string &line, //input
												bool &comment_line,
												std::string &parameter_name,
												std::string &parameter_value);

	void read_Tx_parameters(std::ifstream &file, //input
									std::string &line, //input
									PRMS &prms);
												
	void read_Rx_list(std::ifstream &file, //input
									std::string &line, //input
									PRMS &prms);
	
	void read_freq_list(std::ifstream &file, //input
									std::string &line, //input
									PRMS &prms);
											
	void read_layers_parmeters(std::ifstream &file, //input
												   std::string &line, //input
												   PRMS &prms); 
												
										 
	static inline void replace_tab_by_space(char &c);
	
	//trim a string
	static inline void ltrim(std::string &s);
	static inline void rtrim(std::string &s);
	static inline void trim(std::string &s);
	
	void tokenize(const std::string  &str, 
						  const char &delim,
						  std::vector<std::string> &out);
						  
						  
	//Testing print
	void print_point_vec(const auto &vec) const;
	void print_vec(const auto &vec) const;
	void print_Tx_parmeters(PRMS &prms) const;
};

//declare the inline functions
inline void Run1DCSEM::replace_tab_by_space(char &c){
	if (c == '\t')  c =  ' ';
}

// trim from start (in place)
inline void Run1DCSEM::ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void Run1DCSEM::rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
inline void Run1DCSEM::trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}
#endif