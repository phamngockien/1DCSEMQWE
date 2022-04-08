//
// Class to output the results in a file
//

#ifndef OUTPUT_1D_CSEM_SOLUTION_H
#define OUTPUT_1D_CSEM_SOLUTION_H

#include <vector>
#include <complex>
#include <fstream>
#include <iostream>
#include <exception>
#include <cmath>
#include <string>



#include "Prms.h"

class output_1D_CSEM_solution {
public:

    output_1D_CSEM_solution(PRMS &prms);

    virtual ~output_1D_CSEM_solution();
    
	//function to save the model parameters
	void save_model_parameters() const;
	
	//function to save the TX-RX configurations
	void save_TX_RX_configs() const;
	
	//function to save data in blocks for a sperate frequency
	void save_block(const unsigned int & freq_index) const;
	
	//function to save the EM fields at all receivers induced by the Tx_index source at the freq_index frequency
	void save(const unsigned int & freq_index,
					 const std::vector<std::vector<std::complex<double>>> &EM_fields_1D) const;
								
	private:
	//pointer to prms object
	PRMS *prms_ptr;
};


#endif //OUTPUT_1D_CSEM_SOLUTION_H
