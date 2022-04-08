//
// Class to get EM fields for a given Transmitter at a certain frequency value
//

#include "GET_CSEM1D_FD_QWE.h"


//the code implements for Electric Dipole (dipole type 0), Magnetic Dipole source (source type 1)
GET_CSEM1D_FD_QWE::GET_CSEM1D_FD_QWE(PRMS &prms,
																						   unsigned int & Tx_index,
																						   unsigned int & freq_index){
	//save the address of prms object 
	prms_ptr = &prms;
	Tx_ptr = &Tx_index;
	freq_ptr = &freq_index;
}

GET_CSEM1D_FD_QWE::~GET_CSEM1D_FD_QWE() = default;


void GET_CSEM1D_FD_QWE::solve_1D_CSEM(std::vector<std::vector<std::complex<double>>> &EMfields) {

    setup_1d_CSEM();

    //*******************************************************************
    //for each receiver we implement the following codes
    for (const auto &receiver : prms_ptr->Rx_list) {
        //find the layer has receiver
        RxLayer = 0;
        for (unsigned int j = 1; j < prms_ptr->nLayers; ++j) {
            if (receiver[2] > prms_ptr->z[j]) {
                RxLayer = j;
            }
        }

        //compute x-,y-direction difference between the current receiver and the source
        double dx, dy;
        dx = receiver[0] - prms_ptr->Tx_list[*Tx_ptr][0];
        dy = receiver[1] - prms_ptr->Tx_list[*Tx_ptr][1];
		
        double rho = sqrt(dx * dx + dy * dy);
		singularity = false;
		//prevent rho = 0 leads to explosion (set rho to be unit length)
		if (rho < 1.) {
			rho = 1.;
			singularity = true;
			if (prms_ptr->azimuthTx_list[*Tx_ptr] == 0) dx = rho; //has only x-component for the EM fields
		}
		
		//solve the 1D forward problem
		(prms_ptr->dipole_type == 0) ? Solve_ED(dx, dy, rho, receiver, EMfields) : Solve_MD(dx, dy, rho, receiver, EMfields);
    }
	//***************************************************************************
}

//###############################################################################################
//pre-computing part
//###############################################################################################
void GET_CSEM1D_FD_QWE::setup_1d_CSEM() {

	omega = 2.0 * prms_ptr->PI * prms_ptr->frequencies[*freq_ptr]; 
	
    //find the source layer
    //index of the source layer (index range:0 to nLayers-1)
    TxLayer = 0;
    for (unsigned int j = 1; j < prms_ptr->nLayers; ++j) {
        if (prms_ptr->Tx_list[*Tx_ptr][2] > prms_ptr->z[j]) {
            TxLayer = j;
        }
    }
	//----------------------------------------------------------------------------------------
	
	//find the source_type need to compute
	//first reinitialize false for all source type
	hed = false; ved = false; hmd = false; vmd = false;
	if (prms_ptr->dipTx_list[*Tx_ptr] == 0.) { //only horizontal dipole source: HED (0) or HMD (2) source_type
		(prms_ptr->dipole_type == 0) ? hed = true : hmd = true;
	} else if (prms_ptr->dipTx_list[*Tx_ptr] == -90. || prms_ptr->dipTx_list[*Tx_ptr]  == 90.){  //only vertical dipole source: VED (1) or VMD (3) source_type
		(prms_ptr->dipole_type == 0) ?  ved = true : vmd = true;
	} else { //has both vertical and horizontal dipole source components
		switch (prms_ptr->dipole_type){
			case 0:   hed = true; ved = true; break;
			case 1:   hmd = true; vmd = true; break;
			default: std::cout << " there is error, invalid condition: 0 <= dipole_type <=1 \n"; break;
		}
	}

}

//###############################################################################################
//function to solve for Electric dipole source
void GET_CSEM1D_FD_QWE::Solve_ED(double &dx,
																		const double &dy,
																		const double &rho,
																		const std::vector<double> &receiver,
																		std::vector<std::vector<std::complex<double>>> &EMfields) const
{

    //get the results to fields vector for this receiver
    std::vector<std::complex<double>> fields(6);
		
	if (hed & ved){ //the dipole source has both horizontal and vertical components
		
		//vertical dipole formula works well with dx = dy = 0 and rho != 0, thus compute it first
		//compute the ved fields part - source_type 1
		std::vector<std::complex<double>> VED_fields(6);
		qwe(1, dx, dy, rho, receiver, VED_fields);

		//compute the hed fields part - source_type 0
		std::vector<std::complex<double>> HED_fields(6);
		
		if (prms_ptr->azimuthTx_list[*Tx_ptr] != 0) { // horizontal part of the dipole does not point in x-direction
			
			// convert dx and dy to those in the rotated Ox (dx_hat, dy_hat)
			double dx_hat, dy_hat;
			rotate_coords_comp(dx, dy,  rho, dx_hat, dy_hat);
			
			//compute the EM fields in the rotated coords
			qwe(0, dx_hat, dy_hat, rho, receiver, HED_fields);
			
			// rotate to get EM fields in real (Oxyz) coords
			// because we computed them in the rotated coords before
			horizontal_rotate(HED_fields);
		} else {
			
			//compute the EM fields in the (Oxyz) coords
			qwe(0, dx, dy, rho, receiver, HED_fields);
		}
			
		//add up the HED and VED to get the final responses
		add_horizontal_vertical_souce(HED_fields, VED_fields, fields);
			
		//store the fields at current receiver after sum the horizontal and vertical solutions
		EMfields.push_back(fields);
			
	} else if (hed & (! ved)){ //only HED source
		
		//compute the hed fields part - source_type 0
		if (prms_ptr->azimuthTx_list[*Tx_ptr] != 0) { // horizontal part of the dipole does not point in x-direction
			
			// convert dx and dy to those in the rotated Ox (dx_hat, dy_hat)
			double dx_hat, dy_hat;
			rotate_coords_comp(dx, dy, rho,  dx_hat, dy_hat);
			
			//compute the EM fields in the rotated coords
			qwe(0, dx_hat, dy_hat, rho, receiver, fields);
			
			// rotate to get EM fields in real (Oxyz) coords 
			// because we computed them in the rotated coords before
			horizontal_rotate(fields);

		} else { 	// x-oriented ED
			
			//compute the EM fields in the (Oxyz) coords
			qwe(0, dx, dy, rho, receiver, fields);

		}

		//store the fields at current receiver after sum the horizontal and vertical solutions
		EMfields.push_back(fields);
		
	} else if (ved & (! hed)){ //only VED source - no need to rotate and sum

		//compute the ved fields part - source_type 1
		qwe(1, dx, dy, rho, receiver, fields); 
		
		//store the fields at current receiver
		EMfields.push_back(fields);
	}
}

//###############################################################################################
//function to solve for Magnetic dipole source
void GET_CSEM1D_FD_QWE::Solve_MD(double &dx,
																		const double &dy,
																		const double &rho,
																		const std::vector<double> &receiver,
																		std::vector<std::vector<std::complex<double>>> &EMfields) const
{

    //get the results to fields vector for this receiver
    std::vector<std::complex<double>> fields(6);
		
	if (hmd & vmd){ //the dipole source has both horizontal and vertical components
		
		//vertical dipole formula works well with dx = dy = 0 and rho != 0, thus compute it first
		//compute the vmd fields part - source_type 3
		std::vector<std::complex<double>> VMD_fields(6);
		qwe(3, dx, dy, rho, receiver, VMD_fields);
		
		//compute the hmd fields part - source_type 2
		std::vector<std::complex<double>> HMD_fields(6);
		
		if (prms_ptr->azimuthTx_list[*Tx_ptr] != 0) { // horizontal part of the dipole does not point in x-direction
			
			// convert dx and dy to those in the rotated Ox (dx_hat, dy_hat)
			double dx_hat, dy_hat;
			rotate_coords_comp(dx, dy, rho, dx_hat, dy_hat);
			
			//compute the EM fields in the rotated coords
			qwe(2, dx_hat, dy_hat, rho, receiver, HMD_fields);
			
			// rotate to get EM fields in real (Oxyz) coords
			// because we computed them in the rotated coords before
			horizontal_rotate(HMD_fields);
		}else {
			
			//compute the EM fields in the (Oxyz) coords
			qwe(2, dx, dy, rho, receiver, HMD_fields);			
		}
			
		//add up the HMD and VMD to get the final responses
		add_horizontal_vertical_souce(HMD_fields, VMD_fields, fields);
			
		//store the fields at current receiver after sum the horizontal and vertical solutions
		EMfields.push_back(fields);
			
	} else if (hmd & (! vmd)){ //only HMD source
		
		//compute the hmd fields part - source_type 2
		if (prms_ptr->azimuthTx_list[*Tx_ptr] != 0) { // horizontal part of the dipole does not point in x-direction
			
			// convert dx and dy to those in the rotated Ox (dx_hat, dy_hat)
			double dx_hat, dy_hat;
			rotate_coords_comp(dx, dy, rho, dx_hat, dy_hat);
			
			//compute the EM fields in the rotated coords
			qwe(2, dx_hat, dy_hat, rho, receiver, fields);
			
			// rotate to get EM fields in real (Oxyz) coords
			// because we computed them in the rotated coords before
			horizontal_rotate(fields);
		} else {
			
			//compute the EM fields in the (Oxyz) coords
			qwe(2, dx, dy, rho, receiver, fields);
		}
			
		//store the fields at current receiver 
		EMfields.push_back(fields);
		
	} else if (vmd & (! hmd)){ //only VMD source - no need to rotate and sum
			
		//compute the vmd fields part - source_type 3
		qwe(3, dx, dy, rho, receiver, fields); 
		//store the fields at current receiver
		EMfields.push_back(fields);
	} 
}

//compute the dx_hat, dy_hat in rotated coords when the azimuthTx != 0
void GET_CSEM1D_FD_QWE::rotate_coords_comp(const double &dx,
																						const double &dy,
																						const double &rho,
																						double &dx_hat,
																						double &dy_hat) const
{
		// Angle from Ox to site (degree) - counter clock-wise positive: alpha
		double alpha;	
		// Angle from Ox_hat (rotated Ox) to site (radian)  - counter clock-wise positive: beta
		double beta;
		
		if (singularity) { // dx = dy = 0
			dx_hat = rho;
			dy_hat =  0;
		} else if ((dx == 0) && (dy != 0) ) {
			(dy == rho) ? alpha = 90. : alpha = -90.; // dy = rho or dy = -rho
			beta = (alpha - prms_ptr->azimuthTx_list[*Tx_ptr])  * prms_ptr->PI / (180.);
			
			// convert dx and dy to those in the rotated Ox (dx_hat, dy_hat)
			dx_hat = rho * std::cos(beta);
			dy_hat = rho * std::sin(beta);
		} else if ((dx != 0) && (dy == 0) ) {
			(dx > 0) ? alpha = 0. : alpha = 180.; // dx > 0 or dx < 0
			
			beta = (alpha - prms_ptr->azimuthTx_list[*Tx_ptr])  * prms_ptr->PI / (180.);
			
			// convert dx and dy to those in the rotated Ox (dx_hat, dy_hat)
			dx_hat = rho * std::cos(beta);
			dy_hat = rho * std::sin(beta);
		} else {  // dx != 0; dy != 0
			alpha = std::atan2(dy, dx) * 180. / prms_ptr->PI ;
			
			beta = (alpha - prms_ptr->azimuthTx_list[*Tx_ptr])  * prms_ptr->PI / (180.);
			
			// convert dx and dy to those in the rotated Ox (dx_hat, dy_hat)
			dx_hat = rho * std::cos(beta);
			dy_hat = rho * std::sin(beta);
		}
}

// horizontal rotation to the real (Oxyz)
// only need for the horizontal component of the dipole source
void GET_CSEM1D_FD_QWE::horizontal_rotate(std::vector<std::complex<double>> &fields) const
{
	//change azimuthTx angle to radian
	double rot_ang = prms_ptr->azimuthTx_list[*Tx_ptr] * prms_ptr->PI / 180.;
	
	//some temporary variables
	std::complex<double> tmpEx, tmpEy, tmpHx, tmpHy;
	
	//Fx = Fx' cos(rot_ang) - Fy' sin(rot_ang)
	tmpEx = fields[0] * std::cos(rot_ang) - fields[1] * std::sin(rot_ang);
	tmpHx = fields[3] * std::cos(rot_ang) - fields[4] * std::sin(rot_ang);
	
	//Fy = Fx' sin(rot_ang) + Fy' cos(rot_ang)
	tmpEy =  fields[0] * std::sin(rot_ang) + fields[1] * std::cos(rot_ang);
	tmpHy =  fields[3] * std::sin(rot_ang) + fields[4] * std::cos(rot_ang);
	
	fields[0] = tmpEx; fields[1] = tmpEy;
	fields[3] = tmpHx; fields[4] = tmpHy;
	
}

// add the fields of horizontal and vertical components of the dipole source
// to get the final EM responses
void GET_CSEM1D_FD_QWE::add_horizontal_vertical_souce(const std::vector<std::complex<double>> &horizontal_fields,
																										const std::vector<std::complex<double>> &vertical_fields,
																										std::vector<std::complex<double>> &fields) const
{
	//test the size of the vectors
	 if (fields.size() != horizontal_fields.size() || fields.size() != vertical_fields.size()) 
	{
		std::cout<<"ERROR in computing the EM fields function!!";
		std::cout<< "\n\n";
		throw std::exception();	
	}
	
	//change dipTx angle to radian
	double dip = prms_ptr->dipTx_list[*Tx_ptr] * prms_ptr->PI / 180.;
	
	for (unsigned int k = 0; k < horizontal_fields.size(); ++k){
		fields[k] = horizontal_fields[k] * std::cos(dip) + vertical_fields[k] * std::sin(dip);
	}														   
}


//###############################################################################################
//each receiver computing part
//###############################################################################################

//****************************************************************************************************
//function implements the Quadrature with Extrapolation method (Key,2012 and Weniger, 2003 -pp26)
//with some modifications
//TODO explain clearly the modifications
//returns the fields with fields[0:2]  for Ex, Ey, Ez, and fields[3:5]  for Hx, Hy, Hz
//fields[k] = kernels[k] with k ranging from 0 to 5
void GET_CSEM1D_FD_QWE::qwe(const unsigned int &source_type,
															  const double &dx,
															  const double &dy,
															  const double &rho,
															  const std::vector<double> &receiver,
															  std::vector<std::complex<double>> &fields) const{

    //complex number i
    const std::complex<double> i(0, 1);

    //define the smallest real and biggest real (Weniger, 2003 pp26)
    double realmin = 1e-150;
    double realmax = 1e150;

    //define the number of Kernels
    //Here it is 6 for x,y,z component of E and H fields
    const unsigned int nKernels = 6;

    //initialize some arrays automatically containing zero value for all elements
    //storing the epsilon and etrapolation terms for kenel k and at order n
    //for the recursion coefficients for the Epsilon algorithm
    std::vector<std::vector<std::complex<double>>> E(nKernels);

    // extrapolated result for each order of the expansion
    std::vector<std::vector<std::complex<double>>> extrap(nKernels);

    //for control the convergence
    double relErr; //relative error for each order
    double absErr; //absolute error for each order

    //array to control the convergence of the partial sum
    std::vector<bool> converged(nKernels); //automatically false for all kernels
    bool allconverged; //check if all the kernels are converged

    //vector stores the kernels (results) for the current partial integral
    std::vector<std::complex<double>> this_kernels;

    //the variables need to compute Shanks transform
    std::complex<double> aux1, aux2, diff;		
	
    //The extrapolation loop (epsilon algorithm)
    for (unsigned int n = 0; n < prms_ptr->nIntervals; ++n) {

        //Step 1: compute the current partial integral (Fn) to return this kernels
		
		//test the dipole type if >1 then it is not implemented
		if (source_type > 3 || source_type < 0) 
		{
			std::cout<<"source type must be in the range [0,3] ! ERROR OCCURS IN FINDING SOURCE TYPE LIST";
			std::cout<< "\n\n";
			throw std::exception();	
		}
	
		switch (source_type){
			case 0:  get_kernels_hed(n, dx, dy, rho, receiver, this_kernels); break;
			case 1:  get_kernels_ved(n, dx, dy, rho, receiver, this_kernels); break;
			case 2:  get_kernels_hmd(n, dx, dy, rho, receiver, this_kernels); break; 
			case 3:  get_kernels_vmd(n, dx, dy, rho, receiver, this_kernels); break;
			default: std::cout << " there is error, invalid condition: 0 <= source_type <=3 \n"; break;
		}
		
        //Step 2: Compute Shanks transformation for each kernel
        for (unsigned int k = 0; k < nKernels; ++k) {

            //skip this kernel since it's done
            if (converged[k]) {
                continue;
            }

            if (n == 0) {
                E[k].push_back(this_kernels[k]); //first E_0^0
                extrap[k].push_back(this_kernels[k]); //first extrapolate = first integral evaluation

            } else {

                //E_0^n
                E[k].push_back(E[k][n - 1] + this_kernels[k]);

                // Compute the Shanks transform using the Epsilon algorithm:
                aux2 = 0. + 0. * i;
                for (unsigned int j = n; j >= 1; --j) {
                    aux1 = aux2;
                    aux2 = E[k][j - 1];  //@ kernel k, E(j-1)
                    diff = E[k][j] - aux2;
                    (abs(diff) <= realmin) ? E[k][j - 1] = realmax : E[k][j - 1] = aux1 + 1. / diff;
                }

                //The extrapolated result at order n
                (n % 2 == 0) ? extrap[k].push_back(E[k][0]) : extrap[k].push_back(E[k][1]);
            }

            //Step 3: Analyze for convergence
            if (n > 0) {
                absErr = std::abs(extrap[k][n] - extrap[k][n - 1]);

			//TODO rtol and atol parameters into the INPUT FILE (more user friendly)
                if (absErr < prms_ptr->atol) {
                    converged[k] = true;
                } else {
                    (std::abs(extrap[k][n]) < realmin) ? relErr = realmax : relErr = absErr / std::abs(extrap[k][n]);
                    converged[k] = relErr < prms_ptr->rtol + prms_ptr->atol / std::abs(extrap[k][n]);
                }
            }
        }//end loop over kernel

        //check if all converged then break the extrapolation loop
        allconverged = true;
        for (auto &&c : converged) {
            allconverged = allconverged && c;
        }

        if (allconverged) {
            break;
        }
		
    }//end extrapolation

    //output the E and H fields
    fields.resize(nKernels);
	
    //the results are the last element in the extrapolation vector wrt each kernel
    for (unsigned int k = 0; k < fields.size(); ++k) {
        //NOTE I multiply (1 / (2PI) ) here as I omitted it when evaluating the integrals 
		fields[k] = extrap[k].back() / (2. * prms_ptr->PI);
    }
}
//****************************************************************************************************

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Some common functions that need to compute the potential coefficients
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//compute the wave coefficients u
// for this lambda (this receiver)
// get them in all the layers in the model
void GET_CSEM1D_FD_QWE::u_comp(const double &lamb,
																	std::vector<std::complex<double>> &u,
																	std::vector<std::complex<double>> &expuh) const {

    //complex number i
    const std::complex<double> i(0, 1);

    for (unsigned int k = 0; k < u.size(); ++k) {
        u[k] = sqrt(lamb * lamb + i * omega * prms_ptr->mu0 * prms_ptr->sig[k]);
		expuh[k] = std::exp(-1. * u[k] * prms_ptr->h[k]);
    }
}

//function for computation of a, b (TM mode) or c, d (TE mode) potential coefs
 void GET_CSEM1D_FD_QWE::potential_coefs_comp(const bool &mode, //true for TM mode, false for TE mode
																							const std::vector<std::complex<double>> &u,
																							const std::vector<std::complex<double>> &expuh,
																							Source_amplitude &TM_TE_src, //TM or TE source - input
																							Potential_Coefs &pot_coefs) const //return a, b or c, d
{
	//compute the reflection coefficients (R) on the boundaries above and below the source
    // (p: below, m: above)
	Reflection_Coefs R; 
			
	//true for TM mode, false for TE mode
	(mode) ? TM_reflection_coefs(u, expuh, R) : TE_reflection_coefs(u, expuh, R);

    //*****************************************************************************************************************
    //Compute the potential coefficients for this mode
    //--------------------------------------------------------------------
    //the attenuation from the source to base and top of the source layer
	std::complex<double> expbase, exptop;
	double Tx_height = prms_ptr->z[TxLayer+1] - prms_ptr->Tx_list[*Tx_ptr][2];
	double Tx_depth = prms_ptr->Tx_list[*Tx_ptr][2] - prms_ptr->z[TxLayer];				//note z[TxLayer] <= Tx[2] <= z[TxLayer+1] - Tx[2]
	
	expbase = std::exp(-1. * u[TxLayer] * Tx_height);
	exptop   = std::exp(-1. * u[TxLayer] * Tx_depth);
	
	//compute a[TxLayer] or  c[TxLayer]
	//NOTE Rp[0] here is Rp[TxLayer] in the algorithm because we had set its index as above codes 
	pot_coefs.base = ( TM_TE_src.srcp * expbase + R.Rm[TxLayer] *  TM_TE_src.srcm * exptop * expuh[TxLayer]) 
									* R.Rp[0] / (1. -  R.Rp[0] * R.Rm[TxLayer] * expuh[TxLayer] * expuh[TxLayer]);
	
	//compute b[TxLayer] or  d[TxLayer]	
	pot_coefs.top = (TM_TE_src.srcm * exptop + R.Rp[0] * TM_TE_src.srcp * expbase * expuh[TxLayer]) 
								* R.Rm[TxLayer] / (1. -  R.Rp[0] * R.Rm[TxLayer] * expuh[TxLayer] * expuh[TxLayer]);
	
	//For the receiver located in the layer that is above the source layer
	if (RxLayer < TxLayer){
		std::complex<double> temp_base, temp_top;
		std::complex<double> src_m = TM_TE_src.srcm * exptop;
		for (unsigned int k = TxLayer - 1; k >= RxLayer; --k){
			
			if (k == invalid_unsigned_int) break; // to prevent the -1 value that is set equal to the maximum of int
			
			temp_base = (pot_coefs.base * expuh[k+1] + pot_coefs.top + src_m) 
								/ (1. + R.Rm[k] * expuh[k]);
			temp_top = R.Rm[k] * temp_base;
			pot_coefs.base = temp_base;
			pot_coefs.top = temp_top;
			src_m = 0.; // only on the top of the TxLayer exist the source term
		}
	}
	
	//For the receiver located in the layer that is below the source layer
	//NOTE RTMp is called by its own index
	if (RxLayer > TxLayer){
		std::complex<double> temp_base, temp_top;
		std::complex<double> src_p = TM_TE_src.srcp * expbase;
		for (unsigned int k = TxLayer + 1; k <= RxLayer; ++k){
			unsigned int index =  k - TxLayer; //for call the RTMp with its indices
			temp_top = (pot_coefs.base  + pot_coefs.top * expuh[k-1] + src_p) 
							  / (1. + R.Rp[index] * expuh[k]);
			temp_base = R.Rp[index] * temp_top;
			pot_coefs.base = temp_base;
			pot_coefs.top = temp_top;
			src_p = 0.; // only on the base of the TxLayer exist the source term
		}
	}
}

void GET_CSEM1D_FD_QWE::TM_reflection_coefs(const std::vector<std::complex<double>> &u,
																						const std::vector<std::complex<double>> &expuh,
																						Reflection_Coefs &RTM) const
{
	//compute the TM mode reflection coefficients (RTM) on the boundaries above and below the source
    //initialize Rp, Rm (p: below, m: above)
	RTM.Rp.resize(prms_ptr->nLayers - TxLayer);
	RTM.Rm.resize(TxLayer + 1);
	//initialize rp, rm
	RTM.rp.resize(prms_ptr->nLayers - TxLayer);
	RTM.rm.resize(TxLayer + 1);
	
	//compute rp - NOTE rp[nLayers-1] =0 
	// on z[nLayers] is the base of the [nLayers-1]_th layer
	// so recursively compute from [nLayers-2]_th to [TxLayer] = indices of the layers
	for (unsigned int k = prms_ptr->nLayers - 2; k >= TxLayer ; --k){ 
		
		if (k == invalid_unsigned_int) break; // to prevent the -1 value that is set equal to the maximum of int
		
		unsigned int index = k - TxLayer;
		
		RTM.Rp[index +1] *= expuh[k+1];
		
		RTM.rp[index] = (prms_ptr->sig[k+1] * u[k] - prms_ptr->sig[k] * u[k+1]) 
									/ (prms_ptr->sig[k+1] * u[k] + prms_ptr->sig[k] * u[k+1]);
		
		RTM.Rp[index] =(RTM.rp[index]  + RTM.Rp[index+1] * expuh[k+1]) 
										/ (1. + RTM.rp[index]  * RTM.Rp[index+1] * expuh[k+1]);	
	}
	
	//recursively compute Rm from the upper-most boundary to the top of the source layer
	//compute rm - NOTE rm[0] = 0 on z[0], similarly for Rm
	// so recursively compute from z[1] to z[TxLayer]
	for (unsigned int k = 1; k <= TxLayer; ++k){
		
		RTM.Rm[k-1] *= expuh[k-1];
		
		RTM.rm[k] = (prms_ptr->sig[k-1] * u[k] - prms_ptr->sig[k] * u[k-1]) 
							 / (prms_ptr->sig[k-1] * u[k] + prms_ptr->sig[k] * u[k-1]);
		
		RTM.Rm[k] =  (RTM.rm[k]  + RTM.Rm[k-1] * expuh[k-1]) 
								/ (1. + RTM.rm[k]  * RTM.Rm[k-1] * expuh[k-1]);
	}
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void GET_CSEM1D_FD_QWE::TE_reflection_coefs(const std::vector<std::complex<double>> &u,
																						const std::vector<std::complex<double>> &expuh,
																						Reflection_Coefs &RTE) const
{
	//compute the TE mode reflection coefficients (RTE) on the boundaries above and below the source
    //initialize Rp, Rm (p: below, m: above)
	RTE.Rp.resize(prms_ptr->nLayers - TxLayer);
	RTE.Rm.resize(TxLayer + 1);
	//initialize rp, rm
	RTE.rp.resize(prms_ptr->nLayers - TxLayer);
	RTE.rm.resize(TxLayer + 1);
	
	//compute rp - NOTE rp[nLayers-1] =0 
	// on z[nLayers] is the base of the [nLayers-1]_th layer
	// so recursively compute from [nLayers-2]_th to [TxLayer] = indices of the layers
	for (unsigned int k = prms_ptr->nLayers - 2; k >= TxLayer; --k){
		
		if (k == invalid_unsigned_int) break; // to prevent the -1 value that is set equal to the maximum of int
		
		unsigned int index = k - TxLayer;
		
		RTE.Rp[index+1] *= expuh[k+1];
		
		RTE.rp[index] = (u[k] - u[k+1]) / (u[k] + u[k+1]);
		
		RTE.Rp[index] = (RTE.rp[index]  + RTE.Rp[index+1] * expuh[k+1]) 
									 / (1. + RTE.rp[index]  * RTE.Rp[index+1] * expuh[k+1]);
	}
	
	//recursively compute Rm from the upper-most boundary to the top of the source layer
	//compute rm - NOTE rm[0] = 0 on z[0], similarly for RTEm
	// so recursively compute from z[1] to z[TxLayer]
	for (unsigned int k = 1; k <= TxLayer; ++k){
		
		RTE.Rm[k-1] *= expuh[k-1];
		
		RTE.rm[k] = (u[k] - u[k-1]) / (u[k] + u[k-1]);
		
		RTE.Rm[k] = (RTE.rm[k]  + RTE.Rm[k-1] * expuh[k-1]) 
								/ (1. + RTE.rm[k]  * RTE.Rm[k-1] * expuh[k-1]);
	}
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//compute the vector potentials and theirs derivatives w.r.t z
//NOTE for both TM and TE modes, the potential vector has only z-component
void GET_CSEM1D_FD_QWE::potential_vector(const std::complex<double> &uu,
																				 const std::vector<std::complex<double>> &u,
																				 const std::vector<double> &receiver,
																				 const Source_amplitude &src, //TM or TE source terms
																				 const Potential_Coefs &pot_coefs, //a,b (TM) or c,d (TE)
																				 Potential_Vectors &pot_vec) const  //Az or Fz and their derivatives (TM_pot_vec, TE_pot_vec)
{ 
	
	//the attenuation w.r.t z from the base and top to the receiver position in the receiver layer
	std::complex<double> expbase, exptop;
	//NOTE z[RxLayer] <= Rx[2] <= z[RxLayer+1]
	double Rx_height = prms_ptr->z[RxLayer+1] - receiver[2];
	double Rx_depth  = receiver[2] - prms_ptr->z[RxLayer];
	
	expbase = std::exp(-1. * u[RxLayer] * Rx_height);
	exptop    = std::exp(-1. * u[RxLayer] * Rx_depth);
	
	//Pz
	pot_vec.Pz = pot_coefs.base * expbase + pot_coefs.top * exptop;
	//dPzdz
	pot_vec.dPzdz = u[RxLayer]  * (pot_coefs.base * expbase - pot_coefs.top * exptop);
	//dPz2dz2
	pot_vec.dPz2dz2 = uu * pot_vec.Pz;
	
	//=================================================================================
	//add the source term if the receiver in the same layer with the source
	if (RxLayer == TxLayer){
		
		//the attenuation w.r.t z from the source to the receiver
		std::complex<double> exp_source;
		double dz = std::fabs(receiver[2] - prms_ptr->Tx_list[*Tx_ptr][2]);
		exp_source   = std::exp(-1. * u[TxLayer] * dz);
		
		//receiver is below the source
		if (receiver[2] >= prms_ptr->Tx_list[*Tx_ptr][2]){
			std::complex<double> src_p = src.srcp * exp_source;
			pot_vec.Pz += src_p;
			pot_vec.dPzdz += (-1.) * u[TxLayer] * src_p;
			pot_vec.dPz2dz2 +=  uu * src_p;
		}
		
		//receiver is above the source
		if (receiver[2] < prms_ptr->Tx_list[*Tx_ptr][2]){
			std::complex<double> src_m = src.srcm * exp_source;
			pot_vec.Pz += src_m;
			pot_vec.dPzdz += u[TxLayer] * src_m;
			pot_vec.dPz2dz2 +=  uu * src_m;
		}
		
		/*
		//receiver and source are at the same depth
		//the source term is max at this depth, thus no derivatives
		if (receiver[2] == prms_ptr->Tx_list[*Tx_ptr][2]){
			//std::complex<double> src_m = src.srcm * exp_source;
			pot_vec.Pz += src.srcm;
			pot_vec.dPzdz += 0.;
			pot_vec.dPz2dz2 +=  0.;
		}
		*/
	}
}

							
//==========================================================================================================
//This part is for x-oriented HED source (type 0)
//==========================================================================================================
//function returns the kernels for E and H fields
// with kernels is the partial sum results at interval n
// kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
// with n = 0 to nIntervals-1
//dx, dy are the vertical distance between this receiver and source
void GET_CSEM1D_FD_QWE::get_kernels_hed(const unsigned int &n,
																				  const double &dx,
																				  const double &dy,
																			      const double &rho,
																				  const std::vector<double> &receiver,
																				  std::vector<std::complex<double>> &kernels) const
 {
    //complex number i
    const std::complex<double> i(0, 1);

    //initialize kernels for x, y, z components of E and H for output, respectively
    //each component is a complex number
    kernels.resize(6);

    //set kernel = 0 before summing
    for (auto &kernel : kernels) {
        kernel = 0. + 0. * i;
    }

    //the coefficient to complete the integral evaluation
    double bma = (prms_ptr->xInt[n + 1] - prms_ptr->xInt[n]) / rho;

    unsigned int index;

    //initialize lambda value at the current index
    double lamb;
	
	//for set TM mode (true) or TE mode (false)
	bool mode;
	
	//Here we loop over the quadrature points to compute the integral
    for (unsigned int j = 0; j < prms_ptr->nQuad; ++j) {

        index = n * prms_ptr->nQuad + j;

        //get current lambda value at current index
        lamb = prms_ptr->Bx[index] / rho; 

        //get the wave coefficients u for current lambda
        std::vector<std::complex<double>> u(prms_ptr->sig.size());
		//and the exp(- u(i) * h(i) )
		std::vector<std::complex<double>> expuh(prms_ptr->sig.size());
        u_comp(lamb, u, expuh);	
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//define the amplitude of the TM mode source dipole (ap)
		Source_amplitude TM_src;
		TM_src.srcp = (-1.) * prms_ptr->moment_list[*Tx_ptr] / (2. * lamb * lamb);
		TM_src.srcm = (-1.) * TM_src.srcp; //(ap- = - ap+)

	    //compute the a,b coefficients at this lamb
        //in the receiver layer - RxLayer
		Potential_Coefs TM_coefs;
		mode = true;
        potential_coefs_comp(mode, u, expuh, TM_src, TM_coefs); // for TM mode

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//define the amplitude of the TE mode source dipole (fp)
		//fpp, fpm are for the source term below and above the source
		Source_amplitude TE_src;
		TE_src.srcp = (-1.) * i * omega * prms_ptr->mu0 
								* prms_ptr->moment_list[*Tx_ptr]
								/ (2. * lamb * lamb * u[TxLayer]);
		TE_src.srcm = TE_src.srcp;
		
        //compute the c,d coefficients at this lamb
        //in the receiver layer - RxLayer
        Potential_Coefs TE_coefs;
		mode = false;
        potential_coefs_comp(mode, u, expuh, TE_src, TE_coefs); // for TE mode

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the vector potentials and their derivatives w.r.t z-direction in the receiver layer (i_th)
		//u[RxLayer] ^2
		std::complex<double> uu = u[RxLayer]  * u[RxLayer];
		//TM mode (Az, dAzdz, dAz2dz2 is store in TM_potential struct)
        Potential_Vectors TM_potential;
		potential_vector(uu, u, receiver, TM_src, TM_coefs, TM_potential);
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//TE mode (Fz, dFzdz, dFz2dz2 is store in TE_potential struct)
		Potential_Vectors TE_potential;
		potential_vector(uu, u, receiver, TE_src, TE_coefs, TE_potential);
			
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the potentials for E field (add up TM and TE mode)
        std::complex<double> exj0, exj1, eyj0, eyj1, ezj1;
        ex_hed_j0(lamb, TM_potential, TE_potential, dx, dy, rho, exj0);
        ex_hed_j1(lamb, TM_potential, TE_potential, dx, dy, rho, exj1);
        ey_hed_j0(lamb, TM_potential, TE_potential, dx, dy, rho, eyj0);
        ey_hed_j1(lamb, TM_potential, TE_potential, dx, dy, rho, eyj1);
        ez_hed_j1(lamb, TM_potential, dx, rho, ezj1);

        //compute the potentials for H field (add up TM and TE mode)
        std::complex<double> hxj0, hxj1, hyj0, hyj1, hzj1;
        hx_hed_j0(lamb, TM_potential, TE_potential, dx, dy, rho, hxj0);
        hx_hed_j1(lamb, TM_potential, TE_potential, dx, dy, rho, hxj1);
        hy_hed_j0(lamb, TM_potential, TE_potential, dx, dy, rho, hyj0);
        hy_hed_j1(lamb, TM_potential, TE_potential, dx, dy, rho, hyj1);
        hz_hed_j1(lamb, TE_potential, dy, rho, hzj1);
		
        //now we can get compute the partial integral at the interval i by Gauss quadrature rule
        //for Ex, Ey, Ez
 		kernels[0] += exj0 * prms_ptr->J0xW[index] + exj1 * prms_ptr->J1xW[index];
		kernels[1] += eyj0 * prms_ptr->J0xW[index] + eyj1 * prms_ptr->J1xW[index];
        kernels[2] += ezj1 * prms_ptr->J1xW[index];

        //for Hx, Hy, Hz
        kernels[3] += hxj0 * prms_ptr->J0xW[index] + hxj1 * prms_ptr->J1xW[index];
        kernels[4] += hyj0 * prms_ptr->J0xW[index] + hyj1 * prms_ptr->J1xW[index];
        kernels[5] += hzj1 * prms_ptr->J1xW[index];
       }
	
    //multiply with the coefficient to complete the integral evaluation
    for (auto &kernel : kernels) {
        kernel *= bma;
    }
}

//-------------------------------------------------------------------------------------------------------------------
//compute the vector potential before implementing the quadrature sum
void GET_CSEM1D_FD_QWE::ex_hed_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &exj0) const {
	double rho2 = rho * rho;
	
    exj0 = ( -1. * dx * dx * TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dxdz/sig
	          + dy * dy * TE_pot_vec.Pz ) // -dFdy
			  * lamb * lamb * lamb / rho2 ; 
			  
}

void GET_CSEM1D_FD_QWE::ex_hed_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &exj1) const {
	
	double rho3 = rho * rho * rho;
    
	exj1 = (dx * dx - dy * dy) * lamb * lamb 
				* (TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dxdz/sig
				    + TE_pot_vec.Pz) // -dFdy
				/ rho3;
}

void GET_CSEM1D_FD_QWE::ey_hed_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &eyj0) const {
	double rho2 = rho * rho;
	
    eyj0 = (-1. * dx * dy) * lamb * lamb * lamb 
				* (TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dydz/sig
				    + TE_pot_vec.Pz) // dFdx
				/ rho2;

}

void GET_CSEM1D_FD_QWE::ey_hed_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &eyj1) const {
	double rho3 = rho * rho * rho;
	
    eyj1 = (2. * dx * dy) * lamb * lamb 
				* (TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dydz/sig
				    + TE_pot_vec.Pz) // dFdx
				/ rho3;
}

void GET_CSEM1D_FD_QWE::ez_hed_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const double &dx,
																		const double &rho,
																		std::complex<double> &ezj1) const {

	//complex number i
    const std::complex<double> i(0, 1);
	
    ezj1 = (-1.) * dx * lamb * lamb 
			* (TM_pot_vec.dPz2dz2 / prms_ptr->sig[RxLayer] //dA2dz2/sig
				-  i * omega * prms_ptr->mu0 * TM_pot_vec.Pz ) // -i.w.mu.A
			/rho;
}

//--------------------------------------------------------------------------------------------------------
void GET_CSEM1D_FD_QWE::hx_hed_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hxj0) const {

	//complex number i
    const std::complex<double> i(0, 1);
	double rho2 = rho * rho;
	
	hxj0 = (-1. * dx * dy) * lamb * lamb * lamb 
				* (TM_pot_vec.Pz // dAdy
				    + TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) // dF2dxdz/(i * w * mu0)
				/ rho2;
}

void GET_CSEM1D_FD_QWE::hx_hed_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hxj1) const {

	//complex number i
    const std::complex<double> i(0, 1);
	double rho3 = rho * rho * rho;
	
    hxj1 = (2. * dx * dy) * lamb * lamb 
				* (TM_pot_vec.Pz // dAdy
				    + TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) // dF2dxdz/(i * w * mu0)
				/ rho3;
}

void GET_CSEM1D_FD_QWE::hy_hed_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hyj0) const {
	
	//complex number i
    const std::complex<double> i(0, 1);
	double rho2 = rho * rho;
	
    hyj0 = (  dx * dx * TM_pot_vec.Pz // - dAdx
	          - dy * dy * TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) // dF2dydz/ (i * w * mu0)
			  * lamb * lamb * lamb / rho2 ; 
}

void GET_CSEM1D_FD_QWE::hy_hed_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hyj1) const {
	//complex number i
    const std::complex<double> i(0, 1);
	double rho3 = rho * rho * rho;
	
    hyj1 = (dy * dy - dx * dx) * lamb * lamb 
				* (TM_pot_vec.Pz // - dAdx
				    + TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) // dF2dydz/ (i * w * mu0)
				/ rho3;
}

void GET_CSEM1D_FD_QWE::hz_hed_j1(const double &lamb,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hzj1) const {
	//complex number i
    const std::complex<double> i(0, 1);
	
    hzj1 = (-1.) * dy * lamb * lamb 
			* (TE_pot_vec.dPz2dz2 / (i * omega * prms_ptr->mu0) //dF2dz2/(i * w * mu0)
				- prms_ptr->sig[RxLayer] * TE_pot_vec.Pz ) // -sig.F
			/rho;
}


//==============================================================================================================
//This part is for z-oriented VED source (type 1)
//==============================================================================================================
//function returns the kernels for E and H field
// with kernels is the partial sum results at interval n
// kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
// with n = 0 to nIntervals-1
void GET_CSEM1D_FD_QWE::get_kernels_ved(const unsigned int &n,
																				  const double &dx,
																				  const double &dy,
																				  const double &rho,
																				  const std::vector<double> &receiver,
																				  std::vector<std::complex<double>> &kernels)  const
{
	
    //complex number i
    const std::complex<double> i(0, 1);

    //initialize kernels for x, y, z components of E and H for output, respectively
    //each component is a complex number
    kernels.resize(6);

    //set kernel = 0 before summing
    for (auto &kernel : kernels) {
        kernel = 0. + 0. * i;
    }
	
    //the coefficient to complete the integral evaluation
    double bma = (prms_ptr->xInt[n + 1] - prms_ptr->xInt[n]) / rho;

    //initialize index to inqure the element from the vectors:
    // Bx, J0xW, J1xW
    unsigned int index;

    //initialize lambda value at the current index
    double lamb;
	
	//for set TM mode (true) or TE mode (false)
	bool mode = true;

    for (unsigned int j = 0; j < prms_ptr->nQuad; ++j) {

         index = n * prms_ptr->nQuad + j;

        //get current lambda value at current index
        lamb = prms_ptr->Bx[index] / rho; 

        //get the wave coefficients u for current lambda
        std::vector<std::complex<double>> u(prms_ptr->sig.size());
		//and the exp(- u(i) * h(i) )
		std::vector<std::complex<double>> expuh(prms_ptr->sig.size());
        u_comp(lamb, u, expuh);
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//define the amplitude of the TM mode source dipole (ap)
		Source_amplitude TM_src;
		TM_src.srcp = prms_ptr->moment_list[*Tx_ptr] / (2. * u[TxLayer]);  // m/(2u_j)
		TM_src.srcm = TM_src.srcp; //(ap- = ap+)
		
	    //compute the a,b coefficients at this lamb
        //in the receiver layer - RxLayer
		Potential_Coefs TM_coefs;
        potential_coefs_comp(mode, u, expuh, TM_src, TM_coefs); // for TM mode

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the vector potentials and their derivatives w.r.t z-direction in the receiver layer (i_th)
		//u[RxLayer] ^2
		std::complex<double> uu = u[RxLayer]  * u[RxLayer];
		//TM mode (Az, dAzdz, dAz2dz2 is store in TM_potential struct)
        Potential_Vectors TM_potential;
		potential_vector(uu, u, receiver, TM_src, TM_coefs, TM_potential);

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the potentials for E field
        std::complex<double> exj1, eyj1, ezj0;
        ex_ved_j1(lamb, TM_potential, dx, rho, exj1);
        ey_ved_j1(lamb, TM_potential, dy, rho, eyj1);
		ez_ved_j0(lamb, TM_potential, ezj0);

        //compute the potentials for H field
        std::complex<double> hxj1, hyj1;
        hx_ved_j1(lamb, TM_potential, dy, rho, hxj1);
        hy_ved_j1(lamb, TM_potential, dx, rho,  hyj1);

        //now we can get compute the partial integral at the interval i by quadrature rule
        //for Ex, Ey, Ez
        kernels[0] += exj1 * prms_ptr->J1xW[index];
        kernels[1] += eyj1 * prms_ptr->J1xW[index];
		kernels[2] += ezj0 * prms_ptr->J0xW[index];
        //Hz = 0;

        //for Hx, Hy, Hz
        kernels[3] += hxj1 * prms_ptr->J1xW[index];
        kernels[4] += hyj1 * prms_ptr->J1xW[index];
        //Hz = 0;
    }

    //multiply with the coefficient to complete the integral evaluation
    for (auto &kernel : kernels) {
        kernel *= bma;
    }
}

//------------------------------------------------------------------------------
//compute the vector potential before implementing the quadrature sum

void GET_CSEM1D_FD_QWE::ex_ved_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const double &dx,
																		const double &rho,
																		std::complex<double> &exj1) const {

    exj1 = (-1.) * dx * lamb * lamb * TM_pot_vec.dPzdz
			  / (rho * prms_ptr->sig[RxLayer]);
}

void GET_CSEM1D_FD_QWE::ey_ved_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &eyj1) const {

    eyj1 = (-1.) * dy * lamb * lamb * TM_pot_vec.dPzdz
			  / (rho * prms_ptr->sig[RxLayer]);
}

void GET_CSEM1D_FD_QWE::ez_ved_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		std::complex<double> &ezj0) const {
    //complex number i
    const std::complex<double> i(0, 1);

    ezj0 = (TM_pot_vec.dPz2dz2 / prms_ptr->sig[RxLayer] - i * omega * prms_ptr->mu0 * TM_pot_vec.Pz) * lamb;
}

//--------------------------------------------------------------------------
void GET_CSEM1D_FD_QWE::hx_ved_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hxj1) const {

    hxj1 = (-1.) * dy * lamb * lamb * TM_pot_vec.Pz
			  / rho;
}

void GET_CSEM1D_FD_QWE::hy_ved_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const double &dx,
																		const double &rho,
																		std::complex<double> &hyj1) const {

    hyj1 = dx * lamb * lamb * TM_pot_vec.Pz
			  / rho;
}


//==========================================================================================================
//This part is for x-oriented HMD source (type 2)
//==========================================================================================================
//function returns the kernels for E and H field
// with kernels is the partial sum results at interval n
// kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
// with n = 0 to nIntervals-1
// RxLayer is index of the layer containing this receiver
void GET_CSEM1D_FD_QWE::get_kernels_hmd(const unsigned int &n,
																					const double &dx,
																					const double &dy,
																					const double &rho,
																					const std::vector<double> &receiver,
																					std::vector<std::complex<double>> &kernels)  const
{
    //complex number i
    const std::complex<double> i(0, 1);

    //initialize kernels for x, y, z components of E and H for output, respectively
    //each component is a complex number
    kernels.resize(6);

    //set kernel = 0 before summing
    for (auto &kernel : kernels) {
        kernel = 0. + 0. * i;
    }

    //the coefficient to complete the integral evaluation
    double bma = (prms_ptr->xInt[n + 1] - prms_ptr->xInt[n]) / rho;

    unsigned int index;

    //initialize lambda value at the current index
    double lamb;
	
	//for set TM mode (true) or TE mode (false)
	bool mode;

    for (unsigned int j = 0; j < prms_ptr->nQuad; ++j) {

         index = n * prms_ptr->nQuad + j;

        //get current lambda value at current index
        lamb = prms_ptr->Bx[index] / rho; 

        //get the wave coefficients u for current lambda
        std::vector<std::complex<double>> u(prms_ptr->sig.size());
		//and the exp(- u(i) * h(i) )
		std::vector<std::complex<double>> expuh(prms_ptr->sig.size());
        u_comp(lamb, u, expuh);
		
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//define the amplitude of the TM mode source dipole (ap)
		Source_amplitude TM_src;
		TM_src.srcp = i * omega * prms_ptr->mu0 
								* prms_ptr->sig[TxLayer] * prms_ptr->moment_list[*Tx_ptr] 
								/ (2. * lamb * lamb * u[TxLayer]); //i.w.mu.sig_j.m / (2.lamb^2.u_j)
		TM_src.srcm = TM_src.srcp; //(ap- = ap+)
		
	    //compute the a,b coefficients at this lamb
        //in the receiver layer - RxLayer
		Potential_Coefs TM_coefs;
		mode = true;
        potential_coefs_comp(mode, u, expuh, TM_src, TM_coefs); // for TM mode

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//define the amplitude of the TE mode source dipole (fp)
		//fpp, fpm are for the source term below and above the source
		Source_amplitude TE_src;
		TE_src.srcp = (-1.) * i * omega * prms_ptr->mu0 
								* prms_ptr->moment_list[*Tx_ptr]  
								/ (2. * lamb * lamb); //i.w.mu.sig_j.m / (2.lamb^2.u_j)
		TE_src.srcm = (-1.) * TE_src.srcp; //(fp- = - fp+)
		
        //compute the c,d coefficients at this lamb
        //in the receiver layer - RxLayer
        Potential_Coefs TE_coefs;
		mode = false;
        potential_coefs_comp(mode, u, expuh, TE_src, TE_coefs); // for TE mode

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the vector potentials and their derivatives w.r.t z-direction in the receiver layer (i_th)
		//u[RxLayer] ^2
		std::complex<double> uu = u[RxLayer]  * u[RxLayer];
		//TM mode (Az, dAzdz, dAz2dz2 is store in TM_potential struct)
        Potential_Vectors TM_potential;
		potential_vector(uu, u, receiver, TM_src, TM_coefs, TM_potential);
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//TE mode (Fz, dFzdz, dFz2dz2 is store in TE_potential struct)
		Potential_Vectors TE_potential;
		potential_vector(uu, u, receiver, TE_src, TE_coefs, TE_potential);
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the potentials for E field
        std::complex<double> exj0, exj1, eyj0, eyj1, ezj1;
        ex_hmd_j0(lamb, TM_potential, TE_potential, dx, dy, rho, exj0);
        ex_hmd_j1(lamb, TM_potential, TE_potential, dx, dy, rho, exj1);
        ey_hmd_j0(lamb, TM_potential, TE_potential, dx, dy, rho, eyj0);
        ey_hmd_j1(lamb, TM_potential, TE_potential, dx, dy, rho, eyj1);
        ez_hmd_j1(lamb, TM_potential, dy, rho, ezj1);

        //compute the potentials for H field
        std::complex<double> hxj0, hxj1, hyj0, hyj1, hzj1;
        hx_hmd_j0(lamb, TM_potential, TE_potential, dx, dy, rho, hxj0);
        hx_hmd_j1(lamb, TM_potential, TE_potential, dx, dy, rho, hxj1);
        hy_hmd_j0(lamb, TM_potential, TE_potential, dx, dy, rho, hyj0);
        hy_hmd_j1(lamb, TM_potential, TE_potential, dx, dy, rho, hyj1);
        hz_hmd_j1(lamb, TE_potential, dx, rho, hzj1);

        //now we can get compute the partial integral at the interval i by quadrature rule
        //for Ex, Ey, Ez
        kernels[0] += exj0 * prms_ptr->J0xW[index] + exj1 * prms_ptr->J1xW[index];
        kernels[1] += eyj0 * prms_ptr->J0xW[index] + eyj1 * prms_ptr->J1xW[index];
        kernels[2] += ezj1 * prms_ptr->J1xW[index];

        //for Hx, Hy, Hz
        kernels[3] += hxj0 * prms_ptr->J0xW[index] + hxj1 * prms_ptr->J1xW[index];
        kernels[4] += hyj0 * prms_ptr->J0xW[index] + hyj1 * prms_ptr->J1xW[index];
        kernels[5] += hzj1 * prms_ptr->J1xW[index];
		
    }

    //multiply with the coefficient to complete the integral evaluation
    for (auto &kernel : kernels) {
        kernel *= bma;
    }
}

//-------------------------------------------------------------------------------------------------------------------
//compute the vector potential before implementing the quadrature sum
void GET_CSEM1D_FD_QWE::ex_hmd_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &exj0) const {

    double rho2 = rho * rho;
	
    exj0 = (-1. * dx * dy) * lamb * lamb * lamb 
				* (TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dxdz/sig
				    - TE_pot_vec.Pz) // -dFdy
				/ rho2;
}

void GET_CSEM1D_FD_QWE::ex_hmd_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &exj1) const {
									  
    double rho3 = rho * rho * rho;

    exj1 =  (2. * dx * dy) * lamb * lamb 
				* (TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dxdz/sig
				    - TE_pot_vec.Pz) //  -dFdy
				/ rho3;
}

void GET_CSEM1D_FD_QWE::ey_hmd_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &eyj0) const {
									  
	 double rho2 = rho * rho;

    eyj0 = ( -1. * dy * dy * TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dydz/sig
	          - dx * dx * TE_pot_vec.Pz ) // dFdx
			  * lamb * lamb * lamb / rho2 ; 

}

void GET_CSEM1D_FD_QWE::ey_hmd_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &eyj1) const {
	
	double rho3 = rho * rho * rho;
	
    eyj1 = (dy * dy - dx * dx) * lamb * lamb 
				* (TM_pot_vec.dPzdz / prms_ptr->sig[RxLayer] // dA2dydz/sig
				    - TE_pot_vec.Pz) // dFdx
				/ rho3;
	
}

void GET_CSEM1D_FD_QWE::ez_hmd_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &ezj1) const {

   	//complex number i
    const std::complex<double> i(0, 1);
	
    ezj1 = (-1.) * dy * lamb * lamb 
			* (TM_pot_vec.dPz2dz2 / prms_ptr->sig[RxLayer] //dA2dz2/sig
				-  i * omega * prms_ptr->mu0 * TM_pot_vec.Pz ) // -i.w.mu.A
			/rho;
}

//--------------------------------------------------------------------------------------------------------
void GET_CSEM1D_FD_QWE::hx_hmd_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hxj0) const {
									  
	 //complex number i
    const std::complex<double> i(0, 1);
	double rho2 = rho * rho;

    hxj0 = (  -1. * dy * dy * TM_pot_vec.Pz //  dAdy
	          - dx * dx * TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) // dF2dxdz/ (i * w * mu0)
			  * lamb * lamb * lamb / rho2 ; 
}

void GET_CSEM1D_FD_QWE::hx_hmd_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hxj1) const {
	 //complex number i
	 const std::complex<double> i(0, 1);
	 double rho3 = rho * rho * rho;
	 
    hxj1 = (dy * dy - dx * dx) * lamb * lamb 
				* (TM_pot_vec.Pz // dAdy
				    - TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) //  dF2dxdz/ (i * w * mu0)
				/ rho3;
}

void GET_CSEM1D_FD_QWE::hy_hmd_j0(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hyj0) const {

     //complex number i
    const std::complex<double> i(0, 1);
	double rho2 = rho * rho;
	
	hyj0 = (-1. * dx * dy) * lamb * lamb * lamb 
				* (-1. * TM_pot_vec.Pz // -dAdx
				    + TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) // dF2dydz/(i * w * mu0)
				/ rho2;
}

void GET_CSEM1D_FD_QWE::hy_hmd_j1(const double &lamb,
																		const Potential_Vectors &TM_pot_vec,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hyj1) const {
	//complex number i
	 const std::complex<double> i(0, 1);
	 double rho3 = rho * rho * rho;
   
	hyj1 = (2. * dx * dy) * lamb * lamb 
				* (-1. * TM_pot_vec.Pz // -dAdx
				    + TE_pot_vec.dPzdz / (i * omega * prms_ptr->mu0) ) //  dF2dydz/(i * w * mu0)
				/ rho3;
}

void GET_CSEM1D_FD_QWE::hz_hmd_j1(const double &lamb,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &rho,
																		std::complex<double> &hzj1) const {
	
	//complex number i
    const std::complex<double> i(0, 1);
	
    hzj1 = (-1.) * dx * lamb * lamb 
			* (TE_pot_vec.dPz2dz2 / (i * omega * prms_ptr->mu0) //dF2dz2/(i * w * mu0)
				- prms_ptr->sig[RxLayer] * TE_pot_vec.Pz ) // -sig.F
			/rho;
}


//==============================================================================================================
//This part is for z-oriented VMD source (type 3)
//==============================================================================================================
//function returns the kernels for E and H field
// with kernels is the partial sum results at interval n
// kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
// with n = 0 to nIntervals-1
void GET_CSEM1D_FD_QWE::get_kernels_vmd(const unsigned int &n,
																					const double &dx,
																					const double &dy,
																					const double &rho,
																					const std::vector<double> &receiver,
																					std::vector<std::complex<double>> &kernels) const 
{

    //complex number i
    const std::complex<double> i(0, 1);

    //initialize kernels for x, y, z components of E and H for output, respectively
    //each component is a complex number
    kernels.resize(6);

    //set kernel = 0 before summing
    for (auto &kernel : kernels) {
        kernel = 0. + 0. * i;
    }

    //the coefficient to complete the integral evaluation
    double bma = (prms_ptr->xInt[n + 1] - prms_ptr->xInt[n]) / rho;

    //initialize index to inqure the element from the vectors:
    // Bx, J0xW, J1xW
    unsigned int index;

    //initialize lambda value at the current index
    double lamb;
	
	//for set TM mode (true) or TE mode (false)
	bool mode = false;

    for (unsigned int j = 0; j < prms_ptr->nQuad; ++j) {

         index = n * prms_ptr->nQuad + j;

        //get current lambda value at current index
        lamb = prms_ptr->Bx[index] / rho; 

        //get the wave coefficients u for current lambda
        std::vector<std::complex<double>> u(prms_ptr->sig.size());
		//and the exp(- u(i) * h(i) )
		std::vector<std::complex<double>> expuh(prms_ptr->sig.size());
        u_comp(lamb, u, expuh);
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//define the amplitude of the TE mode source dipole (fp)
		Source_amplitude TE_src;
		TE_src.srcp =  i * omega * prms_ptr->mu0 
								* prms_ptr->moment_list[*Tx_ptr] 
								/ (2. * u[TxLayer]);  // i.w.mu.m/(2u_j)
		TE_src.srcm = TE_src.srcp; //(ap- = ap+)
		
	    //compute the c, d coefficients at this lamb
        //in the receiver layer - RxLayer
		Potential_Coefs TE_coefs;
        potential_coefs_comp(mode, u, expuh, TE_src, TE_coefs); // for TE mode

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the vector potentials and their derivatives w.r.t z-direction in the receiver layer (i_th)
		//u[RxLayer] ^2
		std::complex<double> uu = u[RxLayer]  * u[RxLayer];
		//TM mode (Az, dAzdz, dAz2dz2 is store in TM_potential struct)
        Potential_Vectors TE_potential;
		potential_vector(uu, u, receiver, TE_src, TE_coefs, TE_potential);

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //compute the potentials for E field
        std::complex<double> exj1, eyj1;
        ex_vmd_j1(lamb, TE_potential, dy, rho, exj1);
        ey_vmd_j1(lamb, TE_potential, dx, rho, eyj1);

        //compute the potentials for H field
        std::complex<double> hxj1, hyj1, hzj0;
        hx_vmd_j1(lamb, TE_potential, dx, rho, hxj1);
        hy_vmd_j1(lamb, TE_potential, dy, rho, hyj1);
        hz_vmd_j0(lamb, TE_potential, hzj0);

        //now we can get compute the partial integral at the interval i by quadrature rule
        //for Ex, Ey, Ez
        kernels[0] += exj1 * prms_ptr->J1xW[index];
        kernels[1] += eyj1 * prms_ptr->J1xW[index];
        //Ez = 0;

        //for Hx, Hy, Hz
        kernels[3] += hxj1 * prms_ptr->J1xW[index];
        kernels[4] += hyj1 * prms_ptr->J1xW[index];
        kernels[5] += hzj0 * prms_ptr->J0xW[index];
    }

    //multiply with the coefficient to complete the integral evaluation
    for (auto &kernel : kernels) {
        kernel *= bma;
    }
}


//------------------------------------------------------------------------------
//compute the vector potential before implementing the quadrature sum

void GET_CSEM1D_FD_QWE::ex_vmd_j1(const double &lamb,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &exj1) const {

    exj1 = dy * lamb * lamb * TE_pot_vec.Pz
			  / rho;
}

void GET_CSEM1D_FD_QWE::ey_vmd_j1(const double &lamb,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &rho,
																		std::complex<double> &eyj1) const {

    eyj1 = (-1.) * dx * lamb * lamb * TE_pot_vec.Pz
			  / rho;
}

//--------------------------------------------------------------------------
void GET_CSEM1D_FD_QWE::hx_vmd_j1(const double &lamb,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dx,
																		const double &rho,
																		std::complex<double> &hxj1) const {
	 //complex number i
    const std::complex<double> i(0, 1);
	
    hxj1 = (-1.) * dx * lamb * lamb * TE_pot_vec.dPzdz
			  / (i * omega * prms_ptr->mu0  * rho);
}

void GET_CSEM1D_FD_QWE::hy_vmd_j1(const double &lamb,
																		const Potential_Vectors &TE_pot_vec,
																		const double &dy,
																		const double &rho,
																		std::complex<double> &hyj1) const {
	//complex number i
    const std::complex<double> i(0, 1);
	
    hyj1 = (-1.) * dy * lamb * lamb * TE_pot_vec.dPzdz
			  / (i * omega * prms_ptr->mu0  * rho);
}

void GET_CSEM1D_FD_QWE::hz_vmd_j0(const double &lamb,
																		const Potential_Vectors &TE_pot_vec,
																		std::complex<double> &hzj0) const {
	//complex number i
    const std::complex<double> i(0, 1);
	
    hzj0 = lamb * (TE_pot_vec.dPz2dz2 / (i * omega * prms_ptr->mu0 ) - TE_pot_vec.Pz * prms_ptr->sig[RxLayer]);
}

//test print_vec
void GET_CSEM1D_FD_QWE::print_vec(const auto &vec) const
{
	std::for_each(vec.begin(), vec.end(), [](auto i){std::cout << i << "\n";});
}
