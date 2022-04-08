//
// Class to get EM fields for a given Transmitter at a certain frequency value
//
#ifndef GET_CSEM1D_FD_QWE_H
#define GET_CSEM1D_FD_QWE_H

#include <vector>
#include <complex>
#include <algorithm>    // std::for_each
#include <cmath>
#include <exception>
#include <iostream>

#include <boost/math/special_functions/bessel.hpp>  //for bessel function
#include <boost/math/quadrature/gauss_kronrod.hpp> //for adaptive quadrature rule

#include "Prms.h"

class GET_CSEM1D_FD_QWE {

public:

    GET_CSEM1D_FD_QWE(PRMS &prms,
												 unsigned int & Tx_index,
												 unsigned int & freq_index);

    virtual ~GET_CSEM1D_FD_QWE();
	
	//function for solving 1D problem for each Tx at each frequency
    void solve_1D_CSEM(std::vector<std::vector<std::complex<double>>> &EMfields); //output

private:

	//pointer to prms object
	PRMS *prms_ptr;
	//pointers to current Tx_index and freq_index
	unsigned int *Tx_ptr, *freq_ptr;

	//this prms change w.r.t frequency so I compute it here
	double omega;

	//check the source_type we need to compute
	bool hed, ved, hmd, vmd;

    //index of the source layer (index range:0 to nLayers-1)
    unsigned int TxLayer;

    //index of the the receiver layer (index range:0 to nLayers-1)
    unsigned int RxLayer;
	
	//singularity when rho goes to zeros
	bool singularity;
	
	static const unsigned int invalid_unsigned_int = static_cast<unsigned int>(-1);

	//object to store the source amplitude
	//can be used for either TM or TE mode
	struct Source_amplitude{
		std::complex<double> srcp; //for z >= zs - down-going source term
		std::complex<double> srcm; //for z < zs - up-going source term
	};
	
	//object to store the reflection coefficients
	//can be used for either TM or TE mode
	struct Reflection_Coefs{
		  std::vector<std::complex<double>> Rp, Rm, rp, rm;
	};
	
	//object to store potential coefficients
	//can be used for either TM or TE mode
	struct Potential_Coefs{
		std::complex<double> base; //for a or c
		std::complex<double> top; //for c or d
	};
	
	//object to store the potential vectors and theirs derivatives w.r.t. z
	//can be used for either TM or TE mode
	struct Potential_Vectors{
		 std::complex<double> Pz, dPzdz, dPz2dz2;
	};
	
	//==================================================================

    //###############################################################################################
    //setup the parameters
    void setup_1d_CSEM();

	//function to solve for Electric dipole source
	void Solve_ED(double &dx,
						  const double &dy,
						  const double &rho,
						  const std::vector<double> &receiver,
						  std::vector<std::vector<std::complex<double>>> &EMfields) const;
							
	//function to solve for Electric dipole source
	void Solve_MD(double &dx,
						   const double &dy,
						   const double &rho,
						   const std::vector<double> &receiver,
						   std::vector<std::vector<std::complex<double>>> &EMfields) const;
						   
	//compute the dx_hat, dy_hat in rotated coords when the azimuthTx != 0
	void rotate_coords_comp(const double &dx,
													const double &dy,
													const double &rho,
													double &dx_hat,
													double &dy_hat)const;
							
	// azimuthTx to get the fields in (Oxyz) because 
	// we compute the fileds with the x'-direction is the direction of the dipole source
	// only need for the horizontal component of the dipole source
	void horizontal_rotate(std::vector<std::complex<double>> &fields) const;
	
	// add the fileds of horizontal and vertical components of the dipole source
	void add_horizontal_vertical_souce(const std::vector<std::complex<double>> &horizontal_fields,
															  const std::vector<std::complex<double>> &vertical_fields,
															  std::vector<std::complex<double>> &fields) const;


    //###############################################################################################
    //each receiver computing part
    //###############################################################################################

    //*************************************************************************************************
    //function implements the Quadrature with Extrapolation method (Key,2012 and Weniger, 2003 -pp26)
    //returns the fields with field[0:2] for Ex, Ey, Ez; and fields[3:5] for Hx, Hy, Hz
    // at current receiver position
    void qwe(const unsigned int &source_type,
					const double &dx,
					const double &dy,
					const double &rho,
					const std::vector<double> &receiver,
					std::vector<std::complex<double>> &fields) const;


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Some common functions that need to compute the potential coefficients
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //compute the wave parameter u's for the layer below the source layer (uM)
    // for the layer above the source layer (uN)
    // and the value of u at the receiver layer (uRx)
    void u_comp(const double &lamb,
						  std::vector<std::complex<double>> &u,
						  std::vector<std::complex<double>> &expuh) const;

	//function for computation of a, b (TM mode) or c, d (TE mode) potential coefs
	void potential_coefs_comp(const bool &mode, //true for TM mode, false for TE mode
								  const std::vector<std::complex<double>> &u,
								  const std::vector<std::complex<double>> &expuh,
								  Source_amplitude &TM_TE_src, //TM or TE source - input
								  Potential_Coefs &pot_coefs) const; //return a, b or c, d
							
	//---------------------------------------------------------------------------------------------------
	void TM_reflection_coefs(const std::vector<std::complex<double>> &u,
												const std::vector<std::complex<double>> &expuh,
												Reflection_Coefs &RTM) const;
	
	void TE_reflection_coefs(const std::vector<std::complex<double>> &u,
												const std::vector<std::complex<double>> &expuh,
												Reflection_Coefs &RTE) const;	
												
	//----------------------------------------------------------------------------------------------------------------------------------------------------------
	//compute the vector potentials and theirs derivatives w.r.t z
	void potential_vector(const std::complex<double> &uu,
							       const std::vector<std::complex<double>> &u,
								   const std::vector<double> &receiver,
							       const Source_amplitude &src, //TM or TE source terms
							       const Potential_Coefs &pot_coefs, //a,b (TM) or c,d (TE)
							       Potential_Vectors &pot_vec) const;  //Az or Fz and their derivatives

    //===============================================================================
    //This part is for x-oriented HED source (type 0)
    //===============================================================================
    //function returns the kernels for E and H field by the QWE 
    // with kernels is the partial sum results at interval n
    // kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
    // with n = 0 to nIntervals-1
    void get_kernels_hed(const unsigned int &n,
									const double &dx,
									const double &dy,
									const double &rho,
									const std::vector<double> &receiver,
									std::vector<std::complex<double>> &kernels) const;

    //------------------------------------------------------------------------------
    //compute the vector potential before implementing the quadrature sum
    void ex_hed_j0(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &exj0) const;

    void ex_hed_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &exj1) const;

    void ey_hed_j0(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &eyj0) const;

    void ey_hed_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &eyj1) const;

    void ez_hed_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const double &dx,
						   const double &rho,
						   std::complex<double> &ezj1) const;

    //----------------------------------------------------------------------------
    void hx_hed_j0(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &hxj0) const;

    void hx_hed_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &hxj1) const;

    void hy_hed_j0(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &hyj0) const;

    void hy_hed_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dx,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &hyj1) const;

    void hz_hed_j1(const double &lamb,
						   const Potential_Vectors &TE_pot_vec,
						   const double &dy,
						   const double &rho,
						   std::complex<double> &hzj1) const;
	
	//===============================================================================
    //This part is z-oriented for VED source (type 1)
    //===============================================================================
    //function returns the kernels for E and H field
    // with kernels is the partial sum results at interval n
    // kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
    // with n = 0 to nIntervals-1
    void get_kernels_ved(const unsigned int &n,
                                    const double &dx,
                                    const double &dy,
                                    const double &rho,
                                    const std::vector<double> &receiver,
                                    std::vector<std::complex<double>> &kernels) const;

    //------------------------------------------------------------------------------
    //compute the vector potential before implementing the quadrature sum
    void ex_ved_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const double &dx,
						   const double &rho,
                           std::complex<double> &exj1) const;

    void ey_ved_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const double &dy,
						   const double &rho,
                           std::complex<double> &eyj1) const ;
							 
	 void  ez_ved_j0(const double &lamb,
						     const Potential_Vectors &TM_pot_vec,
							 std::complex<double> &ezj0) const;

    //----------------------------------------------------------------------------
    void hx_ved_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const double &dy,
						   const double &rho,
                           std::complex<double> &hxj1) const ;

    void hy_ved_j1(const double &lamb,
						   const Potential_Vectors &TM_pot_vec,
						   const double &dx,
						   const double &rho,
                           std::complex<double> &hyj1) const;
    //===============================================================================
						   
	
	//===============================================================================
    //This part is for x-oriented HMD source (type 2)
    //===============================================================================
    //function returns the kernels for E and H field
    // with kernels is the partial sum results at interval n
    // kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
    // with n = 0 to nIntervals-1
    void get_kernels_hmd(const unsigned int &n,
                                     const double &dx,
                                     const double &dy,
                                     const double &rho,
                                     const std::vector<double> &receiver,
                                     std::vector<std::complex<double>> &kernels) const;

    //------------------------------------------------------------------------------
    //compute the vector potential before implementing the quadrature sum
    void ex_hmd_j0(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &exj0) const;

	void ex_hmd_j1(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &exj1) const;
								  
	void ey_hmd_j0(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &eyj0) const;
								  
	void ey_hmd_j1(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &eyj1) const;
								  
	void ez_hmd_j1(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &ezj1) const;

    //----------------------------------------------------------------------------
    void hx_hmd_j0(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &hxj0) const ;

    void hx_hmd_j1(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &hxj1) const ;
								  
	void hy_hmd_j0(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &hyj0) const;
								  
	void hy_hmd_j1(const double &lamb,
						    const Potential_Vectors &TM_pot_vec,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &dy,
						    const double &rho,
                            std::complex<double> &hyj1) const ;
								  
	void hz_hmd_j1(const double &lamb,
						    const Potential_Vectors &TE_pot_vec,
						    const double &dx,
						    const double &rho,
                            std::complex<double> &hzj1) const;
    //===============================================================================
	
	
	//===============================================================================
    //This part is for z-oriented VMD source (type 3)
    //===============================================================================
    //function returns the kernels for E and H field
    // with kernels is the partial sum results at interval n
    // kernels[0:2] are for Ex,Ey,Ez; and kernel[3:5] are for Hx, Hy, Hz
    // with n = 0 to nIntervals-1
    void get_kernels_vmd(const unsigned int &n,
									 const double &dx,
									 const double &dy,
									 const double &rho,
									 const std::vector<double> &receiver,
									 std::vector<std::complex<double>> &kernels) const;

    //------------------------------------------------------------------------------
    //compute the vector potential before implementing the quadrature sum
    void ex_vmd_j1(const double &lamb,
							const Potential_Vectors &TE_pot_vec,
							const double &dy,
							const double &rho,
                            std::complex<double> &exj1) const;

    void ey_vmd_j1(const double &lamb,
							const Potential_Vectors &TE_pot_vec,
							const double &dx,
							const double &rho,
							std::complex<double> &eyj1) const;

    //----------------------------------------------------------------------------
    void hx_vmd_j1(const double &lamb,
							const Potential_Vectors &TE_pot_vec,
							const double &dx,
							const double &rho,
							std::complex<double> &hxj1) const;

    void hy_vmd_j1(const double &lamb,
							const Potential_Vectors &TE_pot_vec,
							const double &dy,
							const double &rho,
							std::complex<double> &hyj1) const;

    void hz_vmd_j0(const double &lamb,
							const Potential_Vectors &TE_pot_vec,
							std::complex<double> &hzj0) const;
    //===============================================================================
	
	void print_vec(const auto &vec) const;
	
}; // END OF GET_CSEM1D_FD_QWE

#endif //GET_CSEM1D_FD_QWE_H