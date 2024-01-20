//ABP_ICHE.cpp - ver. 1.0
// Separate AnalogFiltet.h to Analoghilter.h and AnalogFilter.cpp
// last change was made on 20230801 by udata

#define MMAX 50			// Maximum order of POLE,ZERO

#include <string>
#include <complex>
#include "AnalogFilter.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

int main()
{
	const 	std::string b_title	= "Bode plot of Inverse Chebyshev BPF";
	double	f_start 			= 0.1;
	double	f_end   			= 100.0;
	double	m_min				= -200.0;
	double	m_max				= 50.0;
	double	t_start 			= 0.0;
	double	t_end   			= 15.0;

	const 	std::string i_title	= "Impulse response of Inverse Chebyshev BPF";
	double	i_min				= -1.0;
	double	i_max				= 1.0;

	const 	std::string s_title	= "Step response of Inverse Chebyshev BPF";
	double	s_min				= -0.25;
	double	s_max				= 0.25;

	const 	std::string r_title	= "Poles and Zeros  of Inverse Chebyshev BPF";
	double	r_yscale    		= 15.0;

	int							i, impn, stpn;
	std::complex<double> 		pole[MMAX], zero[MMAX], ampi[MMAX], amps[MMAX];

    AnalogFilter ABP_ICHE;
												// Data members
												// For  Filter Type selection.
		ABP_ICHE.isf	= 2;            		// 1=LPF, 2=BPF,  3=HPF, 4=BEF
		ABP_ICHE.is		= 3;					// 1=butterworth, 2=Chebyshev 3=Inv.Cheby.  4=Elliptic

												// For BPF design.
		ABP_ICHE.ws1	=  4.0;					// Define Attenuation of lowest stopband freq.
		ABP_ICHE.ws2	= 20.0;					// Define Attenyation of highest stopband freq.
		ABP_ICHE.attws	= 80.0;  				// Attenuation at wss	[dB]
		ABP_ICHE.wp1	=  8.0;					// Lowest frequency of Pass band
		ABP_ICHE.wp2	= 10.0;					// Highest frequency of Pass band
		ABP_ICHE.attwp	=  1.0;					// Attenuation at wp	[dB]

		// Member Functions.
	   	ABP_ICHE.bpf_inv_chebyshev(pole,zero);
		ABP_ICHE.bode_plot(pole, zero, f_start, f_end, m_min, m_max, b_title);
		ABP_ICHE.func_trans(pole, zero, ampi, &impn, amps, &stpn);
		ABP_ICHE.time_response(impn, i_min, i_max, pole, ampi, t_start, t_end, i_min, i_max, i_title);
		ABP_ICHE.time_response(stpn, s_min, s_max, pole, amps, t_start, t_end, s_min, s_max, s_title);
		ABP_ICHE.root_plot(r_yscale, pole, zero, r_title);
		ABP_ICHE.write_AnalogFilter(impn, stpn, pole, zero);
}