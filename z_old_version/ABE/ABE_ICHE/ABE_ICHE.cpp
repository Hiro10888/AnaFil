//ABE_ICHE.cpp - ver. 2.0
// Separate AnalogFiltet.h to Analoghilter.h and AnalogFilter.cpp
// last change was made on 20230803 by udata

#define MMAX 50			// Maximum order of POLE,ZERO

#include <string>
#include <complex>
#include "AnalogFilter.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

int main()
{
	const 	std::string b_title	= "Bode plot of Inverse Chebyshev BEF";
	double	f_start 			= 0.1;
	double	f_end   			= 100.0;
	double	m_min				= -200.0;
	double	m_max				= 50.0;
	double	t_start 			= 0.0;
	double	t_end   			= 10.0;

	const 	std::string i_title	= "Impulse response of Inverse Chebyshev BEF";
	double	i_min				= -10.0;
	double	i_max				= 10.0;

	const 	std::string s_title	= "Step response of Inverse Chebyshev BEF";
	double	s_min				= -2.0;
	double	s_max				= 2.0;

	const 	std::string r_title	= "Poles and Zeros  of Inverse Chebyshev BEF";
	double	r_yscale    		= 40.0;

	int							impn, stpn;
	std::complex<double> 		pole[MMAX], zero[MMAX], ampi[MMAX], amps[MMAX];

    AnalogFilter ABE_ICHE;
												// Data members
												// For  Filter Type selection.
		ABE_ICHE.isf	= 4;            		// 1=LPF, 2=BPF,  3=HPF, 4=BEF
		ABE_ICHE.is		= 3;					// 1=butterworth, 2=Chebyshev 3=Inv.Cheby.  4=Elliptic

												// For BEF design.
		ABE_ICHE.ws1	= 10.0;					// Define Attenuation of lowest stopband freq.
		ABE_ICHE.ws2	= 12.0;					// Define Attenyation of highest stopband freq.
		ABE_ICHE.attws	= 80.0;  				// Attenuation at wss	[dB]
		ABE_ICHE.wp1	=  6.0;					// Lowest frequency of Pass band
		ABE_ICHE.wp2	= 20.0;					// Highest frequency of Pass band
		ABE_ICHE.attwp	=  2.0;					// Attenuation at wp	[dB]

		// Member Functions.
	   	ABE_ICHE.bef_inv_chebyshev(pole,zero);
		ABE_ICHE.bode_plot(pole, zero, f_start, f_end, m_min, m_max, b_title);
		ABE_ICHE.func_trans(pole, zero, ampi, &impn, amps, &stpn);
		ABE_ICHE.time_response(impn, i_min, i_max, pole, ampi, t_start, t_end, i_min, i_max, i_title);
		ABE_ICHE.time_response(stpn, s_min, s_max, pole, amps, t_start, t_end, s_min, s_max, s_title);
		ABE_ICHE.root_plot(r_yscale, pole, zero, r_title);
		ABE_ICHE.write_AnalogFilter(impn, stpn, pole, zero);
}