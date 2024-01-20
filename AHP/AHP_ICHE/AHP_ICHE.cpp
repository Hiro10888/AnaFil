//AHP_ICHE.cpp - ver. 1.0
// Separate AnalogFiltet.h to Analoghilter.h and AnalogFilter.cpp
// last change was made on 20230802 by udata

#define MMAX 50			// Maximum order of POLE,ZERO

#include <string>
#include <complex>
#include "AnalogFilter.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

int main()
{
	const 	std::string b_title	= "Bode plot of Inverse Chebyshev HPF";
	double	f_start 			= 0.1;
	double	f_end   			= 100.0;
	double	m_min				= -200.0;
	double	m_max				= 50.0;
	double	t_start 			= 0.0;
	double	t_end   			= 5.0;

	const 	std::string i_title	= "Impulse response of Inverse Chebyshev HPF";
	double	i_min				= -50.0;
	double	i_max				= 50.0;

	const 	std::string s_title	= "Step response of Inverse Chebyshev HPF";
	double	s_min				= -1.0;
	double	s_max				= 1.0;

	const 	std::string r_title	= "Poles and Zeros  of Inverse Chebyshev HPF";
	double	r_yscale    		= 20.0;

	int							i, impn, stpn;
	std::complex<double> 		pole[MMAX], zero[MMAX], ampi[MMAX], amps[MMAX];

    AnalogFilter AHP_ICHE;
												// Data members
												// For  Filter Type selection.
		AHP_ICHE.isf	= 3;            		// 1=LPF, 2=BPF,  3=HPF, 4=BEF
		AHP_ICHE.is		= 3;					// 1=butterworth, 2=Chebyshev 3=Inv.Cheby.  4=Elliptic

												// For LPF design.
		AHP_ICHE.wp		= 10.0;		 			// Pass band frequency [rad/sec]
		AHP_ICHE.attwp 	= 5.0;					// Attenuation at wpp  [dB]
		AHP_ICHE.ws		= 7.0;		 			// Stop band frequency [rad/sec]
		AHP_ICHE.attws	= 80.0;  				// Attenuation at wss  [dB]

		// Member Functions.
	   	AHP_ICHE.hpf_inv_chebyshev(pole,zero);
		AHP_ICHE.bode_plot(pole, zero, f_start, f_end, m_min, m_max, b_title);
		AHP_ICHE.func_trans(pole, zero, ampi, &impn, amps, &stpn);
		AHP_ICHE.time_response(impn, i_min, i_max, pole, ampi, t_start, t_end, i_min, i_max, i_title);
		AHP_ICHE.time_response(stpn, s_min, s_max, pole, amps, t_start, t_end, s_min, s_max, s_title);
		AHP_ICHE.root_plot(r_yscale, pole, zero, r_title);
		AHP_ICHE.write_AnalogFilter(impn, stpn, pole, zero);
}