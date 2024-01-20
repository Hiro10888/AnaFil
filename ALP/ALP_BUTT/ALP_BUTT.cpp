// ALP_BUTT.cpp - ver. 1.3
// Separate AnalogFiltet.h to Analoghilter.h and AnalogFilter.cpp
// last change was made on 20230726 by udata

#define MMAX 50			// Maximum order of POLE,ZERO

#include <string>
#include <complex>
#include "AnalogFilter.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

int main()
{
	const 	std::string b_title 	= "Bode plot of Butterworth LPF";
	double	f_start 		= 0.1;
	double	f_end   		= 10.0;
	double	m_min			= -200.0;
	double	m_max			= 50.0;
	double	t_start 		= 0.0;
	double	t_end   		= 10.0;

	const 	std::string i_title	= "Impulse response of Butterworth LPF";
	double	i_min			= -3.0;
	double	i_max			= 3.0;

	const 	std::string s_title	= "Step response of Butterworth LPF";
	double	s_min			= -0.5;
	double	s_max			= 1.5;

	const 	std::string r_title	= "Poles and Zeros  of Butterworth LPF";
	double	r_yscale    	= 10.0;

	int						i, impn, stpn;
	std::complex<double> 		pole[MMAX], zero[MMAX], ampi[MMAX], amps[MMAX];

    AnalogFilter ALP_BUTT;

		ALP_BUTT.isf	= 1;            		// 1=LPF, 2=BPF,  3=HPF, 4=BEF
		ALP_BUTT.is		= 1;					// 1=butterworth, 2=Chebyshev
												// 3=Inv. Cheby.  4=Elliptic
		ALP_BUTT.wp		= 1.0*2.0*M_PI;		 	// Pass band frequency [rad/sec]
		ALP_BUTT.attwp 	= 1.0;					// Attenuation at wpp  [dB]
		ALP_BUTT.ws		= 2.0*2.0*M_PI;		 	// Stop band frequency [rad/sec]
		ALP_BUTT.attws	= 60.0;  				// Attenuation at wss  [dB]

	    ALP_BUTT.setm_butterworth();
	    ALP_BUTT.butterworth(pole,zero);
		ALP_BUTT.bode_plot(pole, zero, f_start, f_end, m_min, m_max, b_title);
		ALP_BUTT.func_trans(pole, zero, ampi, &impn, amps, &stpn);
		ALP_BUTT.time_response(impn, i_min, i_max, pole, ampi, t_start, t_end, i_min, i_max, i_title);
		ALP_BUTT.time_response(stpn, s_min, s_max, pole, amps, t_start, t_end, s_min, s_max, s_title);
		ALP_BUTT.root_plot(r_yscale, pole, zero, r_title);
		ALP_BUTT.write_AnalogFilter(impn, stpn, pole, zero);
}