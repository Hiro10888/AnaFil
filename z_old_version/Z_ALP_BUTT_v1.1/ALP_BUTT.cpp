// ALP_BUTT.cpp - ver. 1.1
// last change was made on 20230723 by udata

#include "AnalogFilter.h"

int main() {
	const 	string b_title 	= "Bode plot of Butterworth LPF";
	double	f_start 		= 0.1;
	double	f_end   		= 10.0;
	double	m_min			= -200.0;
	double	m_max			= 50.0;

	double	t_start 		= 0.0;
	double	t_end   		= 10.0;

	const 	string i_title	= "Impulse response of Butterworth LPF";
	double	i_min			= -3.0;
	double	i_max			= 3.0;

	const 	string s_title	= "Step response of Butterworth LPF";
	double	s_min			= -0.5;
	double	s_max			= 1.5;

	const 	string r_title	= "Poles and Zeros  of Butterworth LPF";
	double	r_yscale    	= 20.0;

	int			    i, impn, stpn;
	makefilter	    fil;
	complex<double> pole[MMAX], zero[MMAX], ampi[MMAX], amps[MMAX];


			fil.isf		= 1;            		/* 1=LPF, 2=BPF,  3=HPF, 4=BEF	 */
			fil.is		= 1;					/* 1=butterworth, 2=Chebyshev 	 */
												/* 3=Inv. Cheby.  4=Elliptic  	 */
			fil.wp		= 1.0*2.0*M_PI;		 	/* Pass band frequency [rad/sec] */
			fil.attwp  	= 1.0;					/* Attenuation at wpp  [dB]      */
			fil.ws 		= 2.0*2.0*M_PI;		 	/* Stop band frequency [rad/sec] */
			fil.attws	= 60.0;  				/* Attenuation at wss  [dB]      */


    Analogfilter analogfilter;

    analogfilter.setm_butterworth(&fil);
    analogfilter.butterworth(&fil,pole,zero);
	analogfilter.bode_plot(&fil, pole, zero, f_start, f_end, m_min, m_max, b_title);
	analogfilter.func_trans(&fil, pole, zero, ampi, &impn, amps, &stpn);
	analogfilter.time_response(impn, i_min, i_max, pole, ampi, t_start, t_end, i_min, i_max, i_title);
	analogfilter.time_response(stpn, s_min, s_max, pole, amps, t_start, t_end, s_min, s_max, s_title);
	analogfilter.root_plot(&fil, r_yscale, pole, zero, r_title);
	analogfilter.write_makefilter(&fil);

}