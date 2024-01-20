// ALP_ICHE.cpp - ver. 1.0
// Working with Analoghilter.h_v1.6 and AnalogFilter.cpp_v1.6
// last change was made on 20230731 by udata

#define MMAX 50			// Maximum order of POLE,ZERO

#include <string>
#include <complex>
#include "AnalogFilter.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

int main()
{
    const   std::string b_title     = "Bode plot of Inv_Chebyshev LPF";
    double  f_start                 = 0.1;
    double  f_end                   = 100.0;
    double  m_min                   =-150.0;
    double  m_max                   = 50.0;
    double  t_start                 = 0.0;
    double  t_end                   = 5.0;
    
    const   std::string i_title     = "Impulse response of Inv_Chebyshev LPF";
    double  i_min                   =-4.0;
    double  i_max                   = 4.0;

    const   std::string  s_title    = "Step response of Inv_Chebyshev LPF";
    double  s_min                   =-0.2;
    double  s_max                   = 1.2;

    const   std::string r_title     = "Poles and Zeros of Inv_Chebyshev LPF";
    double  r_yscale                = 100.0;

    int                             i, impn, stpn;
    std::complex<double>            pole[MMAX], zero[MMAX], ampi[MMAX], amps[MMAX];
    
    AnalogFilter ALP_ICHE;
    
        ALP_ICHE.isf     = 1;           // 1=LPF, 2=BPF,  3=HPF, 4=BEF
        ALP_ICHE.is      = 3;           // 1=butterworth, 2=Chebyshev    
                                        // 3=Inv. Cheby.  4=Elliptic     
        ALP_ICHE.wp      = 10.0;        // Pass band frequency [rad/sec] 
        ALP_ICHE.attwp   = 1.0;         // Attenuation at wpp  [dB]      
        ALP_ICHE.ws      = 20.0;        // Stop band frequency [rad/sec] 
        ALP_ICHE.attws   = 60.0;        // Attenuation at wss  [dB]      

        ALP_ICHE.setm_inv_chebyshev();
        ALP_ICHE.inv_chebyshev(pole, zero);

        ALP_ICHE.bode_plot(pole, zero, f_start, f_end, m_min, m_max, b_title);
        ALP_ICHE.func_trans(pole, zero, ampi, &impn, amps, &stpn);

		ALP_ICHE.time_response(impn, i_min, i_max, pole, ampi, t_start, t_end, i_min, i_max, i_title);
		ALP_ICHE.time_response(stpn, s_min, s_max, pole, amps, t_start, t_end, s_min, s_max, s_title);
        
        ALP_ICHE.root_plot(r_yscale, pole, zero, r_title);
		ALP_ICHE.write_AnalogFilter(impn, stpn, pole, zero);
}