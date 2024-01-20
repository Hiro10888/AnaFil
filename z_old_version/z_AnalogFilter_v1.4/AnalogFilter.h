// AnalogFilter.h - ver. 1.4
// Headder part of Analog Filter Design Class
// Separate AnalogFiltet.h to Analoghilter.h and AnalogFilter.cpp 
// Working with AnalogFilter.cpp and ALP_BUTT.cpp
// Last change was made on 20230726 by udata

#include <string>
#include <complex>

class AnalogFilter 
{
						// For LPF Design
	int 	m;			// Filter order
	double  wc;			// Nutural frequency  			[rad/sec]
	double  attwc;		// Attenuation at wc  			[dB]

						// For Chebyshev Filter Design
	double  fgain;		// Gain caused by Filter Type
	double  eps;		// Epsilon
	double  rmax;		// Radius of outside circule (normalized wc)
	double  rmin;		// Radius of inside circule  (normarized wc)

						// For Transfer Function Definition
	int 	polem;		// Order of denominator polynomial
	int 	zerom;		// Order of numerator polynomial
	double 	gainb;		// Gain caused by denominator 	[LPF]
	double 	gains;		// Gain caused by numerator   	[LPF]

						// For Band pass filter (BPF,BEF) Design
	double 	wp1;		// Lowest frequency of Pass band
	double 	wp2;		// Highest frequency of Pass band
	double 	ws1;		// Define Attenuation of lowest stopband freq.
	double 	ws2;		// Define Attenyation of highest stopband freq.
	double 	bandw;		// Band width of Band pass filter
	double 	w0;			// Center frequency

    public:

	int	isf;			// Filter type				1 = LPF
						// 							2 = BPF
						//							3 = HPF
						//							4 = BEF
	int	is;				// Transfer function		1 = Butterworth
						//							2 = Chebyshev
						//                      	3 = Inverse Chevyshev
						//                      	4 = Elliptic
	double  wp;			// Pass band frequency          [rad/sec]
	double  attwp;		// Attenuation at wp  			[dB]
	double  ws;			// Stop band frequency          [rad/sec]
	double  attws;		// Attenuation at ws  			[dB]


    void setm_butterworth();
	void butterworth(std::complex<double> pole[], std::complex<double> zero[]);
	void bode_plot(std::complex<double> pole[], std::complex<double> zero[],
	    		   double f_start, double f_end, double m_min, double m_max, const std::string& b_title);	
	void func_trans(std::complex<double> pole[], std::complex<double> zero[],
					std::complex<double> ampi[], int *impn, std::complex<double> amps[], int *stpn);
	void time_response(int n, double r_min, double r_max, std::complex<double> pole[], std::complex<double> amp[],
						double t_start, double t_end, double i_min, double i_max, const std::string& i_title);
	void root_plot(double r_yscale, std::complex<double> pole[], std::complex<double> zero[], const std::string& r_title);
	void write_AnalogFilter();

    void getgain(std::complex<double> root[], int m, double *gain);	
    void resp(double w, std::complex<double> root[], int n, std::complex<double> *gain);
	double time_func(double t, std::complex<double> pole[], std::complex<double> amp[], int n);
};
