#define _USE_MATH_DEFINES
#define MMAX 50
#define FDIV 2000
#define TDIV 2000
#define RDIV 300
#define NPOL 11

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

#define DEBUG 0


struct makefilter {
	int	isf;			/* Filter type				1 = LPF					*/
						/* 							2 = BPF					*/
						/*							3 = HPF					*/
						/*							4 = BEF					*/
	int	is;				/* Transfer function		1 = Butterworth			*/
						/*							2 = Chebyshev			*/
						/*                      	3 = Inverse Chevyshev	*/
						/*                      	4 = Elliptic			*/
						/* For LPF Design									*/
	int 	m;			/* Filter order  									*/
	double  wc;			/* Nutural frequency  			[rad/sec]			*/
	double  attwc;		/* Attenuation at wc  			[dB]				*/
	double  ws;			/* Stop band frequency          [rad/sec]			*/
	double  attws;		/* Attenuation at ws  			[dB]				*/
	double  wp;			/* Pass band frequency          [rad/sec]			*/
	double  attwp;		/* Attenuation at wp  			[dB]				*/
						/* For Chebyshev Filter Design						*/
	double  fgain;		/* Gain caused by Filter Type  						*/
	double  eps;		/* Epsilon											*/
	double  rmax;		/* Radius of outside circule (normalized wc)		*/
	double  rmin;		/* Radius of inside circule  (normarized wc)		*/
						/* For Transfer Function Definition					*/
	int 	polem;		/* Order of denominator polynomial					*/
	int 	zerom;		/* Order of numerator polynomial					*/
	double 	gainb;		/* Gain caused by denominator 	[LPF] 				*/
	double 	gains;		/* Gain caused by numerator   	[LPF]  				*/
						/* For Band pass filter (BPF,BEF) Design			*/
	double 	wp1;		/* Lowest frequency of Pass band					*/
	double 	wp2;		/* Highest frequency of Pass band					*/
	double 	ws1;		/* Define Attenuation of lowest stopband freq.		*/
	double 	ws2;		/* Define Attenyation of highest stopband freq.		*/
	double 	bandw;		/* Band width of Band pass filter					*/
	double 	w0;			/* Center frequency									*/
	};

class Analogfilter {
    public:

    //	Calculate Order of Butterworth LPF	(7.18),(7.19),(7.20)
	//	Given  parameter:  par->ws, par->attws, par->wp, par->attwp
    //	Gotten parameter:  par->m,  par->wc,    par->attw
    void setm_butterworth(makefilter *par){
        double 	mm, g;
        g  = pow(10.0,par->attwp/10.0)-1.0;
        mm = log((pow(10.0,par->attws/10.0)-1.0)/g)/(2.0*log(par->ws/par->wp));
        par->m     = (int)(mm+0.99999);
        par->attwc = -10.0*log10(0.5);
        par->wc    = par->wp/pow(g,1.0/(2.0*(double)par->m));

#if DEBUG == 1
        cout << "ws   "   << " = " << par->ws << endl;
        cout << "attws"   << " = " << par->attws << endl;
        cout << "wp   "   << " = " << par->wp << endl;
        cout << "attwp"   << " = " << par->attwp << endl;
        cout << "m    "   << " = " << par->m << endl;
        cout << "wc   "   << " = " << par->wc << endl;
        cout << "attwc"   << " = " << par->attwc << endl;
#endif	

    }

    //	Design of Butterworth LPF	(7.12)
    void butterworth(makefilter *par, complex<double> pole[], complex<double> zero[]){
        int 	i, j;
		double 	paipm, st, r;
	    double gain;
	    par->fgain = 1.0;
	    par->eps   = 0.0;
	    par->rmax  = par->rmin = par->wc;
	    paipm      = M_PI/(2.0*(double)par->m);
	    j=0;
	    for(i=0; i<2*par->m; i++)
	    {
		    st = paipm*(double)(par->m+2*i+1);
		    r  = par->wc*cos(st);if(r>0.0) continue;
		    pole[j].real(r);
		pole[j].imag(par->wc*sin(st));

#if DEBUG == 1
		cout << "pole" << j << " = " << pole[j] << endl;
#endif

		    j++;
	    }
	    par->polem = par->m;
	    getgain(pole, par->polem, &gain);
	    par->gainb = gain;
	    for(i=0; i<par->m; i++)
		{
			zero[i] = complex<double>(0.0, 0.0);
	        par->zerom = 0;
	        par->gains = 1.0;
		}
    }

    //	Gain at DC (w=0)  (double)
    void getgain(complex<double> root[], int m, double *gain)
	{
       	complex<double> k;
       	resp(0.0,root,m,&k);
       	*gain=k.real();

#if DEBUG == 1
		cout << "*gain" << " = " << *gain << endl;
#endif

    }

    //	Gain at frequency w  (Complex)		(7.89)
    void resp(double w, complex<double> root[], int n, complex<double> *gain)
	{
       	complex<double> f;
       	int i;
       	gain->real(1.0); gain->imag(0.0);
       	for(i=0;i<n;i++)
       	{
       		f.real(-root[i].real());
       		f.imag(w-(root[i].imag()));
       		*gain *= f;
       	}
    }


	// bodep_lot
	void bode_plot(makefilter *par, complex<double> pole[], complex<double> zero[],
	    		   double f_start, double f_end, double m_min, double m_max, const string& b_title)
	{
		double w, gain, amp, phs;
    	vector<double> lf(FDIV), f(FDIV), am(FDIV), ph(FDIV);
		complex<double> aa, bb, gg;		

		gain = par->fgain * par->gainb * par->gains;
		for(int i=0; i < FDIV; ++i){
        	lf.at(i) = log10(f_start) + i*log10(f_end/f_start)/FDIV;
        	f.at(i) = pow(10, lf.at(i));
			w    = 2.0 * M_PI * f.at(i);
			resp(w, pole, par->polem, &aa);
			resp(w, zero, par->zerom, &bb);
			gg = bb / aa;
			gg.real(gg.real()*gain);
			gg.imag(gg.imag()*gain);
			amp = 10.0*log10(gg.real()*gg.real() + gg.imag()*gg.imag());
			am.at(i) = amp;

			phs = atan2(gg.imag(),gg.real())*180.0/M_PI;
			ph.at(i) = phs;
/*
#if DEBUG == 1
	cout << "lf" << i << " = " << lf.at(i) << endl;
	cout << "f" << i << " = " << f.at(i) << endl;
	cout << "am" << i << " = " << am.at(i) << endl;
	cout << "ph" << i << " = " << ph.at(i) << endl;    
#endif	
*/
		}

	    plt::figure_size(1200, 780);

    	// Plot a Magnitude
    	// subplot2grid(nrows, ncols, row, col);
    	plt::subplot2grid(2, 1, 0, 0);
    	plt::named_semilogx("Magnitude",f,am,"b");
    	plt::title(b_title);
    	plt::ylabel("Magnitude  [dB]");
    	plt::xlim(f_start, f_end);
    	plt::ylim(m_min, m_max);
    	plt::legend();
    	plt::grid("true");

   		// Plot a Phase
   		plt::subplot2grid(2, 1, 1, 0);
   		plt::named_semilogx("Phase",f,ph,"r");
   		plt::xlabel("Frequency  [Hz]");
    	plt::ylabel("Phase  [deg]");
		plt::xlim(f_start, f_end);
   		plt::ylim(-180, 180);
   		plt::legend();
   		plt::grid("true");

#if DEBUG == 0
    plt::show();
#endif	

#if DEBUG == 1
    // save figure
    const char* filename = "Bode plot of Butterworth LPF";
    std::cout << "Saving result to " << filename << std::endl;;
    plt::save(filename);
#endif

	}


	void func_trans(makefilter *par, complex<double> pole[], complex<double> zero[],
					complex<double> ampi[], int *impn, complex<double> amps[], int *stpn)
	{
	int 	i, j;
	complex<double>	a, b, s;
	double 	r2, g;
	g      = par->fgain * par->gainb * par->gains;
	b.real(1.0);
	b.imag(0.0);
	for(i=0; i<par->polem; i++)
	{
		ampi[i].real(g);
		ampi[i].imag(0.0);
		a.real(-pole[i].real());
		a.imag(-pole[i].imag());
		b=(b*a);
		for(j=0; j<par->polem; j++)
		{
			if(i==j) continue;
			a.real(pole[i].real() - pole[j].real());
			a.imag(pole[i].imag() - pole[j].imag());
			ampi[i] = ampi[i]/a;
		}
		for(j=0; j<par->zerom; j++)
		{
			a.real(pole[i].real() - zero[j].real());
			a.imag(pole[i].imag() - zero[j].imag());
			ampi[i] = ampi[i]*a;
		}
		amps[i] = ampi[i]/pole[i];
	}
	*impn   = par->polem;
	j       = 1;
	s.real(1.0);
	s.imag(0.0);
	for(i=0; i<par->zerom; i++)
	{
		r2 = zero[i].real()*zero[i].real()+zero[i].imag()*zero[i].imag();
		a.real(-zero[i].real());
		a.imag(-zero[i].imag());
		s = s*a;
		if(r2 < 1.0e-8)
		{
			j=0; break;
		}
	}
	if(j == 1)
	{
		pole[par->polem].real(0.0);
		pole[par->polem].imag(0.0);
		amps[par->polem].real(g*s.real()/b.real());
		amps[par->polem].imag(0.0);
	}
	*stpn = *impn+j;	
	}

	/*	Time response (Impulse,Step response, ...etc) 	plot		*/
	void time_response(int n, double r_min, double r_max, complex<double> pole[], complex<double> amp[],
						double t_start, double t_end, double min, double max, const string& title)
	{
		int i = 0;
		double 	t, p;
    	vector<double> tr(TDIV), pr(TDIV);
		t = t_start;
		for (int i=0; i<TDIV; i++)
		{
			t = t_start + ((t_end-t_start)/TDIV)*i;
			p  = time_func(t, pole, amp, n);
			tr.at(i) = t;
			pr.at(i) = p;
		}

	    plt::figure_size(1200, 780);

    	// Plot a Impulse response
    	plt::plot(tr,pr,"b");
    	plt::title(title);
   		plt::xlabel("Time  [sec]");
    	plt::ylabel("Amplitude");
    	plt::xlim(t_start, t_end);
    	plt::ylim(r_min, r_max);		
    	plt::grid("true");

#if DEBUG == 0
    plt::show();
#endif	

#if DEBUG == 1
    // save figure
    const string& filename = title;
    std::cout << "Saving result to " << filename << std::endl;;
    plt::save(filename);
#endif


	}

	/*	Inverse Laplace transform	(7.92), (7.93)			*/
	double time_func(double t, complex<double> pole[], complex<double> amp[], int n)
	{
		int 	i;
		double 	a, si, co, wt, rl, ig;
		rl = ig = 0.0;
		for(i=0; i<n; i++)
		{
			a   = exp(t * pole[i].real());
			wt  = t * pole[i].imag();
			co  = cos(wt);
			si  = sin(wt);
			rl += a * (amp[i].real() * co - amp[i].imag() * si);
			ig += a * (amp[i].imag() * co + amp[i].real() * si);
		}
		return(rl);
	}


	/*	root plot	Pole(White), Zero(Yellow) - plot		*/
	void root_plot(makefilter *par, double r_yscale, complex<double> pole[], complex<double> zero[], const string& title)
	{
		int i, j;
    	vector<double> x(RDIV), y(RDIV), p(par->polem), q(par->polem), r(par->zerom), s(par->zerom);

		for(i=0; i<RDIV; i++)
		{
			x.at(i) = par->rmax * sin(i*2.0*M_PI/300.0);
			y.at(i) = par->rmax * cos(i*2.0*M_PI/300.0);

		}		

		for (int i=0; i<par->polem; i++)
		{
			p.at(i) = pole[i].real();
			q.at(i) = pole[i].imag();
		}

		for (int i=0; i<par->zerom; i++)
		{
			r.at(i) = zero[i].real();
			s.at(i) = zero[i].imag();
		}

		    plt::figure_size(780, 780);

	    	// Plot a Pole, Zero
	    	plt::plot(x,y,":k");
		    plt::named_plot("Pole",p, q, "gx");
			plt::named_plot("Zero",r, s, "ro");			
    		plt::title(title);
	   		plt::xlabel("Real");
	    	plt::ylabel("Imaginary");
    		plt::xlim(-8, 8);
	    	plt::ylim(-8, 8);
			plt::legend();		
	    	plt::grid("true");

#if DEBUG == 0
    plt::show();
#endif	

#if DEBUG == 1
    // save figure
    const string& filename = title;
    std::cout << "Saving result to " << filename << std::endl;;
    plt::save(filename);
#endif

	}


	void write_makefilter(makefilter *par)
	{
	    ofstream ofs("MAKEFILTER.TXT");
		ofs << "isf  "  << " = " << par->isf   << endl;
		ofs << "is   "  << " = " << par->is    << endl;
		ofs << "m	  "  << " = " << par->m	    << endl;
		ofs << "wc   "  << " = " << par->wc    << endl;
		ofs << "attwc"  << " = " << par->attwc << endl;
		ofs << "ws   "  << " = " << par->ws    << endl;
		ofs << "attws"  << " = " << par->attws << endl;
		ofs << "wp   "  << " = " << par->wp    << endl;
		ofs << "attwp"  << " = " << par->attwp << endl;
		ofs << "fgain"  << " = " << par->fgain << endl;
		ofs << "eps  "  << " = " << par->eps   << endl;
		ofs << "rmax "  << " = " << par->rmax  << endl;
		ofs << "rmin "  << " = " << par->rmin  << endl;
		ofs << "polem"  << " = " << par->polem << endl;
		ofs << "zerom"  << " = " << par->zerom << endl;
		ofs << "gainb"  << " = " << par->gainb << endl;
		ofs << "gains"  << " = " << par->gains << endl;
		ofs << "wp1  "  << " = " << par->wp1   << endl;
		ofs << "wp2  "  << " = " << par->wp2   << endl;
		ofs << "ws1  "  << " = " << par->ws1   << endl;
		ofs << "ws2  "  << " = " << par->ws2   << endl;
		ofs << "bandw"  << " = " << par->bandw << endl;
		ofs << "w0   "  << " = " << par->w0    << endl;
	}


/*
void write_makefilter(makefilter *fil)
{
    ofstream ofs("MAKEFILT.BIN", ios::binary);
    ofs.write(reinterpret_cast<char*>(fil), sizeof(makefilter));
}

void write_makefilter(makefilter *fil)
{	FILE	*fp;
	if((fp=fopen("G:\MAKEFILT.BIN","wb")) == NULL)
	{	fprintf(stderr, "!!! Can't open G:\MAKEFILT.BIN to write\n");
			exit(1);
	}
	fwrite(fil, sizeof *fil, 1, fp);
	fclose(fp);
}

void write_pole(makefilter *fil, complex *pole)
{	FILE	*fp;
	if((fp=fopen("G:\POLE.BIN","wb")) == NULL)
	{	fprintf(stderr, "!!! Can't open G:\POLE.BIN to write\n");
			exit(1);
	}
	fwrite(pole, sizeof *pole, fil->polem, fp);
	fclose(fp);
}

void write_zero(makefilter *fil, complex *zero)
{	FILE	*fp;
	if((fp=fopen("G:\ZERO.BIN","wb")) == NULL)
	{	fprintf(stderr, "!!! Can't open G:\ZERO.BIN to write\n");
			exit(1);
	}
	fwrite(zero, sizeof *zero, fil->zerom, fp);
	fclose(fp);
*/




};


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