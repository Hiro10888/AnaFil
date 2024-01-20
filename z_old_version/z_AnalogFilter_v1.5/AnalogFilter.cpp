// AnalogFilter.cpp - ver. 1.5
// Source part of Analog Filter Design Class
// Separate AnalogFiltet.h to Analoghilter.h and AnalogFilter.cpp 
// Working with AnalogFilter.cpp and ALP_BUTT.cpp, ALP_CHEB.cpp
// Last change was made on 20230729 by udata

#ifndef MMAX
#define MMAX 100			// Maximum order of POLE,ZERO
#endif

#ifndef FDIV
#define FDIV 1000		// Frequency axis division value
#endif

#ifndef TDIV
#define TDIV 1000		// Time axis division value
#endif

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>

#include "AnalogFilter.h"
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

#define DEBUG 0


	//	Calculate Order of Butterworth LPF	(7.18),(7.19),(7.20)
	//	Given  parameter:  ws, attws, wp, attwp
	//	Gotten parameter:  m,  wc, attwc 
    void AnalogFilter::setm_butterworth()
	{
        double 	mm, g;
        g  		= pow(10.0,attwp/10.0)-1.0;
        mm 		= log((pow(10.0,attws/10.0)-1.0)/g)/(2.0*log(ws/wp));	// (7.20)
        m     	= (int)(mm+0.99999);
        attwc 	= -10.0*log10(0.5);					// -3.0103 (0.5) at wc
        wc    	= wp/pow(g,1.0/(2.0*(double)m));	//（7.19）
    }


	//	Design of Butterworth LPF	(7.12)
	//	Given Parameters 	m, wc;
	//  Gotten Parameters	pole, zero 
    void AnalogFilter::butterworth(complex<double> pole[], complex<double> zero[])
	{
        int 	i, j;
		double 	paipm, st, r;
	    double gain;
	    fgain = 1.0;
	    eps   = 0.0;
	    rmax  = rmin = wc;
	    paipm      = M_PI/(2.0*(double)m);
	    j=0;
	    for(i=0; i<2*m; i++)
	    {
		    st = paipm*(double)(m+2*i+1);
		    r  = wc*cos(st);if(r>0.0) continue;
		    pole[j].real(r);				// (7.12b)
			pole[j].imag(wc*sin(st));		// (7.12b)

#if DEBUG == 1
		cout << "pole" << j << " = " << pole[j] << endl;
#endif

		    j++;
	    }
	    polem = m;
	    getgain(pole, polem, &gain);
	    gainb = gain;
	    for(i=0; i<m; i++)
		{
			zero[i] = complex<double>(0.0, 0.0);	// zeros do not exist in Butterworth
	        zerom = 0;
	        gains = 1.0;

#if DEBUG == 1
		cout << "zero" << i << " = " << zero[j] << endl;
#endif

		}
    }


	//	Gain at DC (w=0)  (double)
    void AnalogFilter::getgain(complex<double> root[], int m, double *gain)
	{
       	complex<double> k;
       	resp(0.0,root,m,&k);
       	*gain=k.real();

#if DEBUG == 1
		cout << "*gain" << " = " << *gain << endl;
#endif
    }


	//	Gain at frequency w  (Complex)		(7.89) 
    void AnalogFilter::resp(double w, complex<double> root[], int n, complex<double> *gain)
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


//	Calculate Order of Chebyshev LPF
//	Given  parameter:	wp, attwp, ws, attws
//	Gotten parameter:	m, wc, attwc
void AnalogFilter::setm_chebyshev()		// add from ver 1.5
{
	double mm,sq;
	sq    = (pow(10.0,attws/10.0)-1.0)/(pow(10.0,attwp/10.0)-1.0);
	mm    = acosh(sqrt(sq))/acosh(ws/wp);
	m     = (int)(mm+0.99999);
	wc    = wp;
	attwc = attwp;
}


//	Design of Chebyshev LPF
//	Given parameter		m, attwp, wp
//  Gotten Parameters	pole, zero 
void AnalogFilter::chebyshev(complex<double> pole[], complex<double> zero[])	// add from ver 1.5
{

	int 	i,j;
	double 	a, b, r, mxr, mnr, g, paipm, st;
	double 	gain;

	paipm     = M_PI/(2.0*(double)m);
	g         = pow(10.0,attwp/20.0);
	if(m%2==1)	fgain=1.0;
	else 		fgain=1.0/g;
	eps  = sqrt(g*g-1.0);
	b         = 1.0/eps;
	a         = asinh(b)/(double)m;
	rmax = cosh(a);
	rmin = sinh(a);
	mxr       = rmax*wp;
	mnr       = rmin*wp;

	j=0;
	for(i=0; i<2*m; i++)
	{
		st = paipm*(2.0*(double)i+1.0);
		r  = mnr*sin(st);
		if(r>0.0) continue;
		pole[j].real(r);
		pole[j].imag(mxr*cos(st));

#if DEBUG == 1
		cout << "pole" << j << " = " << pole[j] << endl;
#endif

		j++;
	}
	polem = m;
	getgain(pole, polem, &gain);
	gainb = gain;
	for(i=0; i<m;i++)
	{
		zero[i] = complex<double>(0.0, 0.0);
		zerom = 0;
		gains = 1.0;
		rmax  = mxr;
		rmin  = mnr;

#if DEBUG == 1
		cout << "zero" << i << " = " << zero[i] << endl;
#endif

	}
}


	//  Bode plot
	void AnalogFilter::bode_plot(complex<double> pole[], complex<double> zero[],
	    		   double f_start, double f_end, double m_min, double m_max, const string& title)
	{
		double w, gain, amp, phs;
    	vector<double> lf(FDIV), f(FDIV), am(FDIV), ph(FDIV);
		complex<double> aa, bb, gg;		

		gain = fgain * gainb * gains;		// (7.88)
		for(int i=0; i < FDIV; ++i){
        	lf.at(i) = log10(f_start) + i*log10(f_end/f_start)/FDIV;
        	f.at(i) = pow(10, lf.at(i));

			w    = 2.0 * M_PI * f.at(i);
			resp(w, pole, polem, &aa);
			resp(w, zero, zerom, &bb);
			gg = bb / aa;
			gg.real(gg.real()*gain);
			gg.imag(gg.imag()*gain);
			
			amp = 10.0*log10(gg.real()*gg.real() + gg.imag()*gg.imag());
			am.at(i) = amp;

			phs = atan2(gg.imag(),gg.real())*180.0/M_PI;
			ph.at(i) = phs;

#if DEBUG == 1
//		cout << "aa" << i << " = " << aa << endl;
//		cout << "bb" << i << " = " << bb << endl;
//		cout << "gg" << i << " = " << gg << endl;
//		cout << "amp" << i << " = " << amp << endl;
//		cout << "phs" << i << " = " << phs << endl;
#endif

		}

		// Plot a bode plot

	    plt::figure_size(1200, 780);
    	plt::subplot2grid(2, 1, 0, 0);
    	plt::named_semilogx("Magnitude",f,am,"b");
    	plt::title(title);
    	plt::ylabel("Magnitude  [dB]");
    	plt::xlim(f_start, f_end);
    	plt::ylim(m_min, m_max);
    	plt::legend();
    	plt::grid("true");
//		plt::pause(1.01);

   		plt::subplot2grid(2, 1, 1, 0);
   		plt::named_semilogx("Phase",f,ph,"r");
   		plt::xlabel("Frequency  [Hz]");
    	plt::ylabel("Phase  [deg]");
		plt::xlim(f_start, f_end);
   		plt::ylim(-180, 180);
   		plt::legend();
   		plt::grid("true");
//		plt::pause(1.01);

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


	//	Partial function expansion		(7.90), (7.91)
	//	Gotten Parameter for impulse response		ampi[], impn
	//	                 for step    response		amps[], stpn
	void AnalogFilter::func_trans(complex<double> pole[], complex<double> zero[],
					complex<double> ampi[], int *impn, complex<double> amps[], int *stpn)
	{
	int 	i, j;
	complex<double>	a, b, s;
	double 	r2, g;
	g  = fgain * gainb * gains;		// (7.88)
	b.real(1.0);
	b.imag(0.0);
	for(i=0; i<polem; i++)
	{
		ampi[i].real(g);
		ampi[i].imag(0.0);
		a.real(-pole[i].real());
		a.imag(-pole[i].imag());
		b=(b*a);
		for(j=0; j<polem; j++)
		{
			if(i==j) continue;
			a.real(pole[i].real() - pole[j].real());
			a.imag(pole[i].imag() - pole[j].imag());
			ampi[i] = ampi[i]/a;
		}
		for(j=0; j<zerom; j++)
		{
			a.real(pole[i].real() - zero[j].real());
			a.imag(pole[i].imag() - zero[j].imag());
			ampi[i] = ampi[i]*a;
		}
		amps[i] = ampi[i]/pole[i];
	}
	*impn   = polem;
	j       = 1;
	s.real(1.0);
	s.imag(0.0);
	for(i=0; i<zerom; i++)
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
		pole[polem].real(0.0);
		pole[polem].imag(0.0);
		amps[polem].real(g*s.real()/b.real());
		amps[polem].imag(0.0);
	}
	*stpn = *impn+j;	
	}


	void AnalogFilter::time_response(int n, double r_min, double r_max, complex<double> pole[], complex<double> amp[],
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

	    // Plot a Impulse response
		plt::figure_size(1200, 780);
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


	double AnalogFilter::time_func(double t, complex<double> pole[], complex<double> amp[], int n)
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


	//	root plot	Pole(White), Zero(Yellow) - plot
	void AnalogFilter::root_plot(double r_yscale, complex<double> pole[], complex<double> zero[], const string& title)
	{
		int i, j;
    	vector<double> xrmax(361), yrmax(361), xrmin(361), yrmin(361), p(polem), q(polem), r(zerom), s(zerom);

		for(i=0; i<=360; i++)
		{
			xrmax.at(i) = rmax * sin(i*2.0*M_PI/360.0);
			yrmax.at(i) = rmax * cos(i*2.0*M_PI/360.0);
		}		

		for(i=0; i<=360; i++)
		{
			xrmin.at(i) = rmin * sin(i*2.0*M_PI/360.0);
			yrmin.at(i) = rmin * cos(i*2.0*M_PI/360.0);
		}	

		for (int i=0; i<polem; i++)
		{
			p.at(i) = pole[i].real();
			q.at(i) = pole[i].imag();
		}

		for (int i=0; i<zerom; i++)
		{
			r.at(i) = zero[i].real();
			s.at(i) = zero[i].imag();
		}

		    plt::figure_size(780, 780);

	    	// Plot a Pole, Zero
	    	plt::plot(xrmax,yrmax,":k");
	    	plt::plot(xrmin,yrmin,":k");
		    plt::named_plot("Pole",p, q, "gx");
			plt::named_plot("Zero",r, s, "ro");			
    		plt::title(title);
	   		plt::xlabel("Real");
	    	plt::ylabel("Imaginary");
    		plt::xlim(-r_yscale, r_yscale);
	    	plt::ylim(-r_yscale, r_yscale);
			plt::legend();		
	    	plt::grid("true");
			plt::axis("equal");

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


	void AnalogFilter::write_AnalogFilter()
	{
	    ofstream ofs("ANALOGFILTER.TXT");
		ofs << "isf  "  << " = " << isf   << endl;
		ofs << "is   "  << " = " << is    << endl;
		ofs << "m	 "  << " = " << m	  << endl;
		ofs << "wc   "  << " = " << wc    << endl;
		ofs << "attwc"  << " = " << attwc << endl;
		ofs << "ws   "  << " = " << ws    << endl;
		ofs << "attws"  << " = " << attws << endl;
		ofs << "wp   "  << " = " << wp    << endl;
		ofs << "attwp"  << " = " << attwp << endl;
		ofs << "fgain"  << " = " << fgain << endl;
		ofs << "eps  "  << " = " << eps   << endl;
		ofs << "rmax "  << " = " << rmax  << endl;
		ofs << "rmin "  << " = " << rmin  << endl;
		ofs << "polem"  << " = " << polem << endl;
		ofs << "zerom"  << " = " << zerom << endl;
		ofs << "gainb"  << " = " << gainb << endl;
		ofs << "gains"  << " = " << gains << endl;
		ofs << "wp1  "  << " = " << wp1   << endl;
		ofs << "wp2  "  << " = " << wp2   << endl;
		ofs << "ws1  "  << " = " << ws1   << endl;
		ofs << "ws2  "  << " = " << ws2   << endl;
		ofs << "bandw"  << " = " << bandw << endl;
		ofs << "w0   "  << " = " << w0    << endl;
	}
