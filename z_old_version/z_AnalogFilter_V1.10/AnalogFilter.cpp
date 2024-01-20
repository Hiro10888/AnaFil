// AnalogFilter.cpp - ver. 1.10
// Source part of Analog Filter Design Class
// Separate AnalogFiltet.h to Analoghilter.h and AnalogFilter.cpp 
// Working with AnalogFilter.cpp and ALP_BUTT.cpp, ALP_CHEB.cpp
// Working with AnalogFilter.cpp and ALP_BUTT.cpp, ALP_CHEB.cpp, ALP_ICHE.cpp
// Working with AnalogFilter.cpp and ALP_BUTT.cpp, ALP_CHEB.cpp, ALP_ICHE.cpp, ALP_ELLI.cpp
// ... ABP_BUTT.cpp, AHP_BUTT.cpp, ABE_BUTT.cpp
// Last change was made on 2023 by udata

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
	//	Given  parameter:  wp, attwp, ws, attws
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


//	Calculate Order of Chebyshev LPF (7.38)
//	Given  parameter:	wp, attwp, ws, attws
//	Gotten parameter:	m, wc, attwc
void AnalogFilter::setm_chebyshev()				// add from ver 1.5
{
	double mm,sq;
	sq    = (pow(10.0,attws/10.0)-1.0)/(pow(10.0,attwp/10.0)-1.0);
	mm    = acosh(sqrt(sq))/acosh(ws/wp);		// (7.38)
	m     = (int)(mm+0.99999);
	wc    = wp;		// Characteristics of Chebyshev
	attwc = attwp;
}


//	Design of Chebyshev LPF						(7.32), (7.35)
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
	mxr       = rmax*wc;
	mnr       = rmin*wc;

	j=0;
	for(i=0; i<2*m; i++)
	{
		st = paipm*(2.0*(double)i+1.0);			// (7.32) st=ak
		r  = mnr*sin(st);
		if(r>0.0) continue;
		pole[j].real(r);						// (7.35)
		pole[j].imag(mxr*cos(st));				// (7.35)

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


//	Calculate Order of Inverse Chebyshev LPF
//	Given  parameter:	wp, attwp, ws, attws
//	Gotten parameter:	m, wc, attwc
void AnalogFilter::setm_inv_chebyshev()
{
	double 	mm, sq;
	sq=(pow(10.0,attws/10.0)-1.0)/(pow(10.0,attwp/10.0)-1.0);
	mm    = acosh(sqrt(sq))/acosh(ws/wp);
	m     = (int)mm;
	wc    = ws;
	attwc = attws;
}



//	Design of Inverse Chebyshev LPF
//	Given parameter: m, attws, ws
//  Gotten Parameters	pole, zero 
void AnalogFilter::inv_chebyshev(complex<double> pole[], complex<double> zero[])
{
	int 	i, j;
	double 	a, b, r, s, rs, mxr, mnr, paipm, st;
	double 	gain;

	paipm   = M_PI/(2.0*(double)m);
	fgain 	= 1.0;
	eps   	= sqrt(1.0/(pow(10.0,attws/10.0)-1.0));
	b       = 1.0/eps;
	a       = asinh(b)/(double)m;
	rmax  	= cosh(a);
	rmin  	= sinh(a);
	mxr     = rmax/ws;
	mnr     = rmin/ws;

	j=0;
	for(i=0; i<2*m; i++)
	{
		st	= paipm*(2.0*(double)i+1.0);
		r   = mnr*sin(st);
		if(r<0.0) continue;
		s   = mxr*cos(st);
		rs  = r*r+s*s;
		pole[j].real(-r/rs);
		pole[j].imag(s/rs);

#if DEBUG == 1
		cout << "pole" << j << " = " << pole[j] << endl;
#endif

		j++;
	}
	polem = m;
	getgain(pole, polem, &gain);
	gainb = gain;
	j=0;
	for(i=0;i<m/2;i++)
	{
		s	= ws/cos(paipm*(2.0*(double)i+1.0));
		zero[j].real(0.0);
		zero[j++].imag(s);
		zero[j].real(0.0);
		zero[j++].imag(-s);

#if DEBUG == 1
		cout << "zero" << j-2 << " = " << zero[j-2] << endl;
		cout << "zero" << j-1 << " = " << zero[j-1] << endl;
#endif

	}
	zerom = j;
	getgain(zero, zerom, &gain);
	gains = 1.0/gain;
	rmax  = ws;
	rmin  = wp;
}


//	Calculate Order of Elliptic LPF
//	Given  parameter: 	wp, attws, ws;
//	attwp = 		Fixed to 0.5[dB]
//	Gotten parameter:	m,wc(wp),attwc(=0.5)
int AnalogFilter::setm_elliptic()
{

	double lgpr[3][4] =
	{
		 8.3,21.9,36.3,50.6,
		13.9,31.2,48.6,66.1,
		21.5,42.8,64.1,85.5
	};
	int 	i, wcpn, wcn;
	wc    = wp;
	attwc = 0.5;
	wcpn       = ws/wp*10.0;
	if(wcpn!=15 && wcpn!=20 && wcpn!=30)
	{
		printf("Only can choice  ws/wc = 1.5, 2.0, 3.0\n");
		return(0);
	}
	if(wcpn==15) wcn = 0;
	else if(wcpn==20) wcn = 1;
	else if(wcpn==30) wcn = 2;
	for(i=0;i<4;i++)
		if(lgpr[wcn][i]>attws)
		{
			m = i+2;
			return(1);
		}
	m = 6;
	printf("Inpossible to Design\n");
	return(0);	
}


//	Design of elliptic LPF
//	Given  parameter: 	m, wp, attws, ws
//	attwp = 		Fixed to 0.5[dB]
int AnalogFilter::elliptic(complex<double> pole[], complex<double> zero[])
{
	double tblpr_r[3][4][5] =
	{
	 -0.515765, -0.515765, 0.0      , 0.0      , 0.0      ,
	 -0.226430, -0.226430, -0.766952, 0.0      , 0.0      ,
	 -0.126980, -0.126980, -0.460005, -0.460005, 0.0      ,
	 -0.081730, -0.081730, -0.285115, -0.285115, -0.425970,
	 -0.622520, -0.622520, 0.0      , 0.0      , 0.0      ,
	 -0.268935, -0.268935, -0.692120, 0.0      , 0.0      ,
	 -0.150580, -0.150580, -0.442280, -0.442280, 0.0      ,
	 -0.096275, -0.096275, -0.290270, -0.290270, -0.425970,
	 -0.678575, -0.678575, 0.0      , 0.0      , 0.0      ,
	 -0.294710, -0.294710, -0.652630, 0.0      , 0.0      ,
	 -0.164895, -0.164895, -0.431290, -0.431290, 0.0      ,
	 -0.105330, -0.105330, -0.292205, -0.292205, -0.374520
	};

	double tblpr_i[3][4][5] =
	{
	 1.156363, -1.156363, 0.0     , 0.0      , 0.0,
	 1.047807, -1.047807, 0.0     , 0.0      , 0.0,
	 1.021918, -1.021918, 0.510123, -0.510123, 0.0,
	 1.012527, -1.012527, 0.703363, -0.703363, 0.0,
	 1.097387, -1.097387, 0.0     , 0.0      , 0.0,
	 1.037383, -1.037383, 0.0     , 0.0      , 0.0,
	 1.019758, -1.019758, 0.463366, -0.463366, 0.0,
	 1.012300, -1.012300, 0.663885, -0.663885, 0.0,
	 1.046354, -1.046354, 0.0     , 0.0      , 0.0,
	 1.028949, -1.028949, 0.0     , 0.0      , 0.0,
	 1.017885, -1.017885, 0.438017, -0.438017, 0.0,
	 1.011932, -1.011932, 0.641096, -0.641096, 0.0
	};

	double 	tblzr[3][4][4] =
	{
	1.981679,-1.981679,0.0,0.0,
	1.675115,-1.675115,0.0,0.0,
	1.592341,-1.592341,3.478406,-3.478406,
	1.557405,-1.557405,2.331875,-2.331875,
	2.732051,-2.732051,0.0,0.0,
	2.270068,-2.270068,0.0,0.0,
	2.143189,-2.143189,4.922113,-4.922113,
	2.089246,-2.089246,3.250805,-3.250805,	// real(0.0), imag(x10)
	4.181540,-4.181540,0.0,0.0,
	3.439158,-3.439158,0.0,0.0,
	3.233490,-3.233490,7.646633,-7.646633,
	3.145719,-3.145719,5.007684,-5.007684
	};
	
	int pln[4] = {2,3,4,5};
	
	int zrn[4] = {2,2,4,4};
	
	double 	lgpr[3][4] =
	{
	 8.3,21.9,36.3,50.6,
	13.9,31.2,48.6,66.1,
	21.5,42.8,64.1,85.5
	};
	
	int 	lm, wcn, i, wcpn;
	double 	wcp;
	double 	gain;
	rmax  = ws;
	rmin  = wp;
	fgain = 1.0;

	if(m<2 || m>5)
	{
		printf("Filter Order m isl out of range (2<=m<=5)\n");
		return(0);
	}

	lm    = m-2;
	wcpn = ws/wp*10.0;
	if(wcpn!=15 && wcpn!=20 && wcpn!=30)
	{
		printf("ws/wc is only select 1.5, 2.0, 3.0\n");
		return(0);
	}

	if(wcpn==15) wcn = 0;
	else if(wcpn==20) wcn = 1;
	else if(wcpn==30) wcn = 2;
	if(lgpr[wcn][lm]>=attws)
	{
		polem = pln[lm];
		for(i=0;i<pln[lm];i++)
		{
		pole[i].real(tblpr_r[wcn][lm][i]*wp);
		pole[i].imag(tblpr_i[wcn][lm][i]*wp);

#if DEBUG == 1
		cout << "pole" << i << " = " << pole[i] << endl;
#endif

		}
		getgain(pole, polem, &gain);
		gainb = gain;
		zerom = zrn[lm];
		for(i=0;i<zrn[lm];i++)
		{
			zero[i].real(0.0);
			zero[i].imag(tblzr[wcn][lm][i]);
			zero[i].imag(zero[i].imag()*wp);

#if DEBUG == 1
		cout << "zero" << i << " = " << zero[i] << endl;
#endif

		}
		getgain(zero, zerom, &gain);
		gains = 1.0/gain;
		return(1);
	}
	else{
		printf("Order=%d, wc=%lf at ws=%lf\n",m,wp,ws);
		printf("Attenuation is %lf,\n",lgpr[wcn][lm]);
		printf("Requred Attenuation = %lf can't achive\n",attws);
		printf("Increase the Order, or ws is much higher to wc\n");
		return(0);
	}
}


	//	Gain at DC (w=0)  (double)
    void AnalogFilter::getgain(complex<double> root[], int m, double *gain)
	{
       	complex<double> k;
       	resp(0.0,root,m,&k);
       	*gain=k.real();
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


	//  Bode plot
	void AnalogFilter::bode_plot(complex<double> pole[], complex<double> zero[],
	    		   double f_start, double f_end, double m_min, double m_max, const string& title)
	{
		double w, gain, amp, phs;
    	vector<double> lf(FDIV), f(FDIV), am(FDIV), ph(FDIV);
		complex<double> aa, bb, gg;		

		gain = fgain * gainb * gains;		// (7.88)


#if DEBUG == 1
		cout << "gain" << " = " << gain << endl;
#endif

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



	//	Frequency transform from LPF to BPF (Analog filter)
	void AnalogFilter::a_lpftobpf(complex<double> pole[], complex<double> zero[])
	{
		complex<double> a;
		int 	mm,i,j;
		double 	r,g,bw,ww;
		j  = polem - zerom;
		bw = bandw;
		ww = w0;
		for(i=0, mm=polem; i<polem; i++, mm++)
		{
			r       = pole[i].real();
			g       = pole[i].imag();
			a.real(r*r*bw*bw-g*g*bw*bw-4.0*ww*ww);
			a.imag(2.0*r*g*bw*bw);
			a		= sqrt(a);
			pole[i].real((r*bw+a.real())/2.0);
			pole[i].imag((g*bw+a.imag())/2.0);
			pole[mm].real((r*bw-a.real())/2.0);
			pole[mm].imag((g*bw-a.imag())/2.0);
		}
		polem = mm;
		for(i=0, mm = zerom; i<zerom; i++, mm++)
		{
			r       	= zero[i].real();
			g       	= zero[i].imag();
			a.real(r*r*bw*bw-g*g*bw*bw-4.0*ww*ww);
			a.imag(2.0*r*g*bw*bw);
			a			= sqrt(a);
			zero[i].real((r*bw+a.real())/2.0);
			zero[i].imag((g*bw+a.imag())/2.0);
			zero[mm].real((r*bw-a.real())/2.0);
			zero[mm].imag((g*bw-a.imag())/2.0);
		}
		zerom = mm;
		if(j>0)
		{
			mm = zerom;
			for(i=0;i<j;i++)
			{
				zero[i+mm].real(0.0);
				zero[i+mm].imag(0.0);
				(gainb) *= bw;
			}
			zerom += j;
		}else if(j<0)
		{
			mm = polem;
			j  = -j;
			for(i=0; i<j; i++)
			{
				pole[i+mm].real(0.0);
				pole[i+mm].imag(0.0);
				(gains) /= bw;
			}
			polem += j;
		}
	}

	// Given  parameter: ws1, ws2, wp1, wp2, attws, attwp
	// Gotten parameter: bandw, w0, m
	void AnalogFilter::bpf_butterworth(complex<double> pole[], complex<double> zero[])
	{
		ws    = (ws2 - ws1) / (wp2 - wp1);
		wp    = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2 * wp1);
		setm_butterworth();
		butterworth(pole, zero);
		a_lpftobpf(pole, zero);
	}


	// Given  parameter 	ws1, ws2, wp1, wp2, attws, attwp
	// Gotten parameter	bandw, w0, m
	void AnalogFilter::bpf_chebyshev(complex<double> pole[], complex<double> zero[])
	{
		ws    = (ws2 - ws1) / (wp2 - wp1);
		wp    = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2 * wp1);
		setm_chebyshev();
		chebyshev(pole, zero);
		a_lpftobpf(pole, zero);
	}


	// Given  parameter: ws1, ws2, wp1, wp2, attws, attwp
	// Gotten parameter: bandw, w0, m
	void AnalogFilter::bpf_inv_chebyshev(complex<double> pole[], complex<double> zero[])
	{
		ws    = (ws2 - ws1) / (wp2 - wp1);
		wp    = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2 * wp1);
		setm_inv_chebyshev();
		inv_chebyshev(pole,zero);
		a_lpftobpf(pole, zero);
	}


	// Given  parameter: ws1, ws2, wp1, wp2, attws, attwp
	// Gotten parameter: bandw, w0, m
	void AnalogFilter::bpf_elliptic(complex<double> pole[], complex<double> zero[])
	{
		ws    = (ws2 - ws1) / (wp2 - wp1);
		wp    = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2 * wp1);
		setm_elliptic();
		elliptic(pole, zero);
		a_lpftobpf(pole, zero);
	}


	//	Frequency transform from LPF to HPF (Analog Filter)
	void AnalogFilter::a_lpftohpf(complex<double> pole[], complex<double> zero[])
	{
		complex<double> ww;
		int 	i, mm;
		double 	gn;
		ww.real(w0);
		ww.imag(0.0);
		for(i=0; i<polem; i++)
		{
			pole[i] = ww/pole[i];
		}
		for(i=0; i<zerom; i++)
		{
			zero[i] = ww/zero[i];
		}
		gainb = gains = 1.0;
		mm         = polem - zerom;
		if(mm>0)
		{
			for(i=zerom; i<zerom + mm; i++)
				zero[i].real(0.0);
				zero[i].imag(0.0);
			zerom += mm;
		}
		else if(mm<0)
		{
			for(i=polem; i<zerom - mm; i++)
				pole[i].real(0.0);
				pole[i].imag(0.0);
			polem -= mm;
		}
	}


	// Given  parameter: 	ws, wp, attws, attwp
	// Gotten parameter:	w0, m
	void AnalogFilter::hpf_butterworth(complex<double> pole[], complex<double> zero[])
	{
		double 	stcws, stcwp;
		stcws   = ws;
		stcwp   = wp;
		w0 = wp;
		ws = wp/ws;
		wp = 1.0;
		setm_butterworth();
		butterworth(pole, zero);
		a_lpftohpf(pole, zero);
		ws = stcws;
		wp = stcwp;
	}


	// Given  parameter: 	ws, wp, attws, attwp
	// Gotten parameter:	w0, m
	void AnalogFilter::hpf_chebyshev(complex<double> pole[], complex<double> zero[])
	{
		double 	stcws, stcwp;
		stcws   = ws;
		stcwp   = wp;
		w0 = wp;
		ws = wp/ws;
		wp = 1.0;
		setm_chebyshev();
		chebyshev(pole, zero);
		a_lpftohpf(pole, zero);
		ws = stcws;
		wp = stcwp;
	}


	// Given  parameter: ws, wp, attws, attwp
	// Gotten parameter: w0, m
	void AnalogFilter::hpf_inv_chebyshev(complex<double> pole[], complex<double> zero[])
	{
		double 	stcws, stcwp;
		stcws   = ws;
		stcwp   = wp;
		w0 = wp;
		ws = wp/ws;
		wp = 1.0;
		setm_inv_chebyshev();
		inv_chebyshev(pole, zero);
		a_lpftohpf(pole, zero);
		ws = stcws;
		wp = stcwp;
	}


	// Given  parameter: 	ws, wp, attws, attwp
	// Gotten parameter:	w0, m
	void AnalogFilter::hpf_elliptic(complex<double> pole[], complex<double> zero[])
	{
		double 	stcws, stcwp;
		stcws   = ws;
		stcwp   = wp;
		w0 = wp;
		ws = wp/ws;
		wp = 1.0;
		setm_elliptic();
		elliptic(pole, zero);
		a_lpftohpf(pole, zero);
		ws = stcws;
		wp = stcwp;
	}


	// Given  parameter: 	ws1, ws2, wp1, wp2, attws, attwp
	// Gotten parameter:	bandw, w0, m
	void AnalogFilter::bef_butterworth(complex<double> pole[], complex<double> zero[])
	{
		ws    = (wp2 - wp1)/(ws2 - ws1);
		wp    = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2*wp1);
		setm_butterworth();
		butterworth(pole, zero);
		a_lpftobef(pole, zero);
	}


	// Given  Parameter: 	ws1, ws2, wp1, wp2, attws, attwp
	// Gotten Parameter: 	bandw, w0, m
	void AnalogFilter::bef_chebyshev(complex<double> pole[], complex<double> zero[])
	{
		ws    = (wp2-wp1)/(ws2-ws1);
		wp    = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2*wp1);
		setm_chebyshev();
		chebyshev(pole, zero);
		a_lpftobef(pole, zero);
	}


	// Given  parameter: 	ws1, ws2, wp1, wp2, attws, attwp
	// Gotten parameter: 	bandw, w0, m
	void AnalogFilter::bef_inv_chebyshev(complex<double> pole[], complex<double> zero[])
	{
		ws    = (wp2-wp1)/(ws2-ws1);
		wp    = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2*wp1);
		setm_inv_chebyshev();
		inv_chebyshev(pole, zero);
		a_lpftobef(pole, zero);
	}


	// Given  parameter: ws1, ws2, wp1, wp2, attws, attwp
	// Gotten parameter: bandw, w0, m
	void AnalogFilter::bef_elliptic(complex<double> pole[], complex<double> zero[])
	{
		ws	   = (wp2 - wp1)/(ws2 - ws1);
		wp	   = 1.0;
		bandw = wp2 - wp1;
		w0    = sqrt(wp2*wp1);
		setm_elliptic();
		elliptic(pole, zero);
		a_lpftobef(pole, zero);
	}


	//	Frequency transform from LPF to BEF
	void AnalogFilter::a_lpftobef(complex<double> pole[], complex<double> zero[])
	{
		complex<double> a,b,c,gns,gnb;
		int 	mm, i, j, ml;
		double 	r, g, bw, ww;
		j         = polem - zerom;
		bw        = bandw;
		ww        = w0;gnb.real(1.0);
		gnb.imag(0.0);
		for(i=0, mm=polem; i<polem; i++, mm++)
		{
			gnb       = gnb * pole[i];
			gnb.real(-gnb.real());
			gnb.imag(-gnb.imag());
			r         = pole[i].real();
			g         = pole[i].imag();
			c.real(r*2.0);
			c.imag(g*2.0);
			a.real(bw*bw - 4.0*ww*ww*r*r + 4.0*ww*ww*g*g);
			a.imag(-8.0*ww*ww*r*g);
			a         = sqrt(a);
			b.real(bw+a.real());
			b.imag(a.imag());
			pole[i]   = b / c;
			b.real(bw - a.real());
			b.imag(-a.imag());
			pole[mm]  = b / c;
		}
		polem = mm;
		gns.real(1.0);
		gns.imag(0.0);
		for(i=0, mm=zerom; i<zerom; i++, mm++)
		{
			gns       = gns * zero[i];
			gns.real(-gns.real());
			gns.imag(-gns.imag());
			r         = zero[i].real();
			g         = zero[i].imag();
			c.real(r*2.0);
			c.imag(g*2.0);
			a.real(bw*bw - 4.0*ww*ww*r*r + 4.0*ww*ww*g*g);
			a.imag(-8.0*ww*ww*r*g);
			a         = sqrt(a);
			b.real(bw+a.real());
			b.imag(a.imag());
			zero[i]   = b / c;
			b.real(bw-a.real());
			b.imag(-a.imag());
			zero[mm]  = b / c;
		}
		zerom  = mm;
		gainb /= gnb.real();
		gains *= gns.real();
		if(j>0)
		{
			mm = zerom;
			ml = mm + j + j;
			for(i=mm; i<ml; i+=2)
			{
				zero[i].real(0.0);
				zero[i+1].real(0.0);
				zero[i].imag(ww);
				zero[i+1].imag(- ww);
			}
			zerom = ml;
		}else if(j<0)
		{
			mm = polem;
			ml = mm-j-j;
			for(i=mm; i<ml; i+=2)
			{
				pole[i].real(0.0);
				pole[i+1].real(0.0);
				pole[i].imag(ww);
				pole[i+1].imag(- ww);
			}
			polem = ml;
		}
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
