#include "mzUtils.h"

//random collection of useful functions

using namespace std;

namespace mzUtils { 

std::string makeLowerCase(string &s) {
		for (unsigned int i=0; i != s.length(); ++i ) {
			s[i] = std::tolower(s[i]);
		}	
		return s;
}

void split(const string& s, char c, vector<string>& v) {
    string::size_type i = 0;
    string::size_type j = s.find(c);

    while (j != string::npos) {
        v.push_back(s.substr(i, j-i));
        i = ++j;
        j = s.find(c, j);

        if (j == string::npos)
            v.push_back(s.substr(i, s.length( )));
    }
    if ( v.size() == 0) v.push_back(s);
}

char *mystrcasestr(const char *s1, const char *s2) {
	 const char *s = s1;
	 const char *p = s2;
	do {
		if (!*p) {
			return (char *) s1;;
		}
		if ((*p == *s)
			|| (tolower(*((unsigned char *)p)) == tolower(*((unsigned char *)s)))
			) {
			++p;
			++s;
		} else {
			p = s2;
			if (!*s) {
				return NULL;
			}
			s = ++s1;
		}
	} while (1);
}

/* ========================================================================= */
string substituteInQuotedString(const string& s, const string& chars, const string& substitutions ) {
 string result;
 for (string::size_type pos = 0; pos < s.size(); ++pos) {
     char c = s[pos];
     string::size_type subst_pos = chars.find_first_of(c);
     if (subst_pos != string::npos) c = substitutions[subst_pos]; 
     result += c;
 } 
  return result;
}

/**
 * @brief smoothAverage
 * @param y
 * @param s
 * @param smoothWindowLen
 * @param ly
 *
 * @deprecated
 * Issue 42: Not thread safe! Do not use!
 */
void smoothAverage(float *y, float* s, int smoothWindowLen, int ly) {
    if (smoothWindowLen == 0 ) return;

    float *x = new float[smoothWindowLen];
    for(int i=0; i< smoothWindowLen; i++ ){
        x[i] = 1.0/smoothWindowLen;
    }

    conv(smoothWindowLen,-smoothWindowLen/2, x, ly, 0, y, ly, 0, s);

    delete[] x;
}

/**
  *
  * @deprecated
  *
  * @warning
  * Not thread safe!
  */
void conv (int lx, int ifx, float *x, int ly, int ify, float *y, int lz, int ifz, float *z) /*****************************************************************************
Compute z = x convolved with y; i.e.,

           ifx+lx-1
    z[i] =   sum    x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
            j=ifx
******************************************************************************
Input:
lx		length of x array
ifx		sample index of first x
x		array[lx] to be convolved with y
ly		length of y array
ify		sample index of first y
y		array[ly] with which x is to be convolved
lz		length of z array
ifz		sample index of first z

Output:
z		array[lz] containing x convolved with y
******************************************************************************
Notes:
The x samples are contained in x[0], x[1], ..., x[lx-1]; likewise for
the y and z samples.  The sample indices of the first x, y, and z values
determine the location of the origin for each array.  For example, if
z is to be a weighted average of the nearest 5 samples of y, one might
use
	...
	x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
	conv(5,-2,x,lx,0,y,ly,0,z);
	...
In this example, the filter x is symmetric, with index of first sample = -2.

This function is optimized for architectures that can simultaneously perform
a multiply, add, and one load from memory; e.g., the IBM RISC System/6000.
Because, for each value of i, it accumulates the convolution sum z[i] in a
scalar, this function is not likely to be optimal for vector architectures.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/23/91
*****************************************************************************/
{
	int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,i,j,jlow,jhigh;
	float sum;

	x -= ifx;  y -= ify;  z -= ifz;
	for (i=ifz; i<=ilz; ++i) {
		jlow = i-ily;  if (jlow<ifx) jlow = ifx;
		jhigh = i-ify;  if (jhigh>ilx) jhigh = ilx;
		for (j=jlow,sum=0.0; j<=jhigh; ++j)
			sum += x[j]*y[i-j];
		z[i] = sum;
	}
}

/**
 * @brief gaussian1d_smoothing
 * @param ns
 * @param nsr
 * @param data
 *
 * @warning
 * Not thread-safe!
 *
 * @deprecated
 */
void gaussian1d_smoothing (int ns, int nsr, float *data) 
{
//Subroutine to apply a one-dimensional gaussian smoothing

/******************************************************************************
Input:
ns		number of samples in the input data
nsr		width (in samples) of the gaussian for which
		amplitude > 0.5*max amplitude
data		1-D array[ns] of data to smooth

Output:
data		1-D array[ns] of smoothed data
******************************************************************************/
	int is;				/* loop counter */
	float sum=0.0;
	float fcut;
	float r;
	float fcutr=1.0/nsr;
	static int n;
	static int mean;
	static float fcutl;
	static float s[1000];		/* smoothing filter array */
	float *temp;			/* temporary array */

	/* save input fcut */
	fcut=fcutr;

	/* don't smooth if nsr equal to zero */
	if (nsr==0 || ns<=1) return;

	/* if halfwidth more than 100 samples, truncate */
	if (nsr>100) fcut=1.0/100;

	/* allocate space */
	temp=new float[ns];


	/* initialize smoothing function if not the same as the last one used */
	if (fcut != fcutl) {
		fcutl=fcut;

		/* set span of 3, at width of 1.5*exp(-PI*1.5**2)=1/1174 */
		n=(int) (3.0/fcut+0.5);
		n=2*n/2+1;		/* make it odd for symmetry */

		/* mean is the index of the zero in the smoothing wavelet */
		mean=n/2;

		/* s(n) is the smoothing gaussian */
		for (is=1; is<=n; is++) {
			r=is-mean-1;
			r= -r*r*fcut*fcut*3.141;
			s[is-1]=exp(r);
		}

		/* normalize to unit area, will preserve DC frequency at full
		amplitude. Frequency at fcut will be half amplitude */
		for (is=0; is<n; is++) sum +=s[is];
		for (is=0; is<n; is++) s[is] /=sum;
		//cerr << "new gaussian " << nsr << endl;
	}

	/* convolve by gaussian into buffer */
	if (1.01/fcutr>(float)ns) {

		/* replace drastic smoothing by averaging */
		sum=0.0;
		for (is=0; is<ns; is++) sum +=data[is];
		sum /=ns;
		for (is=0; is<ns; is++) data[is]=sum;

	} else {

		/* convolve with gaussian */
		conv (n, -mean, s, ns, -mean, data, ns, -mean, temp);

		/* copy filtered data back to output array */
		for (is=0; is<ns; is++) data[is]=temp[is];
	}

	/* free allocated space */
	delete[] temp;
}

float median(vector <float> y) {
	if (y.empty() ) return(0.0);
	if (y.size() == 1 ) return(y[0]);
	if (y.size() == 2 ) { return(y[0]+y[1])/2; }

	unsigned int n = y.size();
	std::sort(y.begin(), y.end());

	unsigned int i = n / 2;
	float median = 0.0;

	if (n == i * 2)
	    median = (y[i-1] + y[i]) / 2.;
	else
	    median = y[i];

	return (median);
}

float median(float* y, int n) {
	if (n == 0) return(0.0);
	if (n == 1) return(y[0]);
	if (n == 2) return((y[0]+y[1])/2);

	float* tmpy = new float[n];
	memcpy(tmpy,y,n*sizeof(float));
	sort(tmpy,tmpy+n);

	int i = n / 2;
	float median = 0.0;

	if (n == i * 2)
	    median = (tmpy[i-1] + tmpy[i]) / 2.;
	else
	    median = tmpy[i];

	delete[] tmpy;
	//for(int i=0; i < n; i++ ) { cerr << tmpy[i] << " "; }
	//cerr << "median=" << median << endl;
	return (median);
}

/*
  * The following code is public domain.
  * Algorithm by Torben Mogensen, implementation by N. Devillard.
  * This code in public domain.
  *  In place median, no need to allocate memory
*/
float torben_median(const vector<float> &m) {
     int         i, less, greater, equal;
     float min, max, guess, maxltguess, mingtguess;
	 int n = m.size();
	 if (n == 0 ) return 0;
	 if (n == 1 ) return m[0]; 
     min = max = m[0] ;

     for (i=1 ; i<n ; i++) {
         if (m[i]<min) min=m[i];
         if (m[i]>max) max=m[i];
     }
     while (1) {
         guess = (min+max)/2;
         less = 0; greater = 0; equal = 0;
         maxltguess = min ;
         mingtguess = max ;
         for (i=0; i<n; i++) {
             if (m[i]<guess) {
                 less++;
                 if (m[i]>maxltguess) maxltguess = m[i] ;
             } else if (m[i]>guess) {
                 greater++;
                 if (m[i]<mingtguess) mingtguess = m[i] ;
             } else equal++;
         }
         if (less <= (n+1)/2 && greater <= (n+1)/2) break ;
         else if (less>greater) max = maxltguess ;
         else min = mingtguess;
     }
     if (less >= (n+1)/2) return maxltguess;
     else if (less+equal >= (n+1)/2) return guess;
     else return mingtguess;
}

/*
  * Algorithm from N. Wirth’s book, implementation by N. Devillard.
  * This code in public domain.
  */
#define ELEM_SWAP(a,b) {  float t=(a);(a)=(b);(b)=t; }
/*---------------------------------------------------------------------------
    Function :   kth_smallest()
             :   array of elements, # of elements in the array, rank k
    In
             :   one element
    Out
    Job      :   find the kth smallest element in the array
    Notice   :   use the median() macro defined below to get the median.
                 Reference:
                   Author:  Wirth, Niklaus
                    Title:  Algorithms + data structures = programs
                Publisher:  Englewood Cliffs: Prentice-Hall, 1976
     Physical description: 366 p.
                   Series:  Prentice-Hall Series in Automatic Computation
  ---------------------------------------------------------------------------*/
float kth_smallest(float a[], int n, int k)
{
      int i,j,l,m ;
      float x ;
     l=0 ; m=n-1 ;
     while (l<m) {
         x=a[k] ;
         i=l ;
         j=m ;
         do {
             while (a[i]<x) i++ ;
             while (x<a[j]) j-- ;
             if (i<=j) {
                 ELEM_SWAP(a[i],a[j]) ;
                 i++ ; j-- ;
             }
         } while (i<=j) ;
         if (j<k) l=i ;
         if (k<i) m=j ;
     }
     return a[k] ;
}
#define wirth_median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))
/* string2integer */
int string2integer(const std::string& s){
	std::istringstream i(s);
	int x = 0;
	i >> x;
	return x;
}

/* string2float */
float string2float(const std::string& s){
	std::istringstream i(s);
	float x = 0;
	i >> x; 
	return x;
}

/* string2integer */
string integer2string(int x){
	std::stringstream i;
	string s;
	i << x;
	i >> s;
	return s;
}

string float2string(float f, int p) {
	std::stringstream ss;
	ss << setprecision(p) << f; string str; ss >> str;
	return(str);
}
		
float ppmDist(const float mz1, const float mz2) { 
		return ( abs((mz2-mz1)/(mz1/1e6)) ); 
}

double ppmDist(const double mz1, const double mz2) { 
		return ( abs((mz2-mz1)/(mz1/1e6)) ); 
}



float ppmround(const float mz1, const float resolution) {
    //resolution parameter =10  -> one digit after decimal point,
    //                      100 -> two digits after decimal point
    //                      etc..

    return( round(mz1*resolution)/resolution);
}

bool withinXppm( float mz1, float mz2, int ppmWindow ) { 
	if ( mz2 > (mz1 - (mz1/1e6)*ppmWindow) && mz2 < (mz1 + (mz1/1e6)*ppmWindow) ) return(true); 
	else return(false);
}

vector<float> quantileDistribution( vector<float> y ) {
    int ysize = y.size();
    std::sort(y.begin(), y.end());	//sort y
    vector<float> quantiles(101,0);
    for (int i=0; i < 101; i++ ) {
        int pos = (float) i/100 * ysize;
        if (pos < ysize) quantiles[i] = y[pos];
    }
    return(quantiles);
}


float ttest(StatisticsVector<float>& groupA, StatisticsVector<float>& groupB ) {
		int n1 = groupA.size();
		int n2 = groupB.size();
		if ( n1 == 0 && n2 == 0) return 0;       //both empty.. no different
		if ( n1 == 0 || n2 == 0) return 1000;    //one is empty inf difference

		float meanA = groupA.mean();
		float meanB = groupB.mean();

		float stdA  = groupA.stddev(meanA);
		float stdB  = groupB.stddev(meanB);
		if (stdA == 0 ) stdA=1.0;
		if (stdB == 0 ) stdB=1.0;

		float t_test = (meanA-meanB)/sqrt(((stdA*stdA)/n1)+((stdB*stdB)/n2));
                return t_test;
            }


int countBelow(vector<float>& y, float ymax) { 
	vector<float>::iterator itr = lower_bound(y.begin(), y.end(), ymax);
	int lb = itr-y.begin();
	return lb;
/*
	int count=0;
		for (int i=0; i < y.size(); i++) if (y[i] < ymax) count++; 
		return(count);
		*/
}

bool fileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;
  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    blnReturn = true;
  } else {
    blnReturn = false;
  }
  return(blnReturn);
}


int createDir(const char* path) {
  if (isDir(path)) return 0;
   cerr << "Creating path: " << path << endl;
   mode_t old_mask = umask(0);
#ifdef MINGW
   int retval = mkdir(path);
#else
   int retval = mkdir(path, 0771);
#endif
   umask(old_mask);
   return retval;
}

int isFile(const char* path) {
    struct stat sbuf;
    int retval = stat(path, &sbuf);
    return (!retval && (sbuf.st_mode & S_IFREG));
}

int isDir(const char* path) {
    struct stat sbuf;
    int retval = stat(path, &sbuf);
    return (!retval && (sbuf.st_mode & S_IFDIR));
}


float correlation(const vector<float>&x, const vector<float>&y) {
    int n = x.size();
    if (n == 0) return 0;

    if (x == y) return 1;

    float sumx = 0; 		//
    float sumy = 0;
    float sumxy =0;
    float x2 = 0;
    float y2 = 0;

    for (int i = 0; i < n; i++) {
        sumx += x[i];
        sumy += y[i];
        sumxy += x[i]*y[i];
        x2 += x[i]*x[i];
        y2 += y[i]*y[i];
    }

	float var1 = x2-(sumx*sumx)/n;
	float var2 = y2-(sumy*sumy)/n;
	if ( var1 == 0 || var2 == 0 ) return 0;
    return (sumxy -( sumx*sumy)/n) / sqrt((x2-(sumx*sumx)/n)*(y2-(sumy*sumy)/n));
}

tuple<double, double, string> parseMspFragLine(string line){

    int space1 = -1;
    int space2 = -1;
    bool isHasNonWhiteSpace = false;

    for (unsigned int i = 0; i < line.size(); i++) {
        char c = line[i];

        bool isNextWhiteSpace = false;
        if (i < line.size()-1) {
            char c2 = line[i+1];
            isNextWhiteSpace = std::isspace(c2);
        }

        if (!std::isspace(c)) {
            isHasNonWhiteSpace = true;
        } else if (std::isspace(c) && !isNextWhiteSpace) {
            if (space1 == -1){
                space1 = static_cast<int>(i);
            } else if (space2 == -1) {
                space2 = static_cast<int>(i);
                break;
            }
        }
    }

    string fragLabel("");
    double fragMz(-1);
    double fragIntensity(-1);

    if (isHasNonWhiteSpace && space1 != -1) { //avoid exception on bad formatting
        fragMz = stod(line.substr(0, static_cast<unsigned long>(space1)));

        if (space2 == -1) {
            //no label
            fragIntensity = stod(line.substr(static_cast<unsigned long>(space1+1), (line.size()-static_cast<unsigned long>(space1))));
        } else {
            fragIntensity = stod(line.substr(static_cast<unsigned long>(space1+1), static_cast<unsigned long>(space2-space1)));
            fragLabel = line.substr(static_cast<unsigned long>(space2+1), (line.size()-static_cast<unsigned long>(space2)));
        }
    } else {
        //cerr << "Fragment label \"" << line << "\"" << " could not be parsed! Returning (0, 0, \"\")" << endl;
        //bad formatting
    }

    tuple<double, double, string> fragInfo = tuple<double, double, string>(fragMz, fragIntensity, fragLabel);
    return fragInfo;
}

/**
 * peak fitting function
 *
 * @brief
 * for cases where a fit cannot be assessed, return the default values as
 *
 * @param ycoord: reference to intensity vector
 * @param default_sigma: sigma value to be returned if no fit can be determined
 * @param default_gauss_fit_R2: R2 gaussian fit value to be returned if no fit can be determined
 *
 * @return pair<float, float>: <sigma, gauss_R2_fit> as determined by this algorithm
*/
pair<float, float> gaussFit(const vector<float>&ycoord, float default_sigma, float default_gauss_fit_R2) {

        float s = 20;
		float min_s = 0;
        float minR = FLT_MAX;

		//find best fit
        if (ycoord.size()<3) return make_pair(default_sigma, default_gauss_fit_R2);
		vector<float>yobs=ycoord;
		int ysize=yobs.size();
		int midpoint  = int(ysize/2);
		//find maximum point ( assuming it somewhere around midpoint of the yobs);
		float ymax = max(max(yobs[midpoint], yobs[midpoint-1]),yobs[midpoint+1]);
		float ymin = min( yobs[0], yobs[ysize-1]);

		//initialize x vector, values centered around 0, forxample  -2, -1, 0, 1, 2
		int xinit = int(ysize/2)*-1;
		vector<float>x(ysize,0);
		int greaterZeroCount=0;

		for(int i=0; i<ysize; i++ ) {
				x[i] = xinit+i; 
				if ( yobs[i] > ymin ) greaterZeroCount++;
				yobs[i] = (yobs[i]-ymin)/(ymax-ymin); 
				if(yobs[i]<0) yobs[i]=0; 
		}

		/*
		cerr << "fitting yobs:" << endl; for(int i=0; i < ysize; i++ ) 
				cerr << setprecision(2) << yobs[i] << ", "; cerr << endl;
		*/
	
        bool converged = false; 
		int ittr = 0; 

        if (greaterZeroCount <= 3 ) return make_pair(default_sigma, default_gauss_fit_R2);
		while (!converged ) {
				if ( ittr++ > 20 ) break;
                float Rsqr=0;
				for(int i=0; i < ysize; i++ )  { Rsqr += POW2(exp(-0.5*POW2(x[i]/s)) - yobs[i]); }
         //       cerr << "\t\ts=" << s << " Rsqr=" << Rsqr << endl;
				if ( Rsqr < minR || ittr == 0 ) { minR = Rsqr; min_s = s; }
                else if ( Rsqr >= minR ) break;
				s /= 1.25;
		}

//		*sigma = min_s;
//		*R2 = minR/(ysize*ysize);	//corrected R2
//        //cerr << "fit() s=" << *sigma << " R2=" << *R2 << endl;

        return make_pair(min_s, minR/(ysize*ysize));
	}


inline unsigned long factorial(int n) { 
		long p=1; while(n>1) p*=n--; return p; }

		/*
long nchoosek(int n, int k) { 
	if(k==n) return 1;
	if(k==0) return 1;
	if(k>n) return 0;
	return (factorial(n)/(factorial(n-k)*factorial(k))); 
}
*/
long long nchoosek(int n, int k) {

		int n_k = n - k;

		if (k < n_k)
		{
				k = n_k;
				n_k = n - k;
		}

        long long  nchsk = 1;
		for ( int i = 1; i <= n_k; i++)
		{
				nchsk *= (++k);
				nchsk /= i;
		}
		return nchsk;
}

string cleanFilename(const string& filename) { 

		string outstring=filename;
		std::string::size_type pos =outstring.find_last_of("/");

		if (pos != std::string::npos) { 
			outstring=outstring.substr(pos+1, outstring.length());
		}

		pos=outstring.find_last_of("\\");
		if (pos != std::string::npos) { 
				outstring=outstring.substr(pos+1, outstring.length());
		}

		pos=outstring.find_last_of(".");
		if (pos != std::string::npos) { 
			outstring=outstring.substr(0,pos);
		}
		return outstring;
}

std::vector<double> naturalAbundanceCorrection(int nC, std::vector<double>& M) {
    //inputs
    int n = nC;		//number of carbonms
    int nr_13C=n-1; //number of labeled carbons

    //uncorrected values
    vector<double>C(n,0); // output vector with corrected values

    for(int k=0; k < nr_13C-1; k++ ) {
        double contamination=0;

        for(int i=0; i<k; i++ ) {
            contamination += pow( 0.011,k-i) * pow(0.989, n-k) * nchoosek(n-i,k-i) * C[i];
        }

        C[k] = (M[k]-contamination) / pow(0.989,n-k);
        if(C[k] < 1e-4) C[k]=0;
    }

    return C;
}


const long double 	Pi = 3.1415926535897932384626433832795028841968;

double beta(double x, double y) { 
	
	if (x >= 2 && y >= 2 ) { //approximation
		return sqrt(2*Pi)*pow(x,x-0.5)*pow(y,y-0.5)/pow(x+y,(x+y-0.5));
	}

	//integral form
	double dt=0.01;
	double sum=0;
	for(double t=0.0001; t<=0.9999; t+=dt) {
			sum += pow(t,(x-1))*pow((1-t),(y-1))*dt;
	}
	return sum;
}

double gamma(double z) { 
	//integral form
	double dt=0.0001;
	double sum=0;
	for(double t=0.000001; t<=10; t+=dt) {
		sum += pow(t,z-1)*exp(-t)*dt;
	}
	return sum;
}

double betaPDF(double x, double a, double b) { 
	return pow(x,a-1)*pow(1-x,b-1)/beta(a,b);
}

double pertPDF(double x, double min, double mode, double max ) {
	double a = 6*(mode-min)/(max-min);
	double b = 6*(max-mode)/(max-min);
	return pow(x-min,a-1)*pow(max-x,b-1)/(beta(a,b)*pow(max-min,a+b-1));
}


void tridiagonal ( int n, float *c, float *a, float *b, float *r )

{
     int i;

     for ( i = 0; i < n-1; i++ ) {
         b[i] /= a[i];
         a[i+1] -= c[i]*b[i];
     }

     r[0] /= a[0];
     for ( i = 1; i < n; i++ )
         r[i] = ( r[i] - c[i-1] * r[i-1] ) / a[i];

     for ( i = n-2; i >= 0; i-- )
         r[i] -= r[i+1] * b[i];
}


void cubic_nak ( int n, float *x, float *f, float *b, float *c, float *d )

/*
     PURPOSE:
          determine the coefficients for the 'not-a-knot'
          cubic spline for a given set of data


     CALLING SEQUENCE:
          cubic_nak ( n, x, f, b, c, d );


     INPUTS:
          n		number of interpolating points
          x		array containing interpolating points
          f		array containing function values to
                        be interpolated;  f[i] is the function
                        value corresponding to x[i]
          b		array of size at least n; contents will
                        be overwritten
          c		array of size at least n; contents will
                        be overwritten
          d		array of size at least n; contents will
                        be overwritten


     OUTPUTS:
          b		coefficients of linear terms in cubic
                        spline
          c		coefficients of quadratic terms in
                        cubic spline
          d		coefficients of cubic terms in cubic
                        spline

     REMARK:
          remember that the constant terms in the cubic spline
          are given by the function values being interpolated;
          i.e., the contents of the f array are the constant
          terms

          to evaluate the cubic spline, use the routine
          'spline_eval'
*/

{
     float *h,
            *dl,
            *dd,
            *du;
     int i;

     h  = new float [n];
     dl = new float [n];
     dd = new float [n];
     du = new float [n];

     for ( i = 0; i < n-1; i++ )
         h[i] = x[i+1] - x[i];
     for ( i = 0; i < n-3; i++ )
         dl[i] = du[i] = h[i+1];

     for ( i = 0; i < n-2; i++ ) {
         dd[i] = 2.0 * ( h[i] + h[i+1] );
         c[i]  = ( 3.0 / h[i+1] ) * ( f[i+2] - f[i+1] ) -
                 ( 3.0 / h[i] ) * ( f[i+1] - f[i] );
     }
     dd[0] += ( h[0] + h[0]*h[0] / h[1] );
     dd[n-3] += ( h[n-2] + h[n-2]*h[n-2] / h[n-3] );
     du[0] -= ( h[0]*h[0] / h[1] );
     dl[n-4] -= ( h[n-2]*h[n-2] / h[n-3] );

     tridiagonal ( n-2, dl, dd, du, c );

     for ( i = n-3; i >= 0; i-- )
         c[i+1] = c[i];
     c[0] = ( 1.0 + h[0] / h[1] ) * c[1] - h[0] / h[1] * c[2];
     c[n-1] = ( 1.0 + h[n-2] / h[n-3] ) * c[n-2] - h[n-2] / h[n-3] * c[n-3];
     for ( i = 0; i < n-1; i++ ) {
         d[i] = ( c[i+1] - c[i] ) / ( 3.0 * h[i] );
         b[i] = ( f[i+1] - f[i] ) / h[i] - h[i] * ( c[i+1] + 2.0*c[i] ) / 3.0;
     }

     delete [] h;
     delete [] du;
     delete [] dd;
     delete [] dl;
 }

float spline_eval ( int n, float *x, float *f, float *b, float *c,
                     float *d, float t )

/*
     PURPOSE:
          evaluate a cubic spline at a single value of
          the independent variable given the coefficients of
          the cubic spline interpolant (obtained from
          'cubic_nak' or 'cubic_clamped')


     CALLING SEQUENCE:
          y = spline_eval ( n, x, f, b, c, d, t );
          spline_eval ( n, x, f, b, c, d, t );


     INPUTS:
          n		number of interpolating points
          x		array containing interpolating points
          f		array containing the constant terms from
                        the cubic spline (obtained from 'cubic_nak'
                        or 'cubic_clamped')
          b		array containing the coefficients of the
                        linear terms from the cubic spline
                        (obtained from 'cubic_nak' or 'cubic_clamped')
          c		array containing the coefficients of the
                        quadratic terms from the cubic spline
                        (obtained from 'cubic_nak' or 'cubic_clamped')
          d		array containing the coefficients of the
                        cubic terms from the cubic spline
                        (obtained from 'cubic_nak' or 'cubic_clamped')
          t		value of independent variable at which
                        the interpolating polynomial is to be
                        evaluated


     OUTPUTS:
          y		value of cubic spline at the specified
                        value of the independent variable
*/

{
     int i=1;
     int found=0;

     while ( !found && ( i < n-1 ) ) {
           if ( t < x[i] )
              found = 1;
           else
              i++;
     }

     t = f[i-1] + ( t - x[i-1] ) * ( b[i-1] + ( t - x[i-1] ) * ( c[i-1] + (t - x[i-1] ) * d[i-1] ) );
     return ( t );
}


/** Decompress an STL string using zlib and return the original data. */
std::string decompress_string(const std::string& str)
{
    std::string outstring;

#ifdef ZLIB
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (inflateInit(&zs) != Z_OK)
        throw(std::runtime_error("inflateInit failed while decompressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();

    int ret;
    char outbuffer[32768];

    // get the decompressed bytes blockwise using repeated calls to inflate
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = inflate(&zs, Z_NO_FLUSH);

        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }

    } while (ret == Z_OK);

    inflateEnd(&zs);

   if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib decompression: (" << ret << ") "
            << zs.msg;
       // throw(std::runtime_error(oss.str()));
    }

#endif
    return outstring;
}

bool ends_with(std::string const & value, std::string const & ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}


// for correlattion just reverse the order of sequence y;
float  crossCorrelationZ(vector<float>&xvector, vector<float>& yvector, float offset=.45)
{
	int lenx = xvector.size();
	int leny = yvector.size();
	int lenz= lenx+leny-1;

	//cerr << "LEN: " << lenx << "\t" << leny << endl;
	float meanx=0; float meany=0;
	for (float xi: xvector) meanx+=xi;
	for (float yi: yvector) meany+=yi;
	meanx /= lenx;
	meany /= leny;

	float obs_corr=0;
	float exp_corr=0;
	float exp_n=0; 
	
	for (int i=offset*lenz;i<lenz*(1-offset);i++) {
		int yp=(leny-i-1);
		if (yp < 0) yp=0;
		if (yp >= leny) yp = (leny-1);

		int xp=(i-leny+1);
		if (xp < 0) xp = 0;
		if (xp >= lenx) xp = (lenx-1);

		bool isShifted=true;
		if (xp == 0 and yp == 0) isShifted=false;

//		cerr << i << "\t" << xp << ":" << yp << endl;
		float xy=0; float xx=0; float yy=0;
		while(xp<lenx & yp<leny) {
			float x = xvector[xp]-meanx;
			float y = yvector[yp]-meany;
			xy +=  x*y;
			xx +=  x*x;
			yy +=  y*y;
			xp++;
			yp++;
		}
		float r = xy / sqrt(xx*yy);
		if(isShifted) { exp_n++; exp_corr+=r; } else { obs_corr=r; }
	}

	//cerr << exp_n << endl;
	if ( exp_n > 0 ) {
		exp_corr /= exp_n;
		//float exp_corr = E.mean();
		//float exp_std  = E.stddev(exp_corr);
		//if (exp_std <= 0.05)  exp_std = 0.05;
		float score = (obs_corr - exp_corr)/0.05;
		if (score > 100) score = 100;
		return(score);
	} else {
		return(0);
	}
}

long mzToIntKey(const double mz, const long multFactor){
    return static_cast<long>(round(mz*multFactor));
}

double intKeyToMz(const long intKey, const long multFactor){
    return static_cast<double>(intKey) / static_cast<double>(multFactor);
}

//Issue 270: implement match parsimony
vector<vector<int>> simpleParsimonyReducer(vector<vector<int>> input){

    if (input.empty() || input.size() == 1) {
        return input;
    }

    sort(input.begin(), input.end(), [](const vector<int>& lhs, const vector<int>& rhs){
        return lhs.size() < rhs.size();
    });

    unordered_set<unsigned int> idsToSkip{};

    for (unsigned int i = 0; i < input.size(); i++) {

        vector<int> ith = input[i];
        unordered_set<int> ithIds;

        for (unsigned int j = i+1; j < input.size(); j++) {

            vector<int> jth = input[j];

            //removal requires same element in group with more members
            if (jth.size() > ith.size()) {
                vector<int> intersection;
                std::set_intersection(ith.begin(), ith.end(), jth.begin(), jth.end(), back_inserter(intersection));
                for (auto x : intersection) {
                    ithIds.insert(x);
                }
            }

            //skip this id
            if (ithIds.size() == ith.size()) {
                idsToSkip.insert(i);
                break;
            }

        }
    }

    unsigned int counter = 0;
    vector<vector<int>> output(input.size()-idsToSkip.size());
    for (unsigned int i = 0; i < input.size(); i++) {
        if (idsToSkip.find(i) == idsToSkip.end()) {
            output[counter] = input[i];
            counter++;
        }
    }

    return output;
}

vector<string> getMzSampleFilesFromDirectory(const char* path){

    vector<string> fileNames;

    #ifdef _WIN32
    cerr << "TODO: mzUtils::getMzSampleFilesFromDirectory() not available on windows!" << endl;
    abort();
    #else
    glob_t glob_result;
    string mzMLglob = std::string(path) + "/*.mzML*";
    string mzXMLglob = std::string(path) + "/*.mzXML*";
    vector<string> globpatterns { mzMLglob, mzXMLglob };
    for( string pattern : globpatterns ) {
        glob(pattern.c_str(), GLOB_TILDE,nullptr,&glob_result);
        for(unsigned int i=0; i<glob_result.gl_pathc; ++i) fileNames.push_back(glob_result.gl_pathv[i]);
    }
    #endif

    return fileNames;

    /*
    vector<string> fileNames;
    QDirIterator dirIterator(QString(path), QStringList({"*.mzML", "*.mzXML"}));
    while (dirIterator.hasNext()){
        fileNames.push_back(dirIterator.next().toStdString());
    }
    return fileNames;
    */
}

unordered_map<string, string> decodeParameterMap(string encodedParams, string delimiter){

    unordered_map<string, string> decodedMap = {};

    if (encodedParams.size() < 3) return decodedMap;

    //remove starting/ending brackets, if they are included
    if (encodedParams[0] == '{'){
        encodedParams.erase(0, 1);
    }
    if (encodedParams[encodedParams.size()-1] == '}') {
        encodedParams.erase(encodedParams.size()-1);
    }

    unsigned long posPrevious = 0;
    unsigned long posCurrent = 0;

    while ((posCurrent = encodedParams.find(delimiter, posPrevious)) != string::npos) {

        string encodedParam = encodedParams.substr(posPrevious, posCurrent-posPrevious);
        posPrevious = posCurrent + delimiter.length();

        unsigned long equalCoord = encodedParam.find("=");

        string paramKey = encodedParam.substr(0, equalCoord);
        string paramVal = encodedParam.substr(equalCoord + 1, encodedParam.size());

        decodedMap.insert(make_pair(paramKey, paramVal));
    }

    return decodedMap;

} //decodeParameterMap()

vector<string> decodeParameterVector(string encodedParams, string delimiter) {
    vector<string> decodedVector{};

    if (encodedParams.size() < 3) return decodedVector;

    //remove starting/ending brackets, if they are included
    if (encodedParams[0] == '{'){
        encodedParams.erase(0, 1);
    }
    if (encodedParams[encodedParams.size()-1] == '}') {
        encodedParams.erase(encodedParams.size()-1);
    }

    unsigned long posPrevious = 0;
    unsigned long posCurrent = 0;

    while ((posCurrent = encodedParams.find(delimiter, posPrevious)) != string::npos) {

        string encodedParam = encodedParams.substr(posPrevious, posCurrent-posPrevious);
        posPrevious = posCurrent + delimiter.length();

        decodedVector.push_back(encodedParam);

    }

    return decodedVector;

} //decodeParameterVector()

} //namespace end

