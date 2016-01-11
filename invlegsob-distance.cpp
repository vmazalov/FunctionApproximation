#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include<cstdlib>
using namespace std;

vector<double> arclen(const vector<double> &x, const vector<double> &y) {
    vector<double> a;
    double s = 0;
    a.push_back(s);
    for (int i=0; i<x.size()-1; i++) {
        s += sqrt((x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]));
        a.push_back(s);
    }
    return a;
}

vector<double> invariant1(const vector<double> &X, const vector<double> &Y) {
  int n = X.size();
  double Y10 = 0, X10 = 0;
  vector<double> I1;
  I1.push_back(-0.5*X[0]*Y[0]);
  for (int i=1; i<n; i++) {
    double x = X[i];// -X[0];
    double y = Y[i];// -Y[0];
    double dx = X[i]-X[i-1];
    double dy = Y[i]-Y[i-1];
    X10 += y * dx;
    Y10 += x * dy;
    I1.push_back(Y10 - 0.5 * x * y);
  }
  return I1;
}

vector<double> invariant1deg2(const vector<double> &X, const vector<double> &Y) {
  int n = X.size();
  double Y10 = 0, X10 = 0;
  vector<double> I1;
  for (int i=2; i<n-1; i++) {
    double x = X[i];// -X[0];
    double y = Y[i];// -Y[0];
    double dx = (-X[i+2]+8*X[i+1]-8*X[i-1]+X[i-2])/12.0;
    double dy = (-Y[i+2]+8*Y[i+1]-8*Y[i-1]+Y[i-2])/12.0;
    X10 += y * dx;
    Y10 += x * dy;
    I1.push_back(Y10 - 0.5 * x * y);
  }
  return I1;
}


vector<double> invariant2(const vector<double> &X, const vector<double> &Y) {
  int n = X.size();
  double Y11 = 0, Y20 = 0;
  vector<double> I2;
  for (int i=1; i<n; i++) {
    double x = (X[i]+X[i-1]) * 0.5 -X[0];
    double y = (Y[i]+Y[i-1]) * 0.5 -Y[0];
    double dy = Y[i]-Y[i-1];
    Y11 += x * y * dy;
    Y20 += x * x * dy;
    I2.push_back(Y11 * x - 0.5 * Y20 * y - 1.0/6 * x * x * y * y);
  }
  return I2;
}

vector<double> invariant3(const vector<double> &X, const vector<double> &Y) {
  int n = X.size();
  double Y11 = 0, Y20 = 0;
  vector<double> I3;
  for (int i=1; i<n; i++) {
    double x = (X[i]+X[i-1]) * 0.5 - X[0];
    double y = (Y[i]+Y[i-1]) * 0.5 - Y[0];
    double dy = Y[i]-Y[i-1];
    Y11 += x * y * dy;
    Y20 += x * x * dy;
    I3.push_back(Y20 * x + 2 * Y11 * y - 1.0/3 * x*x*x * y - 2.0/3 * x * y*y*y);
  }
  return I3;
}

    
vector<double> moments(const vector<double> &X, const vector<double> &T, int n) {
    vector<double> TP(T);
    vector<double> M;
    int N = T.size();
    double dT = T[N-1] - T[0];
    if (dT == 0) return vector<double>();
    dT = 1/dT; 
    double t = 1.0;
    for (int k=0; k<=n; k++) {
        double s = (X[N-2] + X[N-1]) * TP[N-1] - (X[0] + X[1]) * TP[0];
        for (int i=1; i<N-1; i++)
            s += (X[i-1] - X[i+1]) * TP[i];
        t *= dT;
        M.push_back(s * t / (2.0*k + 2));
        for (int i=0; i<N; i++)
            TP[i] *= T[i];
    }
    return M;    
}

vector<double> momentsToExpansion(const vector<double> &M, const vector<vector<double> > &P) {
    vector<double> e;
    for (int i=0; i<M.size(); i++) {
        double s = 0;
        for (int j=0; j<=i; j++)
            s += M[j] * P[i][j];
        e.push_back(s);
    }
    return e;
}

vector<double> scaleTranslate(const vector<double> &Ex, const vector<double> &Ey, const vector<double> &origin_scale, int start_index) {
    vector<double> e;
    if (start_index == 0) {
        e.push_back((Ex[0] - origin_scale[0]) * origin_scale[2]);
        e.push_back((Ey[0] - origin_scale[1]) * origin_scale[2]);
    }
    for (int i=1; i<Ex.size(); i++) {
        e.push_back(Ex[i] * origin_scale[2]);
        e.push_back(Ey[i] * origin_scale[2]);
    }
    return e;
}

void scaleTranslate(vector<double> &X, vector<double> &Y, const vector<double> &origin_scale) {
	for (int i=0; i<X.size(); i++) {
		X[i] -= origin_scale[0];
		Y[i] -= origin_scale[1];
		X[i] *= origin_scale[2];
		Y[i] *= origin_scale[2];
	}
}

vector<double> getOriginScale(const vector<double> &Ex, const vector<double> &Ey) {
    vector<double> origin_scale(3);
    origin_scale[0] = Ex[0];
    origin_scale[1] = Ey[0];
    double w = 0;
    for (int i=1; i<Ex.size(); i++)
        w += Ex[i]*Ex[i] + Ey[i]*Ey[i];
    if (w <= 0)
        w = 1;
    else
        w = 1/sqrt(w);
    origin_scale[2] = w;
    return origin_scale;
}


class LegendreSobolev {
    vector<double> e0N;
    vector<double> e1N;
    vector<vector<double> > LSPdiffN;
    vector<vector<double> > LSPN;
    int LSPn;
    double mu;
    
public:  
  int getn() { return LSPn; }
  double getmu() { return mu; }

    LegendreSobolev(const char* filename) {
        filebuf fb;
        fb.open(filename,ios::in);
        istream inp(&fb);
        inp >> LSPn;
        inp >> mu;
        double x;
        for (int i=0; i<=LSPn; i++) {
            vector<double>c;
            for (int j=0; j<=i; j++) {
                inp >> x;
                c.push_back(x);
            }
            LSPdiffN.push_back(c);
        }
        for (int i=0; i<=LSPn; i++) {
            inp >> x;
            e0N.push_back(x);
        }
        for (int i=0; i<=LSPn; i++) {
            inp >> x;
            e1N.push_back(x);
        }
        for (int i=0; i<=LSPn; i++) {
            vector<double>c;
            for (int j=0; j<=i; j++) {
                inp >> x;
                c.push_back(x);
            }
            LSPN.push_back(c);
        }
    }    

    vector<double> expansion(vector<double> &X, vector<double> &T, int n) {
        if (n > LSPn) return vector<double>();
        
        vector<double> M(moments(X,T,n));
	if (M.size() == 0) return vector<double>(n+1);
        vector<double> e(momentsToExpansion(M, LSPdiffN));

        int N = X.size();
        for (int i=0; i<e.size(); i++) 
            e[i] += mu * (X[N-1] * e1N[i] - X[0] * e0N[i]);
        return e;
    }

  // Evaluate k-th normalized L-S polynomial at point x
  double evaluate(int k, double x) {
    double r = LSPN[k][k];
    for (int i=k; i>0; i--)
      r = r * x + LSPN[k][i-1];
    return r;
  }

  double evaluate(vector<double> c, double x) {
    double r = 0;
    for (int i=0; i<c.size(); i++)
      r += c[i] * evaluate(i,x);
    return r;
  }
};

void rotate(vector<double> &X, vector<double> &Y, double alpha) {
	double cos_alpha = cos(alpha);
	double sin_alpha = sin(alpha);
	for (int i=0; i<X.size(); i++) {	
		double x = cos_alpha * X[i] + sin_alpha * Y[i];
		double y = -sin_alpha * X[i] + cos_alpha * Y[i];
		X[i] = x;
		Y[i] = y;
	}
}

vector<double> radius(const vector<double> &X, const vector<double> &Y) {
	vector<double> r;
	for (int i=0; i<X.size(); i++) {
	  double x = X[i];
	  double y = Y[i];
	  r.push_back(sqrt(x*x + y*y));
	}
	return r;
}

vector<double> sqradius(const vector<double> &X, const vector<double> &Y) {
	vector<double> r;
	for (int i=0; i<X.size(); i++)
		r.push_back(X[i]*X[i] + Y[i]*Y[i]);
	return r;
}

void read_invariant_coeffs(string filename, vector<vector<vector<double> > > &Gamma) {
  ifstream f(filename.c_str());
  vector<vector<double> > blk;
  while (true) {
    string s;
    getline(f,s);
    if (!f.good()) break;
    if (s=="") continue;
    istringstream inp(s);
    vector<double> ln;
    while (true) {
      double x;
      inp >> x;
      if (!inp.good()) break;
      ln.push_back(x);
    }
    blk.push_back(ln);
    if (blk.size() == ln.size()) {
      Gamma.push_back(blk);
      blk.clear();
    }
  }
}

vector<double> einvar_from_exy(vector<double> &e, vector<vector<vector<double> > > &Gamma) {
  /*
  vector<double> e(2);
  for (int i=0; i<exy.size(); i++)
    e.push_back(exy[i]);
  */

  vector<double> r;
  for (int k=0; k<e.size()/2; k++) {
    double ck=0;
        for (int i=0; i<e.size(); i+=2) 
          for (int j=1; j<e.size(); j+=2) 
	    //        for (int i=0; i<=2*k; i+=2) 
	    //          for (int j=1; j<=2*k; j+=2) 

	ck += e[i] * e[j] * Gamma[k][i/2][j/2];
    r.push_back(ck);
  }
  return r;
}

vector<double> eradius_from_exy(vector<double> &e, vector<vector<vector<double> > > &Gamma) {
  vector<double> r;
  for (int k=0; k<e.size()/2; k++) {
    double ck=0;
    for (int i=0; i<e.size(); i+=2) 
      for (int j=0; j<e.size(); j+=2) 
	ck += (e[i] * e[j] + e[i+1] * e[j+1]) * Gamma[k][i/2][j/2];
    r.push_back(ck);
  }
  return r;
}


double randomDouble()
{
  double(rand()) / RAND_MAX;
} 

int main (int argc, char* argv[]) {
  if (argc != 6) {
    cout << "Usage: legsobjoin file-with-leg-sob-polys(input) xy-file(file with separate strokes, input) order(input) angle(input) leg-sob-coeff-file(output) " << endl;
    cout << "Requires files radius-squared.dat and invariant1.dat in the current directory." << endl;
    cout << "Use maple <leg_sob_matr.mpl to generate them." << endl;
    return -1;
  }

    LegendreSobolev LS(argv[1]);
    filebuf fb;


    vector<vector<vector<double> > > GammaRad;
    read_invariant_coeffs("radius-squared.dat", GammaRad);

    vector<vector<vector<double> > > GammaI1;
    read_invariant_coeffs("invariant1.dat", GammaI1);

    /*
    for (int i=0; i<Gamma.size(); cout << "\n", i++)
      for (int j=0; j<Gamma.size(); cout << "\n", j++)
	for (int k=0; k<Gamma.size(); k++)
	  cout << Gamma[i][j][k] << " " ;
    */

    
      

    fb.open(argv[2],ios::in);
    istream fin(&fb);
    
    int n;
    istringstream ssdim(argv[3]);
    ssdim >> n;

    if (n<1) {
      cout << "Error: order must be positive" << endl;
      return -1;
    }
    if (n>LS.getn()) {
      cout << "Error: order is too large. The known LSP coefficients allow to compute expansions only up to order " << LS.getn() << endl;
      return -1;
    }

    double alpha;
    istringstream ssalpha(argv[4]);
    ssalpha >> alpha;

    filebuf fout;
    fout.open(argv[5],ios::out);  
    ostream out(&fout);
    
    string line;
    while (getline(fin,line)) {
        istringstream inp(line);
        string charname;
        inp >> charname;
        charname = charname.substr(0,charname.length()-1);

        int k;
        inp >> k;
        int L = 0;
        vector<double> X, Y;
        for (int i=0; i<k; i++) {
            getline(fin, line);
            istringstream inp(line);
            string penupdown; inp >> penupdown;
            if (penupdown == "PENUP") continue;
	    double prevx=1e16, prevy=1e16;
            while(!inp.eof()) {
                double x,y;
                inp >> x >> y;
		//		if (prevx != x || prevy != y) {
		  X.push_back(x);
		  Y.push_back(y);
		  prevx = x; prevy = y;
		  //		}
            }
            L ++;
        }

    if (X.size() >= 2) {        
    	    vector<double> al(arclen(X,Y));

            vector<double> Ex(LS.expansion(X,al,n));
            vector<double> Ey(LS.expansion(Y,al,n));           

	        out << charname << ": " << X.size();

    	    for (int i=0; i<X.size(); i++) 
    	      out
    		 << " " << X[i] << " " << LS.evaluate(Ex,al[i]/al[al.size()-1]) << " " << Y[i] << " " << LS.evaluate(Ey,al[i]/al[al.size()-1]);
    	    out << endl;


	} else {
	    out << charname << ":";
	    for (int i=0; i<=n; i++) 
	      out << " 0 0";
	    out << endl;
	}
    }
    
    return 0;
}
