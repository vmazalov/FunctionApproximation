#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
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

vector<vector<vector<double> > > Gamma;

void read_Gamma(char* filename) {
  ifstream f(filename);
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

vector<double> e1_from_exy(vector<double> &e) {
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

int main (int argc, char* argv[]) {
  if (argc != 6) {
    cout << "Usage: legsobjoin file-with-leg-sob-polys(input) Gamma-file xy-file(file with separate strokes, input) order(input) leg-sob-coeff-file(output) " << endl;
    return -1;
  }

    LegendreSobolev LS(argv[1]);
    filebuf fb;

    read_Gamma(argv[2]);

    /*
    for (int i=0; i<Gamma.size(); cout << "\n", i++)
      for (int j=0; j<Gamma.size(); cout << "\n", j++)
	for (int k=0; k<Gamma.size(); k++)
	  cout << Gamma[i][j][k] << " " ;
    */

    
      

    fb.open(argv[3],ios::in);
    istream fin(&fb);
    
    int n;
    istringstream ssdim(argv[4]);
    ssdim >> n;

    if (n<1) {
      cout << "Error: order must be positive" << endl;
      return -1;
    }
    if (n>LS.getn()) {
      cout << "Error: order is too large. The known LSP coefficients allow to compute expansions only up to order " << LS.getn() << endl;
      return -1;
    }

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
		if (prevx != x && prevy != y) {
		  X.push_back(x);
		  Y.push_back(y);
		  prevx = x; prevy = y;
		}
            }
            L ++;
        }

	vector<double> exy, er, e1, e2, e3, e1quick, e1deg2;
        if (X.size() >= 2) {        
	  //  rotate(X,Y,0.923456);
	    	  vector<double> al(arclen(X,Y));

		  /** Test: replace arc length by time
				 vector<double> al;
			  for (int i=0; i<X.size(); i++)
			    al.push_back(1.0 * i / (X.size()-1));
		  //	  	  for (int i=0; i<al.size(); i++) 
	  //	  	    cout << al[i] << " ";
	  //	  	  cout << endl;
	  */


            vector<double> Ex(LS.expansion(X,al,n));
	    /*
	    cout << Ex.size() << endl;
	    for (int i=0; i<X.size(); i++) 
	      cout << " " << X[i] << " " << LS.evaluate(Ex,i*1.0/(X.size()-1)) << "   ";
	    cout << endl;
	    */
            vector<double> Ey(LS.expansion(Y,al,n));           

            vector<double> origin_scale(getOriginScale(Ex,Ey));

	    /*
	    scaleTranslate(X,Y,origin_scale);
	    */



	    Ex[0] = 0; Ey[0] = 0;
	    for (int i=1; i<Ex.size(); i++) 
	      Ex[i]*=origin_scale[2], Ey[i]*=origin_scale[2];

	    /*
	    const int NPOINTS=100000;
	    vector<double> sX, sY;
	    for (int i=0; i<=NPOINTS; i++) {
	      sX.push_back(LS.evaluate(Ex,i*1.0/NPOINTS));
	      sY.push_back(LS.evaluate(Ey,i*1.0/NPOINTS));
	    }

	   
	    vector<double> sal;
	    for (int i=0; i<sX.size(); i++)
	      sal.push_back(1.0 * i / (sX.size()-1));


	    vector<double> sEx(LS.expansion(sX,sal,n));
	    vector<double> sEy(LS.expansion(sY,sal,n));
	    */

	    /* *** test Ex =? sEx ***
	       for (int i=0; i<Ex.size(); i++)
	      cout << Ex[i] << " " << sEx[i] << "   ";
	    cout << endl;
	    */

	    /*
	    vector<double> originv1(invariant1(X,Y));
	    vector<double> eoriginv1(LS.expansion(originv1,al,n));
	    cout << "Eoriginv1: ";
	    for (int i=0; i<eoriginv1.size(); i++)
	      cout << eoriginv1[i] << " ";
	    cout << endl;
	    */

	    /*
	    vector<double> sinv1(invariant1(sX,sY));
	    vector<double> esinv1(LS.expansion(sinv1,sal,n));
	    cout << "Eresamplinv1: ";
	    for (int i=0; i<esinv1.size(); i++)
	      cout << esinv1[i] << " ";
	    cout << endl;
	    */

	    for (int i=0; i<Ex.size(); i++) {
	      exy.push_back(Ex[i]);
	      exy.push_back(Ey[i]);
	    }
	    exy[0] = -LS.evaluate(Ex,0);
	    exy[1] = -LS.evaluate(Ey,0);

	    e1quick = e1_from_exy(exy);

	    /*
	    cout << "Ebyformula: ";
	    for (int i=0; i<e1quick.size(); i++)
	      cout << setprecision(16) <<  e1quick[i] << " ";
	    cout << endl;
	    */

	    scaleTranslate(X,Y,origin_scale);

	    vector<double> rad(sqradius(X,Y));
	    er = LS.expansion(rad,al,n);

	    

	    //	    er.erase(er.begin());





	    //	    exy = scaleTranslate(Ex,Ey,origin_scale,1);
	    //	    exy = scaleTranslate(Ex,Ey,origin_scale,0);
	    

	    //	    cout << origin_scale[0] << " " << origin_scale[1] << " " << origin_scale[2] << endl;


	    
	    /*
	    for (int i=0; i<exy.size(); i+=2) {
	      double x = exy[i] * 0.2 + exy[i+1] * 0.4;
	      double y = exy[i] * (-1.0) + exy[i+1] * 3.0;
	      exy[i] = x, exy[i+1] = y;
	    }
	    */


	    
	    
	    /*	    
	    vector<double> inv1(invariant1(X,Y));
	    e1 = LS.expansion(inv1,al,n-1);
	    vector<double>inv1deg2(invariant1deg2(X,Y));
	    vector<double>aldeg2;
	    for (int i=2; i<X.size()-1; i++)
	      aldeg2.push_back(al[i]);
	    e1deg2 = LS.expansion(inv1deg2,aldeg2,n-1);

	    */

	    /*	    
	    vector<double> inv2(invariant2(X,Y));
	    e2 = LS.expansion(inv2,al,n-1);
	    */

	    //	    vector<double> inv3(invariant2(X,Y));
	    //	    e3 = LS.expansion(inv3,al,n-1);

	} else {
	  exy = vector<double>(2*n);
	  er = vector<double>(n);
	  e1 = vector<double>(n);
	  e1deg2 = vector<double>(n);
	  e1quick = vector<double>(n);
	  e2 = vector<double>(n);
	  //	  e3 = vector<double>(n);
	}
    
	out << charname << ":";
	//            out << " " << L * L;
	for (int i=0; i<n; i++) 
	  //	  out << " " << er[i];
	  out 
	     << setprecision(18) 
	    //	    	    << " " << exy[i*2] << " " << exy[i*2+1]
	    //	    << " " << er[i] << " " << e1[i] << " " << e2[i];
	    //	     << " " << e1[i] << " " << e1deg2[i] << " " << e1quick[i];
	    //	    	  << " " << e1quick[i];
	        << " " << er[i] << " " << e1quick[i];
	    //	    << " " << er[i];
	out << endl;
    }
    
    return 0;
}
