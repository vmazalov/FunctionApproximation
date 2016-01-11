#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include<map>
#include<set>
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

double sq_eucl_dist(vector<double> &p, vector<double> &q) {
  double r = 0;
  for (int i=0; i<p.size(); i++)
    r += (p[i]-q[i]) * (p[i]-q[i]);
  return r;
}

// Manhattan distance \sum |p_i - q_i|
double manh_dist(vector<double> &p, vector<double> &q) {
  double r = 0;
  for (int i=0; i<p.size(); i++)
    r += abs(p[i]-q[i]);
  return r;
}

// Chebyshev distance
double cheb_dist(vector<double> &p, vector<double> &q) {
  double r = 0;
  for (int i=0; i<p.size(); i++)
    r = max(r,abs(p[i]-q[i]));
  return r;
}

// Chebyshev-Manhattan distance
double cheb_manh_dist(vector<double> &p, vector<double> &q) {
  double r = 0;
  for (int i=0; i<p.size(); i+=2)
    r += max(abs(p[i]-q[i]),abs(p[i+1]-q[i+1]));
  return r;
}

// Chebyshev Inner distance
double cheb_inner_distance(vector<double> &p, vector<double> &q, vector<double> &x) {
  double r = 0;
  //numberic integration with trapezoid rule
  vector<double> dv;
  for (int i=0;i<p.size();i++)
    dv.push_back(p[i]-q[i]);

  for (int i=1; i<p.size()-1; i++)
    if (1-x[i]*x[i]>0 && 1-x[i-1]*x[i-1]>0)
      r += (dv[i]*dv[i]/sqrt(1-x[i]*x[i]) + dv[i-1]*dv[i-1]/sqrt(1-x[i-1]*x[i-1]))*(x[i] - x[i-1])/2;
  return r;
}

// Legendre distance
double legendre_distance(vector<double> &p, vector<double> &q, vector<double> &x) {
  double r = 0;
  //numeric integration with trapezoid rule
  vector<double> dv;
  for (int i=0;i<p.size();i++)
    dv.push_back(p[i]-q[i]);
  for (int i=1; i<dv.size(); i++){
    r += (dv[i]*dv[i] + dv[i-1]*dv[i-1])*(x[i] - x[i-1])/2;//0000000;
    //cout << p[i] << " " << q[i] << " "  << p[i-1] << " " << q[i-1] << " " << x[i]<< " " << x[i-1] << " r=" << r  << endl;
  }
  //r = r/10000000;
  return r;
}

// LS distance
double ls_distance(vector<double> &p, vector<double> &q, vector<double> &x, double mu) {
  double r = legendre_distance(p,q,x);
  //numberic integration with trapezoid rule
  vector<double> dv;
  for (int i=0;i<p.size();i++)
    dv.push_back(p[i]-q[i]);

  for (int i=2; i<p.size(); i++)
    r += mu*(x[i] - x[i-1])*( (dv[i] - dv[i-1])*(dv[i] - dv[i-1])/((x[i]-x[i-1])*(x[i]-x[i-1])) + (dv[i-1] - dv[i-2])*(dv[i-1] - dv[i-2])/((x[i-1]-x[i-2])*(x[i-1]-x[i-2])))/2;
  return r;
}



int main (int argc, char* argv[]) {
  if (argc != 3) {
    cout << "Usage: chnn-manh file-with-coordinates(input) mu(input)" << endl;
    return -1;
  }

    filebuf fb;
    fb.open(argv[1],ios::in);
    istream fls(&fb);

    double mu;
    istringstream ssmu(argv[2]);
    ssmu>> mu;
   
   // Inputs LS coeff vectors of all characters from file-with-leg-sob-coeffs
    string line;
    int dim;
    double err = 0.0;
    int cCh = 0;
    
    while (getline(fls,line)) {
      istringstream inp(line);
      string charname;
      inp >> charname;
      inp >> dim; 
      if (dim <= 2)
         continue;
      charname = charname.substr(0,charname.length()-1);
      vector<double> X(dim);
      vector<double> Y(dim);
      vector<double> Xa(dim);            
      vector<double> Ya(dim);
      for (int i=0; i<dim; i++) {
	     inp >> X[i];
	     inp >> Xa[i];
	     inp >> Y[i];
	     inp >> Ya[i];     
      }
      vector<double> al(arclen(X,Y));
      for (int i=0; i< al.size();i++)
	al[i] /= al[al.size()-1];
      err += (cheb_inner_distance(X,Xa,al) + cheb_inner_distance(Y,Ya,al))/2;
      //cout <<"err " << err << endl;
      cCh++;
    }
    
    cout << err/cCh <<  endl;
  
    return 0;
}
