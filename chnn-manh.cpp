#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include<map>
#include<set>
#include<lapackpp.h>
using namespace std;

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


// Inner product
double inner_prod(vector<double> &p, vector<double> &q) {
  double r = 0;
  for (int i=0; i<p.size(); i++)
    r += p[i] * q[i];
  return r;
}

// Difference of two vectors, i.e. p - q
vector<double> vect_diff(vector<double> &p, vector<double> &q) {
  vector<double> r(p.size());
  for (int i=0; i<p.size(); i++)
    r[i] = p[i] - q[i];
  return r;
}

// Product of scalar 'c' and vector 'p'
vector<double> vect_cmul(double c, const vector<double> &p) {
  vector<double> r(p.size());
  for (int i=0; i<p.size(); i++)
    r[i] = c * p[i];
  return r;
}

bool debug = false;

// The main function
// Computes distance from point 'u0' to its projection on the convex hull of points in 'p'
// Assuming that this convex hull is a simplex
// If the points happen to be not in general position, perturbs them slightly
// The perturbations are chosen randomly, until the corresponding determinant is nonzero
// In principle, this may run forever; in practice it doesn't.
// See
// http://www.orcca.on.ca/TechReports/2009/TR-09-02.html
// Tie-Breaking for Curve Multiclassifiers
// Oleg Golubitsky and Stephen M. Watt 
double project(vector<double> u0, vector<vector<double> > p) {
  if (p.size() == 0) return -1;
  if (p.size() == 1) return sqrt(sq_eucl_dist(u0,p[0]));

  int d = p.size() - 1;
  LaGenMatDouble A(d,d);
  vector<vector<double> > q(d);
  for (int i=0; i<d; i++) {
    q[i] = vect_diff(p[i+1],p[0]);
    A(i,i) = inner_prod(q[i],q[i]);
    for (int j=0; j<i; j++) 
      A(i,j) = A(j,i) = inner_prod(q[i],q[j]);
  }

  //  cout << A << endl;

  vector<double> u(vect_diff(u0,p[0]));
  LaVectorDouble B(d);
  for (int i=0; i<d; i++)
    B(i) = inner_prod(u,q[i]);

  LaVectorDouble X(d);  
  srand (1);
  while (true) {
    bool excep = false;
    try {   // Attention! This loop may run forever. If it does, try to change the constant 1e-12 below
      LaLinearSolve(A,X,B);
    } catch (LaException e) {
      if (debug)  cout << "Matrix A is singular." << endl;
      for (int i=0; i<d; i++)
	for (int j=0; j<d; j++) {
	  int s = (rand() % 2) * 2 - 1;
	  A(i,j) += s * 1e-11; // Change this constant if the loop runs forever
	}
      //      cout << i << " " << j << " " << s << endl;
      //      cout << A << endl;
      excep = true;
    }
    if (!excep) break;
  }

  //  cout << "X=" << X << endl;

  double s = 0;
  bool inside = true;
  for (int i=0; inside && i<d; i++) {
    double coef = X(i);
    if (coef < 0) inside = false;
    s += coef;    
  }
  if (s > 1) inside = false;

  //  cout << "Inside=" << inside << endl;

  if (inside) {
    vector<double> r(vect_cmul(X(0),q[0]));
    for (int i=1; i<d; i++) {
      vector<double> v(vect_cmul(-X(i),q[i]));
      r = vect_diff(r,v);
    }
    return sqrt(sq_eucl_dist(r,u));
  } else {
    /*

    double d = 1e16;
    for (int i=0; i<p.size(); i++) {
      vector<vector<double> > t;
      for (int j=0; j<p.size(); j++) 
	if (j!=i) 
	  t.push_back(p[j]);
      d = min(d, project(u0,t));
    }
    //    cout << "d=" << d << endl;
    return d;			      
    */
    
    vector<vector<double> > t;
    if (1-s >= 0) t.push_back(p[0]);
    for (int i=0; i<d; i++)
      if (X(i) >= 0) 
	t.push_back(p[i+1]);

    return project(u0,t);
  }
}

// Computes the distance from the vector of LS coeffs of character labelled 'u'
// to the convex hull of LS coeff vectors of characters whose labels are in 'p'
// 'ls' is supposed to store the LS coeff vectors of all characters involved.
double proj_dist(string &u, vector<string> &p, map<string,vector<double> > &ls) {
  vector<vector<double> > pv;
  for (int i=0; i<p.size(); i++)
    pv.push_back(ls[p[i]]);
  return project(ls[u],pv);
}

// Extracts the number of strokes from the class name.
// This number is supposed to follow the last dash ('-') in the class name
int nstrokes(string classname) {
  int p = classname.rfind("-");
  istringstream inp(classname.substr(p+1));
  int n;
  inp >> n;
  return n;
}


int main (int argc, char* argv[]) {
  if (argc != 7) {
    cout << "Usage: chnn-manh file-with-leg-sob-coeffs(input) dimension(input) file-with-training-samples(input) file-with-test-samples(input) k(#neighbors for conv hull) T(top classes selected by manh)" << endl;
    return -1;
  }

    filebuf fb;
    fb.open(argv[1],ios::in);
    istream fls(&fb);

    int dim;
    istringstream ssdim(argv[2]);
    ssdim >> dim;

    filebuf fb_trn;
    fb_trn.open(argv[3],ios::in);
    istream ftrn(&fb_trn);

    filebuf fb_tst;
    fb_tst.open(argv[4],ios::in);
    istream ftst(&fb_tst);

    int k;
    istringstream ssk(argv[5]);
    ssk >> k;

    int T;
    istringstream ssT(argv[6]);
    ssT >> T;

    map<string,vector<double> > ls; // maps each character label to its LS coeff vector

    // Inputs LS coeff vectors of all characters from file-with-leg-sob-coeffs
    string line;
    while (getline(fls,line)) {
      istringstream inp(line);
      string charname;
      inp >> charname;
      charname = charname.substr(0,charname.length()-1);
      vector<double> v(dim);
      for (int i=0; i<dim; i++) 
	inp >> v[i];
      ls[charname] = v;
    }

    map<string,string> cls;  // maps each character label to the corresponding class
    map<string,vector<string> > trn;  // maps each class label to the list of corresponding training characters 
    map<int,set<string> > classnames; // maps each number of strokes to the set of all classes with this many strokes

    // Inputs the above 3 maps from file-with-training-samples
    while(getline(ftrn,line)) {
      istringstream inp(line);
      string charname, classname;
      inp >> charname >> classname;
      trn[classname].push_back(charname);
      cls[charname] = classname;
      int ns = nstrokes(classname);
      if (classnames.count(ns) == 0) classnames[ns] = set<string>();
      classnames[ns].insert(classname);
    }

    // Inputs training samples from file-with-test-samples
    // For each training sample, ranks all classes by CHNN distance to this sample
    // Outputs the result to standard output
    while (getline(ftst,line)) {
        istringstream inp(line);
        string charname, classname;
        inp >> charname >> classname;

	//	debug = (charname == "LAVIOLA-7-116");
	// if (debug) cout << "Debugging " << charname << endl;

	int ns = nstrokes(classname); // number of strokes in the current test sample

	vector<string> v(classnames[ns].begin(),classnames[ns].end()); // list of classes with the same number of strokes
	int ncl = v.size();
	int t = ncl;

	vector<pair<double,string> > d; 
	// elements of 'd' are pairs <distance, class_name>
	// where distance is the Manhattan distance from the current test sample to a training sample from class_name
	// all training samples are considered
	for (int i=0; i<t; i++) {
	  vector<string> chars = trn[v[i]];
	  for (int j=0; j<chars.size(); j++)
	    d.push_back(make_pair(manh_dist(ls[chars[j]],ls[charname]),chars[j]));
	}

	sort(d.begin(),d.end()); // sort by distance increasingly

	map<string,vector<string> > cand;
	// maps the closest 'T' classes to the list of closest 'k' training samples from them
	int ncand = 0;
	for (int i=0; i<d.size(); i++) {
	  string chr = d[i].second;
	  string cl = cls[chr];
	  if (!cand.count(cl)) {
	    cand[cl] = vector<string>();
	    ncand++;
	  }
	  if (cand[cl].size() < k) cand[cl].push_back(chr);
	  if (ncand>=T) break;
	}



	vector<pair<string,vector<string> > > vcand(cand.begin(),cand.end());
	vector<pair<double,string> > rcand;
	for (int i=0; i<vcand.size(); i++) {
	  //	  double rho1 = conv_hull_dist(charname,vcand[i].second,ls);
	  double rho2 = 0;
	  //	  if (i<T) rho2 = proj_dist(charname,vcand[i].second,ls);
	  rho2 = proj_dist(charname,vcand[i].second,ls);
	  // rho2 = manh_dist(ls[charname],ls[vcand[i].second[0]]);
	  //	  if (rho1 != rho2) cout << rho1 << " " << rho2 << endl;
	  //	  else cout << "Equal rho's!!!" << endl;
	  rcand.push_back(make_pair(rho2,vcand[i].first));
	}
	sort(rcand.begin(),rcand.end());

	set<string> winners;
	cout << charname;
	for (int r=0; r<rcand.size() && r<T; r++) {
	  string winner = rcand[r].second;
	  cout << "   " << winner << " " << rcand[r].first;
	  winners.insert(winner);
	}

	for (int r=0; r<d.size(); r++) {
	  string winner = cls[d[r].second];
	  if (!winners.count(winner)) {
	    cout << "   " << winner << " " << -2;
	    winners.insert(winner);
	  }
	}

	for (int r=0; r<ncl; r++)
	  if (!winners.count(v[r]))
	    cout << "   " << v[r] << " -1";
	cout << endl;
    }
    return 0;
}
