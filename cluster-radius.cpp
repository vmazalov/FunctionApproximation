/* Compute cluster radius: average and standard deviation 
Reads from files with LS coefs, Train and Test samples*/
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<map>
#include<set>
using namespace std;

typedef map<string,vector<string> > StringVectorStringMap;
typedef map<string,vector<double> > StringVectorDoubleMap;

bool findChar(StringVectorStringMap &cls, string charname){
  bool res = false;
  StringVectorStringMap::iterator pos;
      for (pos = cls.begin(); pos != cls.end(); ++pos) {
          vector<string> chars(pos->second);
          for (int i=0;i<chars.size(); i++) {
        	  if (chars[i] == charname){
        		  return true;
        	  }
          }
      }

  return res;
}

StringVectorDoubleMap ComputeCentersOfClasses(StringVectorStringMap &cls, StringVectorDoubleMap &chs, int D) {
	StringVectorDoubleMap res;
	StringVectorStringMap::iterator pos;
	for (pos = cls.begin(); pos != cls.end(); ++pos) {
		  vector<double> clsAv(D);
		  int c=0;
		  vector<string> chars(pos->second);
		  for (int i=0;i< D ; i++) {
			  vector<double> xy(chs[chars[i]]);
			  for (int j=0;j<xy.size();j++){
				  clsAv[j] += xy[j];
			  }
			  c++;
		  }
		  for (int i=0;i<clsAv.size();i++)
			  clsAv[i] /= c;

		  res[pos->first] = clsAv;
	}
	return res;
}

int main (int argc, char* argv[]) {
  if (argc != 5) {
    cout << "Usage: legsobjoin leg-sob-coefs(input) trn-cv0(input) cv0(input) Dimension(input) " << endl;
    return -1;
  }

    filebuf fb;
    filebuf fb2;
    filebuf fb3;

    fb.open(argv[1],ios::in);
    istream fin(&fb);
    
    fb2.open(argv[2],ios::in);
    istream trn0(&fb2);

    fb3.open(argv[3],ios::in);
    istream cv0(&fb3);

    int D;
	istringstream ssT(argv[4]);
	ssT >> D;

    string line;
    StringVectorStringMap cls;
    while (getline(trn0,line)) {
        istringstream inp(line);
        string charname;
        string classname;
        inp >> charname;
        inp >> classname;
        cls[classname].push_back(charname);
    }

    while (getline(cv0,line)) {
        istringstream inp(line);
        string charname;
        string classname;
        inp >> charname;
        inp >> classname;
        cls[classname].push_back(charname);
    }
    
    //map of charname to the vector of coefs
    StringVectorDoubleMap chs;
    while (getline(fin,line)) {
        istringstream inp(line);
        string charname;
        inp >> charname;
        charname = charname.substr(0,charname.length()-1);
        if (!findChar(cls,charname))  continue;
        vector<double> xy;
        while(!inp.eof()) {
          double x,y;
          inp >> x >> y;
		  xy.push_back(x);
		  xy.push_back(y);
	    }
	    chs[charname]=xy;
    }

    StringVectorDoubleMap centers = ComputeCentersOfClasses(cls, chs, D);
    StringVectorDoubleMap::iterator pos;
    	for (pos = centers.begin(); pos != centers.end(); ++pos) {
    		cout << pos->first;
    		vector<double> cen(pos->second);
    		for (int i=0;i<cen.size();i++)
    			cout << " " << cen[i];
    		cout << endl;
    	}
    return 0;
}
