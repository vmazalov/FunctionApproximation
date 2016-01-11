
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealMatrixImpl;

import API.InkAPI;

public class ChnnManh {
    
    private static double project(double[] u0, ArrayList<double[] > p) {
	if (p.size() == 0) return -1;
	if (p.size() == 1) return Math.sqrt(InkAPI.SqEucleadDistance(u0,p.get(0)));
	
	int d = p.size() - 1;
	double[][] mat = new double[d][d];
	
	ArrayList<double[]> q = new ArrayList<double[]>(); 
	for (int i=0; i<d; i++) {
	    q.add(InkAPI.vect_diff(p.get(i+1),p.get(0)));
	    mat[i][i] = InkAPI.inner_prod(q.get(i),q.get(i));
	    for (int j=0; j<i; j++) 
		mat[i][j] =mat[j][i] = InkAPI.inner_prod(q.get(i),q.get(j));
	}
	RealMatrix matrix = new RealMatrixImpl(mat);
	
	double[] u = InkAPI.vect_diff(u0,p.get(0));
	
	double[] B = new double[d];
	for (int i=0;i<d;i++)
	    B[i] = InkAPI.inner_prod(u, q.get(i));
	
	double[] X = new double[d];
	Random randomGenerator = new Random();
	while (true) {
	    boolean excep = false;
	    try {   // Attention! This loop may run forever. If it does, try to change the constant 1e-12 below
		X=  matrix.solve( B );
	    } catch (InvalidMatrixException e) {
		//       System.out.println("Matrix A is singular!");
		for (int i=0; i<d; i++)
		    for (int j=0; j<d; j++) {
			int s = (randomGenerator.nextInt() % 2) * 2 - 1;  
			mat[i][j] += s * 1e-11; // Change this constant if the loop runs forever
		    }
		matrix.setSubMatrix(mat, 0, 0);
		excep = true;
	    }
	    if (!excep) break;
	}
	
	double s = 0;
	boolean inside = true;
	for (int i=0; inside && i<d; i++) {
	    double coef = X[i];
	    if (coef < 0) inside = false;
	    s += coef;    
	}
	if (s > 1) inside = false;
	
	if (inside) {
	    double[] r = InkAPI.vect_cmul(X[0],q.get(0));
	    for (int i=1; i<d; i++) {
		double[] v = InkAPI.vect_cmul(-X[i],q.get(i));
		r = InkAPI.vect_diff(r,v);
	    }
	    return Math.sqrt(InkAPI.SqEucleadDistance(r,u));
	} else {
	    ArrayList<double[]> t = new ArrayList<double[] > ();
	    if (1-s >= 0) t.add(p.get(0));
	    for (int i=0; i<d; i++)
		if (X[i] >= 0) 
		    t.add(p.get(i+1));
	    
	    return project(u0,t);
	}
    }
    
    // Computes the distance from the vector of LS coeffs of character labelled 'u'
    // to the convex hull of LS coeff vectors of characters whose labels are in 'p'
    // 'ls' is supposed to store the LS coeff vectors of all characters involved.
    static double proj_dist(String u, ArrayList<String> p, Map<String,double[] > ls) {
	ArrayList<double[] > pv = new ArrayList<double[] >();
	for (String el: p)
	    pv.add(ls.get(el));
	return project(ls.get(u),pv);
    }
    
    // Extracts the number of strokes from the class name.
    // This number is supposed to follow the last dash ('-') in the class name
    static int nstrokes(String classname) {
	return Integer.parseInt(classname.substring(classname.lastIndexOf("-")+1));
    }
    
    public static void main(String[] args) throws FileNotFoundException {
	if (args.length != 6) {
	    System.out.println("Usage: chnn-manh file-with-leg-sob-coeffs(input) dimension(input) file-with-training-samples(input) file-with-test-samples(input) k(#neighbors for conv hull) T(top classes selected by manh)");
	    System.exit(0);
	}
	
	Map<String,double[] > ls = new HashMap<String,double[] >(); // maps each character label to its LS coeff vector
	Scanner fls = new Scanner(new File(args[0]));
	int dim = Integer.parseInt(args[1]);
	Scanner ftrn = new Scanner(new File(args[2]));
	Scanner ftst = new Scanner(new File(args[3]));
	int k = Integer.parseInt(args[4]);
	int T = Integer.parseInt(args[5]);
	while (fls.hasNextLine()) {
	    Scanner line = new Scanner(fls.nextLine());
	    String charname = line.next();
	    charname = charname.substring(0,charname.length()-1);
	    double[] v = new double[dim];
	    for (int i=0; i<dim; i++) 
		v[i] = line.nextDouble();
	    ls.put(charname, v);
	}
	    
	Map<String,String> cls = new HashMap<String,String>();  // maps each character label to the corresponding class
	Map<String,ArrayList<String> > trn = new HashMap<String,ArrayList<String> >();  // maps each class label to the list of corresponding training characters 
	Map<Integer,Set<String> > classnames = new HashMap<Integer,Set<String> >(); // maps each number of strokes to the set of all classes with this many strokes

	// Inputs the above 3 maps from file-with-training-samples
	while(ftrn.hasNextLine()) {
	    Scanner line = new Scanner(ftrn.nextLine());
	    String charname = line.next();
	    String classname = line.next();
	    if (trn.get(classname) == null)
		trn.put(classname, new ArrayList<String>());
	    trn.get(classname).add(charname);
	    cls.put(charname, classname);
	    int ns = nstrokes(classname);
	    if (classnames.get(ns) == null) 
		classnames.put(ns, new HashSet<String>());
	    classnames.get(ns).add(classname);
	}
	    
	// Inputs training samples from file-with-test-samples
	// For each training sample, ranks all classes by CHNN distance to this sample
	// Outputs the result to standard output
	while (ftst.hasNext()) {
	    Scanner line = new Scanner(ftst.nextLine());
	    String charname = line.next();
	    String classname = line.next();
	    int ns = nstrokes(classname); // number of strokes in the current test sample

	    Set<String> v = classnames.get(ns); // list of classes with the same number of strokes
	    ArrayList<Pair<Double, String> > d = new ArrayList<Pair<Double, String> >(); //
	    // elements of 'd' are pairs <distance, class_name>
	    // where distance is the Manhattan distance from the current test sample to a training class from class_name
	    // all training samples are considered
	    for (String el:v) {
		ArrayList<String> chars = trn.get(el);
		for (String ch:chars)
		    d.add(new Pair<Double, String>(InkAPI.manh_dist(ls.get(ch),ls.get(charname)),ch));
	    }
	    Collections.sort(d, new PairsComparator());
	    Map<String,ArrayList<String> > cand = new HashMap<String,ArrayList<String> >();
	    // maps the closest 'T' classes to the list of closest 'k' training samples from them
	    int ncand = 0;
	    for (Pair<Double, String> pair: d) {
		String chr = (String)pair.getSecond();
		String cl = cls.get(chr);
		if (cand.get(cl) == null) {
		    cand.put(cl, new ArrayList<String>());
		    ncand++;
		}
		if (cand.get(cl).size() < k) 
		    cand.get(cl).add(chr);
		if (ncand>=T) 
		    break;
	    }

	    ArrayList<Pair<Double,String> > rcand = new ArrayList<Pair<Double,String> >();
	    for(Map.Entry<String,ArrayList<String>> entry : cand.entrySet()){
		double rho2 = 0;
		rho2 = proj_dist(charname,entry.getValue(),ls);
		rcand.add(new Pair(rho2,entry.getKey()));
	    }
	    Collections.sort(rcand, new PairsComparator());
	    Set<String> winners = new HashSet<String>();
	    System.out.print(charname);
	    for (int r=0; r<rcand.size() && r<T; r++) {
		String winner = rcand.get(r).getSecond();
		System.out.print("   " + winner + " " + (double)Math.round(1000000*rcand.get(r).getFirst())/1000000);
		winners.add(winner);
	    }
	    System.out.println();
	}
    }

}
