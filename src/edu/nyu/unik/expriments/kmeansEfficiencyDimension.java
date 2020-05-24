package edu.nyu.unik.expriments;

import java.io.IOException;

import au.edu.rmit.trajectory.clustering.kmeans.kmeansAlgorithm;
import edu.wlu.cs.levy.cg.KeyDuplicateException;
import edu.wlu.cs.levy.cg.KeySizeException;

public class kmeansEfficiencyDimension {

	public kmeansEfficiencyDimension() {
		// TODO Auto-generated constructor stub
	}

	public static void main(String[] args) throws IOException, KeySizeException, KeyDuplicateException {
		// TODO Auto-generated method stub
		int[] kvalue = new int[]{10, 50, 100, 200, 400, 600, 800, 1000};//10, 100, 1000
		int[] scales = new int[] {1000, 5000, 10000, 50000, 100000, 500000, 1000000};
		int[] capacities = new int[] {10, 20, 30, 40, 50, 60};//capacity
		int[] dimensions = new int[] {10, 20, 30, 40, 50};//dimension
		kvalue = new int[]{10};//10, 100, 1000
		scales = new int[] {1000};
		capacities = new int[] {30};
	//	dimensions = new int[] {2};// test the parameters
		int testTime = 1;//test one time
		kmeansAlgorithm<?> runkmeans = new kmeansAlgorithm<>(args);
		for(int dimension: dimensions)
		for(int capacity: capacities)
		for(int scale: scales) {
			//set trajectory number 
			runkmeans.setDimension(dimension);
		//	runkmeans.setCapacity(capacity);
		//	runkmeans.setScale(scale);
			runkmeans.experiments(kvalue, testTime);
		}
	//	runkmeans.staticKmeans(false, true, false);//index sign, bound, scan whole tree again
	//	runkmeans.staticKmeans(true, false, true);//test the functionality		
		//run script to draw the figures using gnuplot
		
	}
}