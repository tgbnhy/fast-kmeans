package edu.nyu.unik.expriments;

import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;

import com.yammer.metrics.core.HealthCheck.Result;

import au.edu.rmit.trajectory.clustering.kpaths.KPathsOptimization;
import au.edu.rmit.trajectory.clustering.validity.KMeansCVIs;

public class kpathEfficiency {

	public kpathEfficiency() {
		// TODO Auto-generated constructor stub
	}
	
	/*
	 * it will read the log files and conclude the file for figure ploting
	 */
	public void ResultProcessing() {
		//Increase k to get the running time
		//increase k to get the iterations time
		//Increase data scale
	}
	
	public static int [] generateSeeds(int kMax, int traNum) {
		int []centroids = new int[kMax];
		Random xRandom = new Random();
		for(int j =0 ; j<kMax; j++) {
			centroids[j] = xRandom.nextInt(traNum);
			System.out.print(centroids[j]+",");
		}
		System.out.println();
		return centroids;
	}
	
	/*
	 * this function will test the optimizations provided for kpaths with EBD by changing k, data scale
	 */
	public static void compareOwnMethods(String[] args) throws IOException {
		KPathsOptimization run2 = new KPathsOptimization(args);
		run2.runExperiments(args);
	//	run2.runExperimentsElbow(args);
	}
	
	/*
	 * this function will test kpaths with 6 distance measures by changing k
	 */
	public static void compareDifferentDistanceMeasure(String[] args)throws IOException, SQLException, InterruptedException {
		PrintStream originalOut = System.out;
		boolean indexPivot = true;
		boolean indexInverte = true;
        String indexfilename = "./logs/validity/indices_"+args[5]+"_"+args[2];
        PrintStream indexfileOut = new PrintStream(indexfilename);
		System.setErr(indexfileOut);
        
		for(int expTimes = 0; expTimes< 10; expTimes++) {
			String groupname = "./seeds/"+expTimes+"_"+args[2];//set the same random seeds first, random multiple times
			PrintStream seedOut = new PrintStream(groupname);
			System.setOut(seedOut);
			int []centroids = generateSeeds(100, Integer.valueOf(args[2]));
			System.setOut(originalOut);
			for(int k = 10; k < 60; k+=10) {//test different k
				args[1] = Integer.toString(k);
				for(int dis = 0; dis<6; dis++) {//test different distance functions
					String filename = "./logs/efficiency/"+args[5]+"_"+dis+"_"+args[1]+"_"+args[2]+"_"+expTimes;
					PrintStream fileOut = new PrintStream(filename);
					System.setOut(fileOut);			
					KPathsOptimization run2 = new KPathsOptimization(args, dis, centroids);
					run2.staticKpath(false);
					double maxDis[] = run2.getmaxDistance();
					KMeansCVIs effectiveness = new KMeansCVIs(k, dis, run2.getCENTERS(), run2.getdatamap(), run2.getPointcoordiate(), maxDis);
					String effectiveResult = effectiveness.calculaIndices();
				//	System.err.println(effectiveResult);
					
					System.setOut(originalOut);
					filename = "./logs/validity/"+args[5]+"_"+dis+"_"+args[1]+"_"+args[2]+"_"+expTimes;
					fileOut = new PrintStream(filename);
					System.setOut(fileOut);
					run2.printResults();	//out put the results to the files for computing the validity
					System.setOut(originalOut);
				}
			}
		}
	}
	/*
	 * this class will mainly test all the clustering algorithm, and store the running time in the file and results
	 */
	public static void main(String[] args) throws IOException, SQLException, InterruptedException {
	//	compareDifferentDistanceMeasure(args);
		compareOwnMethods(args);
		
	//	run2.runIndexbuildQueue(10, 20); // the first parameter is radius, the second is the capacity
		
	//	1565595 
	//	mindex.runkpath(args);
	//	KMeansHMTree assignment = new KMeansHMTree(args);
	//	assignment.staticKpath();
	}
}
