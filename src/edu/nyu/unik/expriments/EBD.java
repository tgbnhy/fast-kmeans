package edu.nyu.unik.expriments;

import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import au.edu.rmit.trajectory.clustering.kpaths.CoOccurrence;
import au.edu.rmit.trajectory.clustering.kpaths.KPathsOptimization;
import au.edu.rmit.trajectory.clustering.kpaths.Util;
import scala.reflect.api.Trees.NewExtractor;

public class EBD {

	public EBD() {
		// TODO Auto-generated constructor stub
	}
	
	/*
	 * conduct the top-k search, check with EDR by changing k and EDR
	 * compute the distance one by one
	 * output the running time
	 */
	public static void testEffectiveness(Map<Integer, int[]> datamap, int n, Map<Integer, double[]> idPoints, int k) {
		double [][]matrixEDR = testEfficiency(datamap, n, 1, idPoints);
		double [][]matrixEBD = testEfficiencySorted(datamap, n, 0, idPoints);
		//compute the percentage based on the matrix
		int count = 0, countall = 0, countNotEqual=0;
		double error =0;
		Multiset<Double> distribution = HashMultiset.create();
		for ( int i =0 ; i < n; i++) {
			for(int j =0 ; j< n; j++) {
				if(matrixEDR[i][j]==matrixEBD[i][j])
					count++;
				else {
					double error_dis = (double)(matrixEDR[i][j]-matrixEBD[i][j])/matrixEDR[i][j];
					error += error_dis;
					distribution.add(Math.round(error_dis * 1000.0)/1000.0);	
					countNotEqual++;
				}
				countall++;
			}
		}
		double percentage = 0;
		double percentage1 = 0;
		int centroid[] = new int[k];
		Random xrandom = new Random();
		for(int i=0; i<k; i++) {
			centroid[i] = xrandom.nextInt(n);
		}
		double countAss = 0;
		double sumEDRdis=0;
		double sumEBDdis=0;
		
		for(int i = 0; i<n; i++) {
			double minEBD = Double.MAX_VALUE;
			int jEBD = 0;
			double minEDR = Double.MAX_VALUE;
			int jEDR = 0;
			int[] trajectoryB = datamap.get(i);
			for(int j = 0; j<k; j++) {
				int[] trajectoryA = datamap.get(centroid[j]);
				//wrong, both trajectories should be sorted, we can get it from the matrix;
			//	double EBD = Util.Intersection(trajectoryA, trajectoryB, trajectoryA.length, trajectoryB.length);
				double EBD = matrixEBD[i][centroid[j]];
				if(minEBD>EBD) {
					minEBD = EBD;
					jEBD = j;
				}
			//	double EDR = Util.EDR(trajectoryA, trajectoryB);
				double EDR = matrixEDR[i][centroid[j]];
				if(minEDR>EDR) {
					minEDR = EDR;
					jEDR = j;
				}
			//	System.out.println(EBD+"\t"+EDR);
			}
			sumEDRdis += minEDR;
			sumEBDdis += minEBD;
			if(jEBD != jEDR)
				countAss++;
		//	assignEBD.put(i, jEBD);
		//	assignEDR.put(i, jEDR);
			
			
		//	System.out.println(datamap.get(i).length);
			double a[] = matrixEDR[i];
		//	System.out.println(Arrays.toString(a));
			Map<Integer, Double> aMap = new HashMap<>();
			int counter = 0;
			double sum = 0;
			for(double dis: a) {
				sum += dis;
				aMap.put(counter++, dis);
			}
			
		//	System.out.println(aMap.values());
			 Map<Integer, Double> result2 = new LinkedHashMap<>();
			 aMap.entrySet().stream()
		                .sorted(Map.Entry.<Integer, Double>comparingByValue())
		                .forEachOrdered(x -> result2.put(x.getKey(), x.getValue()));
		//	System.out.println(result2.values());
			
			ArrayList<Integer> ListEBD = new ArrayList<>();
			int cou = 0;
			for(int j: result2.keySet()) {
				ListEBD.add(j);
				if(++cou>=k){
					break;
				}
			}
			
			double b[] = matrixEBD[i];
			TreeMap<Integer, Double> bMap = new TreeMap<>();
			counter = 0;
			ArrayList<Integer> ListEDR = new ArrayList<>();
			double sumebd = 0;
			for(double dis: b) {
				sumebd += dis;
				bMap.put(counter++, dis);
			}
			Map<Integer, Double> result3 = new LinkedHashMap<>();
			bMap.entrySet().stream()
		                .sorted(Map.Entry.<Integer, Double>comparingByValue())
		                .forEachOrdered(x -> result3.put(x.getKey(), x.getValue()));
			cou = 0;
			for(int j: result3.keySet()) {
				ListEDR.add(j);
				if(++cou>=k)
					break;
			}
			ListEBD.retainAll(ListEDR);//intersection
			percentage += ListEBD.size()/(double)k;
			percentage1 += sumebd/sum;
		}
		percentage /= n;// the similarity search precision.
		percentage1 /= n;//the sum distance of matrix
		double percentage2 = 1 - countAss/n;// the assignment
		double percentage3 = sumEBDdis/sumEDRdis;// the assignment distance
		double percentage4 = (double)count/countall;
		double percentage5 = error/countNotEqual;
		System.out.print("precion: ");
		System.out.printf("%.3f", percentage);
		System.out.print(" ");
		System.out.printf("%.3f", percentage1);
		System.out.print(" ");
		System.out.printf("%.3f", percentage2);
		System.out.print(" ");
		System.out.printf("%.3f", percentage3);
		System.out.print(" ");
		System.out.printf("%.3f", percentage4);
		System.out.print(" ");
		System.out.printf("%.3f", percentage5);
		System.out.println();
		for(double a: distribution.elementSet()) {
			System.out.printf(a+"\t"+distribution.count(a)+"\n");
		}
	}
	
	/*
	 * build the distance matrix
	 * compute the distance one by one
	 * output the running time
	 */
	public static double [][] testEfficiency(Map<Integer, int[]> datamap, int n, int simi, Map<Integer, double[]> idPoints) {
		double [][]matrix = new double[n][n];
		Long time1 = System.nanoTime();
		for ( int i =0 ; i < n; i++) {
			int[] trajectoryA = datamap.get(i);
			for(int j =0 ; j< n; j++) {
				if(i>j) {
					matrix[i][j] = matrix[j][i];
					continue;
				}
				if(i==j) {
					matrix[i][j] = 0;
					continue;
				}					
				int[] trajectoryB = datamap.get(j);
				if(simi == 0) {
					matrix[i][j] = Util.Intersection(trajectoryA, trajectoryB, trajectoryA.length, trajectoryB.length);
				}else if(simi == 1) {
				//	System.out.print(i+" "+j);
					matrix[i][j] = Util.EDR(trajectoryA, trajectoryB);
				}else if(simi == 2) {
					matrix[i][j] = Util.DTW(trajectoryA, trajectoryB, idPoints);
				}else if(simi == 3) {
					matrix[i][j] = Util.Frechet(trajectoryA, trajectoryB, idPoints);
				}else if(simi == 4) {
					matrix[i][j] = Util.Hausdorff(trajectoryA, trajectoryB, idPoints);
				}else if(simi == 5) {
					matrix[i][j] = Util.ERP(trajectoryA, trajectoryB, idPoints);
				}
			}
		}
		Long time2 = System.nanoTime();
		System.out.print("matrix consumes time: ");
		System.out.printf("%.3f", (time2-time1)/1000000000.0);
		System.out.println();
		return matrix;
	}
	
	/*
	 * build the distance matrix
	 * compute the distance one by one
	 * output the running time
	 */
	public static double [][] testEfficiencySorted(Map<Integer, int[]> datamap, int n, int simi, Map<Integer, double[]> idPoints) {
		double [][]matrix = new double[n][n];
		Long time1 = System.nanoTime();
		for ( int i =0 ; i < n; i++) {
			int[] trajectoryA = datamap.get(i);
			Arrays.sort(trajectoryA);
			for(int j =0 ; j< n; j++) {
				if(i>j) {
					matrix[i][j] = matrix[j][i];
					continue;
				}
				if(i==j) {
					matrix[i][j] = 0;
					continue;
				}
				int[] trajectoryB = datamap.get(j);
				Arrays.sort(trajectoryB);
				matrix[i][j] = Util.Intersection(trajectoryA, trajectoryB, trajectoryA.length, trajectoryB.length);
			}
		}
		Long time2 = System.nanoTime();
		System.out.print("matrix consumes time: ");
		System.out.printf("%.3f", (time2-time1)/1000000000.0);
		System.out.println();
		return matrix;
	}
	
	/*
	 * test the performance related to EBD
	 */
	public static void main(String[] args) throws IOException, SQLException, InterruptedException {
		/*
		 * test efficiency
		 */
		KPathsOptimization run2 = new KPathsOptimization(args, 1);//0: sorted, >0: unsorted
		run2.loadDataOnly(false);
		Map<Integer, int[]> datamap = run2.getdatamap();
		CoOccurrence occurrence = new CoOccurrence();
		/*
		 * Before running occurrence, we need to set KPathsOptimization(args, 1);
		 */
		System.out.println(datamap.size());
		occurrence.buildMatrix(datamap, datamap.size()+"_cooccurance.txt");//

		
		Map<Integer, double[]> idPoints = run2.getPointcoordiate();
		int tranum = run2.getTraNum();
		String filename = "./logs/ebd/"+args[5]+"_"+0+"_"+tranum+".log";
		PrintStream fileOut = new PrintStream(filename);
		System.setOut(fileOut);
//		testEfficiency(datamap, tranum, 0, idPoints);
		
		run2 = new KPathsOptimization(args, 1);// unsorted for dynamic programing
		run2.loadDataOnly(false);
		datamap = run2.getdatamap();//unsorted map
		idPoints = run2.getPointcoordiate();
		run2.getTraNum();
/*		for(int i= 1 ; i< 6; i++) {
			filename = "./logs/ebd/"+args[5]+"_"+i+"_"+tranum+".log";
			fileOut = new PrintStream(filename);
			System.setOut(fileOut);
			testEfficiency(datamap, tranum, i, idPoints);
		}*/
				
		for(int k= 10; k< 60; k+= 10) {
			filename = "./logs/ebd/"+args[5]+"_"+k+"_"+tranum+".effective";
			fileOut = new PrintStream(filename);
			System.setOut(fileOut);
			testEffectiveness(datamap, tranum, idPoints, k);
		}
	}
}
