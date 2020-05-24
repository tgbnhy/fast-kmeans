package au.edu.rmit.trajectory.clustering.kpaths;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

import org.checkerframework.framework.qual.FromByteCode;

import au.edu.rmit.trajectory.clustering.TrajectoryMtree;
import au.edu.rmit.trajectory.clustering.validity.KMeansCVIs;
import au.edu.rmit.trajectory.clustering.visualization.mapv;
import java_cup.internal_error;
import org.jgrapht.*;
import org.jgrapht.alg.interfaces.ShortestPathAlgorithm.SingleSourcePaths;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
/*
 * this class will mainly optimize the class process on assignment and refinement
 */
public class KPathsOptimization<T> extends KPaths {
	private Thread runingThread;//multi thread
	private String threadName;
	TrajectoryMtree centroidindex;//the tree used to group the centroid in the first
	protected Map<Integer, Double> center_drift = null;
	protected double[] group_drift = null;
	protected int []centroids = null;
	protected int groupIteration = 1;// index building iteration to form small group
	
	protected Map<Integer, ArrayList<Integer>> groupsBuilt;
	
	String groupFilename;// the file used to store the leaf file	
	
	public KPathsOptimization(String datapath) {
		super(datapath);
		threadName = datapath;
		iterationStops = true;
	}
	
	public KPathsOptimization(String []datapath) {
		super(datapath);
		iterationStops = true;
		invertedIndex = false;
		distanceModel = 0;
	}
	
	public KPathsOptimization(String []datapath, int dismodel) {
		super(datapath);
		iterationStops = true;
		invertedIndex = false;
		distanceModel = dismodel;
		if(distanceModel>0)
			invertedIndex = false;
	}
	
	public Map<Integer, ArrayList<Integer>> getGroupsBuilt(){
		return groupsBuilt;
	}
	/*
	 * given centroids pool in the beginning, this is for experiment only
	 */
	public KPathsOptimization(String []datapath, int dismodel, int []centroids) {
		super(datapath);
		iterationStops = true;
		invertedIndex = false;
		distanceModel = dismodel;
		if(distanceModel>0)
			invertedIndex = false;
		this.centroids = centroids;
	}
	
	public void initializeClustersRandom(int k){
		if(centroids==null) {//if there is no given centroids
			Random rand = new Random();
			for(int t=0; t<k; t++) {
				int  n = rand.nextInt(trajectoryNumber);
				int[] cluster = datamap.get(n);
				ClusterPath cl = new ClusterPath(cluster, n);
				CENTERS.add(cl);
			}
		}else {
			for(int t=0; t<k; t++) {
				int  n = centroids[t];
				int[] cluster = datamap.get(n);
				ClusterPath cl = new ClusterPath(cluster, n);
				CENTERS.add(cl);		
			}
		}
	}
	/*
	 * 
	 */
	public void initializeClusterskMeansPlusPlus(int k){
		int fartraid = -1;
		if(centroids==null) {//if there is no given centroids
			Random rand = new Random();
			for(int t=0; t<k; t++) {
				if(t==0) {
					int n = rand.nextInt(trajectoryNumber);
					int[] cluster = datamap.get(n);
					ClusterPath cl = new ClusterPath(cluster, n);
					CENTERS.add(cl);
				}else {//select the farest.
					int[] cluster = datamap.get(fartraid);
					ClusterPath cl = new ClusterPath(cluster, fartraid);
					CENTERS.add(cl);
				}
				double maxdis = 0;
				for(int id: datamap.keySet()) {				// update the nearest centroids for every point
					double mindis = Double.MAX_VALUE;
					for(ClusterPath ce: CENTERS) {
						double dis = Util.Intersection(datamap.get(id), ce.getTrajectoryData());
						if(dis<mindis) {
							mindis = dis;
						}
					}
					if(mindis>maxdis) {
						maxdis = mindis;
						fartraid = id;
					}
				}
			//	System.out.println(fartraid);
			}
		}else {
			for(int t=0; t<k; t++) {
				int  n = centroids[t];
				int[] cluster = datamap.get(n);
				ClusterPath cl = new ClusterPath(cluster, n);
				CENTERS.add(cl);		
			}
		}
	}
	
	public void initializeClustersRandomStreaming(int k){
		ArrayList<Integer> arrayList = new ArrayList<>(datamap.keySet());
		
		if(CENTERS.isEmpty() ) {//if there is no given centroids
			Random rand = new Random();
			for(int t=0; t<k; t++) {
				int  n = rand.nextInt(trajectoryNumber);
				int traid = arrayList.get(n);
				int[] cluster = datamap.get(traid);
				ClusterPath cl = new ClusterPath(cluster, traid);
				CENTERS.add(cl);
			}
		}else {
			ArrayList<Integer> candidates = new ArrayList<>();
			for(ClusterPath aClusterPath: CENTERS) {
				candidates.add(aClusterPath.getTrajectoryID());
			}
			CENTERS = new ArrayList<ClusterPath>();	
			for(int t=0; t<k; t++) {
				int traid = candidates.get(t);
				int[] cluster = datamap.get(traid);
				ClusterPath cl = new ClusterPath(cluster, traid);
				CENTERS.add(cl);
			}
		}
	}
	
	public void initializeClustersRandom1(int k, boolean graph, ArrayList<Integer> nonImportantRoads){
		if(centroids==null) {//if there is no given centroids
			Random rand = new Random();
			for(int t=0; t<k; t++) {
				int  n = rand.nextInt(trajectoryNumber);
				int[] cluster = datamap.get(n);
				ClusterPath cl = new ClusterPath(cluster, n, graph, nonImportantRoads);
				CENTERS.add(cl);		
			}
		}else {
			for(int t=0; t<k; t++) {
				int  n = centroids[t];
				int[] cluster = datamap.get(n);
				ClusterPath cl = new ClusterPath(cluster, n, graph, nonImportantRoads);
				CENTERS.add(cl);		
			}
		}
	}
	
	public void start() {
		System.out.println("Starting " + threadName);
		if (runingThread == null) {
			runingThread = new Thread(this, threadName);
			runingThread.start();
		}
	}
	
	
	/*
	 * divide k clusters into t groups when k not equals to t
	 */
	protected void groupInitialClusters(int t, int k) throws IOException {
		group = new HashMap<>();
		centerGroup = new HashMap<>();
		if(t==k) {// when k is small, we do not divide into too many groups
			for(int i = 0;i<k; i++) {
				ArrayList<Integer> a = new ArrayList<>();
				a.add(i);
				group.put(i, a);
				centerGroup.put(i, i);
			}
		}else {
			String LOG_DIR = "./seeds/Groupcentroid.txt";
			PrintStream fileOut = new PrintStream(LOG_DIR);
			System.setOut(fileOut);	
			for(int i = 0; i<k; i++) {
				System.out.print((i+1)+"\t");
				int[] aa = CENTERS.get(i).getTrajectoryData();
				int counter = 0;
				for(int id: aa) {
					System.out.print(id);
					if(counter++ < aa.length-1)
						System.out.print(",");
				}
				System.out.println();
			}
			String LOG_DIR1 = "./seeds/Groupcentroid1.txt";
			PrintStream fileOut1 = new PrintStream(LOG_DIR1);
			System.setOut(fileOut1);
			String args[] = new String[5];
			args[0] = LOG_DIR;
			args[1] = Integer.toString(t);
			args[2] = Integer.toString(k);
			args[3] = edgefile;
			args[4] = graphfile;
			KPathsOptimization<?> run2 = new KPathsOptimization(args);
			run2.staticKpath(false);//why there is bug
			group = run2.getGroupsBuilt();
			for(int i = 0;i<t; i++) {
				ArrayList<Integer> centers = group.get(i);
				for(int center: centers) {
					centerGroup.put(center, i);
				}
			}
		}
	}
	
	/*
	 * update the lower bound of trajectory toward group i
	 */
	protected void updateSingleLowerBound(Map<Integer, double[]> trajectoryBounds, int traid, int group_i, double newbound) {
		double [] bounds = trajectoryBounds.get(traid);
		bounds[group_i+2] = newbound;
		trajectoryBounds.put(traid, bounds);
	}
	
	/*
	 * compute the drift between current center and previous center
	 */
	protected void computeDrift(int k, int t) {
		group_drift = new double[t];
		if(center_drift.isEmpty())//if it is the first time to compute
			for(int i=0; i<k; i++) {
				int [] alist =  CENTERS.get(i).getTrajectoryData();;
				int [] blist = PRE_CENS.get(i).getTrajectoryData();
				double dis = computeRealDistance(alist, blist, 0);
				center_drift.put(i, dis);
			}
		for(int group_i = 0; group_i < t; group_i++) {		//choose the minimum one as the group
			ArrayList<Integer> centers = group.get(group_i);
			double max_drift = 0;
			for(int centerid: centers) {
				if(max_drift<center_drift.get(centerid)) {
					max_drift = center_drift.get(centerid);
				}
			}
			group_drift[group_i] = max_drift;
		}
	}
	
	/*
	 * compute the real distance use different model, DTW, Frechet, EDR, and EBD
	 */
	public double computeRealDistance(int []tra, int[] clustra, int idx) {
		if(tra==null) {
			long startTime = System.nanoTime();							
			tra = datamap.get(idx);//  read the trajectory data
			long endtime = System.nanoTime();
			runrecord.addIOTime((endtime-startTime)/1000000000.0);
		}
		numCompute++;//it shows the number of distance computations
		long Time1 = System.nanoTime();
		double min_dist = 0;
		if (distanceModel == 0) {//means the metric
			min_dist = (double)Intersection(tra, clustra, tra.length, clustra.length);			
			if(maxDistance[0]<min_dist) {
				maxDistance[0] = min_dist;//this is for the normalization
			}
		}
		else if( distanceModel == 1) {
			min_dist = Util.EDR(tra, clustra);
			if(maxDistance[1]<min_dist)
				maxDistance[1] = min_dist;
		}else if (distanceModel == 2){
			min_dist = Util.DTW(tra, clustra, idPoints);
			if(maxDistance[2]<min_dist)
				maxDistance[2] = min_dist;
		}else if (distanceModel == 3){
			min_dist = Util.Frechet(tra, clustra, idPoints);
			if(maxDistance[3]<min_dist)
				maxDistance[3] = min_dist;
		}else if (distanceModel == 4){
			min_dist = Util.Hausdorff(tra, clustra, idPoints);
			if(maxDistance[4]<min_dist)
				maxDistance[4] = min_dist;
		}else if (distanceModel == 5){
			min_dist = Util.ERP(tra, clustra, idPoints);
			if(maxDistance[5]<min_dist)
				maxDistance[5] = min_dist;
		}
		long Time2 = System.nanoTime();
		runrecord.addsimiComputationTime((Time2 - Time1) / 1000000000.0);
		return min_dist;
	}

	/*
	 * update the histogram when tra moves from oldcenter to the newcenter
	 */
	public void accumulateHistogramGuava(int []tra, int idx, int newCenter, int oldCenter) {
		if(tra==null) {//never scan
			tra = datamap.get(idx);// the trajectory data is read in the first time in this iteration
		}
		ClusterPath newCluster = CENTERS.get(newCenter);//update the entry in min_id						
		newCluster.updateHistorgramGuava(tra, idx);// add the new trajectory in
		ClusterPath oldCluster = CENTERS.get(oldCenter);
		oldCluster.removeHistorgramGuava(tra, idx);//remove the old trajectory out
	}
	
	/*
	 * update the trajectories in each clusters
	 */
	public void updateCenters(Map<Integer, ArrayList<Integer>> idxNeedsIn, Map<Integer, ArrayList<Integer>> idxNeedsOut) {
		long Time1 = System.nanoTime();
		for(int idx: idxNeedsIn.keySet()) {
			ArrayList<Integer> idxs = idxNeedsIn.get(idx);
			ClusterPath newCluster = CENTERS.get(idx);
			newCluster.mergeTrajectoryToCluster(idxs);
		}
		for(int idx: idxNeedsOut.keySet()) {
			ArrayList<Integer> idxs = idxNeedsOut.get(idx);
			ClusterPath newCluster = CENTERS.get(idx);
			newCluster.removeTrajectoryToCluster(idxs);
		}
		long Time2 = System.nanoTime();
		runrecord.addHistorgramTime((Time2-Time1)/1000000000.0);			
	}
	
	public double getMinimumLowerbound(double [] bounds, int groupNumber) {
		double lowerboud = Double.MAX_VALUE;
		for(int group_j=0; group_j<groupNumber; group_j++) {//get the minimum lower bound of all group
			double lowerboud_temp = Math.abs(bounds[group_j+2] - group_drift[group_j]);
			if(lowerboud_temp<lowerboud)
				lowerboud = lowerboud_temp;
		}
		return lowerboud;
	}
	
	public boolean checkInvertedIndex(Set<Integer> candilist, int idx) {
		long startTime1 = System.nanoTime();	
		boolean indexcheck = candilist.contains(idx);
		long endtime1 = System.nanoTime();
		runrecord.addIOTime((endtime1-startTime1)/1000000000.0);
		return indexcheck;
	}
	
	/*
	 * compute the distance between any two centers for bound computation
	 */
	public void computeInterCentorid(int k, ArrayList<ClusterPath> Center, Map<Integer, int[]> clustData) {
		for(int i=0; i<k; i++) {
			innerCentoridDis[i] = new double[k];
			int []a = clustData.get(i);		
			double min = Double.MAX_VALUE;
			for(int j=0; j<k; j++) {				
				if(i!=j) {
					int []b = clustData.get(j);
					double distance = (double)Intersection(a, b, a.length, b.length);
					innerCentoridDis[i][j] = distance;
					if(distance<min) {
						min = distance;
					}					
				}
			}
			interMinimumCentoridDis[i] = min;
		}
		for (int i = 0; i < k; i++) {
			innerCentoridDisGroup[i] = new double[group.size()];
			for (int groupid : group.keySet()) {
				ArrayList<Integer> arrayList = group.get(groupid);
				double min = Double.MAX_VALUE;
				for (int centerid : arrayList) {
					if (innerCentoridDis[i][centerid] < min) {
						min = innerCentoridDis[i][centerid];
					}
				}
				innerCentoridDisGroup[i][groupid] = min;
			}
		}
	}
	
	/*
	 *  assign based on previous center to save time on IO and computation, group can be eliminated
	 */
	public void assignByTriangleFeaturesGroup(int k, int groupNumber) {
		Set<Integer> candidateofAllclusters = new HashSet<Integer>();
		Map<Integer, int[]> clustData = new HashMap<Integer, int[]>();
		Map<Integer, ArrayList<Integer>> idxNeedsIn = new HashMap<>();//it stores all the idxs of trajectories that move in
		Map<Integer, ArrayList<Integer>> idxNeedsOut = new HashMap<>();
		for (int j = 0; j < k; j++) {// combine the inverted index for prunning
			int []clustra = CENTERS.get(j).getTrajectoryData();
			clustData.put(j, clustra);
			long startTime1 = System.nanoTime();
			if(invertedIndex) {
				Set<Integer> candilist = CENTERS.get(j).creatCandidateList(edgeIndex, datamap);//generate the candidate list			
				Collections.addAll(candidateofAllclusters, candilist.toArray(new Integer[0]));
			}
			long endtime1 = System.nanoTime();
			runrecord.addIOTime((endtime1-startTime1)/1000000000.0);
		}
		double sumDist = 0;
		if(distanceModel==0)
			computeInterCentorid(k, CENTERS, clustData);//compute the inter centroid bound martix
		numeMovedTrajectories= 0;
		for (int group_i = 0; group_i < groupNumber; group_i++) {//check each group
			ArrayList<Integer> centers = group.get(group_i);//get the belonging 
			for (int centerID:centers){//check each center in the group
				Set<Integer> tralist = CENTERS.get(centerID).getClusterTrajectories();				
				for (int idx : tralist) { // check every trajectory in the center to assign which integrate the group filtering and local filtering				
					//check whether this trajectory is a pivot, use the radius to see whether the whole node can be assigned to a node directly
					int newCenterId = centerID;//initialize as the original center
					int centroid = CENTERS.get(centerID).getTrajectoryID();
					if(idx == centroid)// the centroid itself.
						continue;
					int[] tra = null;
					int tralength = traLength.get(idx); // the length of trajectory is read
					double min_dist = Double.MAX_VALUE;// to record the best center's distance					
					if (invertedIndex && !checkInvertedIndex(candidateofAllclusters, idx)) { // if it is never contained by any list, we can assign it to the cluster with minimum length						
						indexFil += k;
						for(int group_j=0; group_j<groupNumber; group_j++) {
							int[] clustra = clustData.get(group_j);
							double dist = Math.max(tralength, clustra.length);
							if (min_dist > dist) {
								min_dist = dist; 
								newCenterId = group_j;
							}
						}
					} else {//check whether we need to change the center by comparing the bounds															
						Set<Integer> canlist = CENTERS.get(centerID).getCandidateList();	
						int[] clustra = clustData.get(centerID);
						if(invertedIndex) {
							if (checkInvertedIndex(canlist, idx)) {
								min_dist = computeRealDistance(tra, clustra, idx);//compute the distance with new center						
							}else {// do not need to read as no overlap
								indexFil++;
								min_dist = Math.max(tralength, clustra.length);	
							}
						}else {
							min_dist = computeRealDistance(tra, clustra, idx);//compute the distance with new center	
						}
						if(distanceModel == 0 && assBoundSign) {//check the bound one by one							
						double [] bounds = trajectoryBounds.get(idx);
						double lowerbound = getMinimumLowerbound(bounds, groupNumber);	// bound from drift		
						double newupperbound = min_dist;// tighten the upper bound
						newCenterId = centerID;
						double centroidBound = interMinimumCentoridDis[centerID]/2.0;
						lowerbound = Math.max(lowerbound, centroidBound);//global bounds
						if(lowerbound < newupperbound){	//cannot not pass the global filtering
							for(int group_j=0; group_j<groupNumber; group_j++) {
								if( group_j == group_i)	//skip current group
									continue;
								double localbound = Math.max((bounds[group_j+2]-group_drift[group_j]), innerCentoridDisGroup[centerID][group_j]/2.0);
								if( localbound < min_dist) {//cannot pass the group filtering of bound 							
									ArrayList<Integer> centerCandidates = group.get(group_j);
									double second_min_dist_local = Double.MAX_VALUE;
									for(int center_j: centerCandidates) {// goto the local filtering on center in a group, by checking the candidate list and bounds												
										if(innerCentoridDis[centerID][center_j]/2.0 > min_dist) {//pass the inner centroid bound prunning
											numFillocal++;// local pruning
											continue;
										}
										canlist = CENTERS.get(center_j).getCandidateList();// get the candidate list of each cluster										
										clustra = clustData.get(center_j);
										double dist = 0;
										if(invertedIndex) {
											if (checkInvertedIndex(canlist, idx)) {
												dist = computeRealDistance(tra, clustra, idx);
											} else {
												indexFil++;
												dist = Math.max(tralength, clustra.length);
											}
										}else {
											dist = computeRealDistance(tra, clustra, idx);
										}
										if (min_dist > dist) {
											min_dist = dist; // maintain the one with min distance, and second min distance
											newCenterId = center_j;
										}
										if(second_min_dist_local>dist) {
											second_min_dist_local = dist;
										}
									}
									updateSingleLowerBound(trajectoryBounds, idx, group_j, second_min_dist_local);
								}else {
									numFilGroup += group.get(group_j).size();//pruned 
									updateSingleLowerBound(trajectoryBounds, idx, group_j, bounds[group_j+2] - group_drift[group_j]);
								}
							}
						}else {
							numFilGlobal += k;//k centroids are all pruned
						}
					}else {//brute force if there is no index
						min_dist = Double.MAX_VALUE;
						for(int group_j=0; group_j<groupNumber; group_j++) {
							clustra = clustData.get(group_j);
							double dist = computeRealDistance(tra, clustra, idx);
							if (min_dist > dist) {
								min_dist = dist; // maintain the one with min distance, and second min distance
								newCenterId = group_j;
							}
						}							
					}
					}	
					sumDist += min_dist;
					if(newCenterId != centerID) {// the trajectory moves to other center, this should be counted into the time of refinement.
						numeMovedTrajectories++;
						long Time1 = System.nanoTime();		
						ArrayList<Integer> idxlist;
						if(idxNeedsIn.containsKey(newCenterId))
							idxlist = idxNeedsIn.get(newCenterId);
						else
							idxlist = new ArrayList<Integer>();
						idxlist.add(idx);
						idxNeedsIn.put(newCenterId, idxlist);// temporal store as we cannot add them the trajectory list which will be scanned later, batch remove later
						if(idxNeedsOut.containsKey(centerID))
							idxlist = idxNeedsOut.get(centerID);
						else
							idxlist = new ArrayList<Integer>();
						idxlist.add(idx);
						idxNeedsOut.put(centerID, idxlist);// temporal store, batch remove later
						if(distanceModel == 0)
							accumulateHistogramGuava(tra, idx, newCenterId, centerID);	// update the histogram directly
						long Time2 = System.nanoTime();
						runrecord.addHistorgramTime((Time2-Time1)/1000000000.0);
					}
				}
			}
		}
		System.out.println("the sum distance after assignment is: "+sumDist+", #moved trajectories "+numeMovedTrajectories);
		long Time1 = System.nanoTime();		
		updateCenters(idxNeedsIn, idxNeedsOut);
		long Time2 = System.nanoTime();
		runrecord.addHistorgramTime((Time2-Time1)/1000000000.0);
	}
	
	/*
	 * the framework using Lloyd's algorithm 
	 */
	public int runKPath(int k, String folder, Set<Integer> candidateset) throws IOException {
		int groupNumber = k;
		if(k>50)
			groupNumber = k/10;
        trajectoryBounds = new HashMap<>();
        center_drift = new HashMap<Integer, Double>();
		groupInitialClusters(groupNumber, k); //	Step 1: divide k centroid into t groups
		singleKpath(k, 0, true, groupNumber, folder, candidateset); // 	Step 2, generate the initial center
		computeDrift(k, groupNumber);// Step 3.1 compute the drift using PRE_CENS and CENTERS
		int t = 1;
		double minimum = Double.MAX_VALUE;
		for(; t < TRY_TIMES; t++){
			if(graphPathExtraction){
				printCluterTrajectory(k, t, folder);
			}else {
				printCluterTraID(k, t, folder);
			}
			long startTime1 = System.nanoTime();		
			assignByTriangleFeaturesGroup(k, groupNumber);
			long endtime = System.nanoTime();
			runrecord.addAssignmentTime((endtime-startTime1)/1000000000.0);
			System.out.print("assign time cost: ");
			System.out.printf("%.3f", (endtime-startTime1)/1000000000.0);
			System.out.println("s");
			long startTime = System.nanoTime();
			double overallDis = 0;
			for(int i=0; i<k; i++) {
				double drfit = 0;
				if(distanceModel==0 && refinementSign) {
					if(graphPathExtraction){
						drfit = CENTERS.get(i).extractNewPathCPEP(forwardGraph, backwardGraph, i);// test the optimal
					}else {
						drfit = CENTERS.get(i).extractNewPathGuava(datamap, runrecord, traLength, trajectoryHistogram); //update the centroid of each cluster
					}
					center_drift.put(i, drfit);
				}else {
					CENTERS.get(i).extractNewPathStupid(datamap, runrecord, distanceModel, idPoints);
				}
				overallDis += CENTERS.get(i).getSumDistance();
			}
			if(distanceModel == 0) {
				computeDrift(k, groupNumber);// 	Step 3.1 compute the drift using PRE_CENS and CENTERS
			}
			endtime = System.nanoTime();
			runrecord.addRefinementTime((endtime-startTime)/1000000000.0);			
			System.out.print("iteration "+(t+1)+", the sum distance is "+overallDis+", time cost: ");
			if(overallDis<minimum)
				minimum = overallDis;
			System.out.printf("%.3f", (endtime-startTime1)/1000000000.0);
			System.out.println("s\n");
			if(timeToEnd()) {//all center does not change any more
				runrecord.setIterationtimes(t+1);
			//	t = TRY_TIMES+1;
				break;//convergence
			}
		}		
		System.out.println("\n#filtered trajectories by bound, inverted index: "+(numFilGroup+numFilGlobal)+", "+ indexFil);
		System.out.println("#moved trajectories: "+numeMovedTrajectories);
		System.out.println("#computation: "+numCompute);
	//	System.err.println(minimum);
		return t;
	}
	
	
	/*
	 * use the k-paths to build the index, when a group has a radius less than the threshold use the capacity as k to run kpath
	 */
	public void runIndexbuildQueue(int radius, int capacity) throws IOException {
		k = capacity;
		graphPathExtraction = false;
		mtreebuild = false;
		second = true;
		assBoundSign = true;
		refinementSign = true;
		invertedIndex = true;
		Initialization();
		GroupedTrajectory = new HashSet<>();
		loadData(datafile, trajectoryNumber, edgefile);	// load the data and create index
		createTrajectoryHistogram(datamap, trajectoryNumber);  // build inverted index if there is no index		
		initializeClustersRandom(k);
		groupFilename = datafile+"_"+Integer.toString(trajectoryNumber)+"_"+Integer.toString(radius)+"_"+Integer.toString(capacity)+"_index";
		String pivotname = groupFilename+".all";
		pivotGroup = new HashMap<>();
		Queue<Set<Integer>> queue = new LinkedList<>();
		queue.add(datamap.keySet());
		int numIteration = groupIteration;
		int firstIteration = 0;
		while(!queue.isEmpty()) {
			Set<Integer> candidates = queue.poll();
			CENTERS = new ArrayList<ClusterPath>();	
			interMinimumCentoridDis = new double[k];
			innerCentoridDis = new double[k][];					
			initializeClustersRandom(k, candidates);// initialize the center		
			System.out.print(candidates.size()+":");
			runKPath(capacity, folder, candidates);// run k path to divide into k groups
			String content = "";
			for(int i = 0; i<capacity; i++) {
				ClusterPath node = CENTERS.get(i);
				int nodeCapacity = node.getClusterTrajectories().size();
				System.out.print(nodeCapacity+";");
				content += Integer.toString(node.getTrajectoryID()) + ",";
			}
			System.out.println();

			for(int i = 0; i<capacity; i++) {
				ClusterPath node = CENTERS.get(i);
				int pivot = node.getTrajectoryID();//father
				candidates = node.getClusterTrajectories();
				int nodeCapacity = candidates.size();
				
				int nodeRadius = getRadius(candidates, pivot);// too slow, this is not necessary
				
				content = nodeRadius + ":" + pivot + ",";
				ArrayList<Integer> list = new ArrayList<>();
				list.add(pivot);// the first element
				for (int idx : candidates) {
					if (idx != pivot) {
						list.add(idx);
						content += Integer.toString(idx) + ",";
					}
				}
				if(nodeCapacity <= capacity || nodeRadius <= radius ) {// leaf node
					if(firstIteration == 0 && !content.equals("0:0,"))
						Util.write(pivotname, content+"\n");//write all the contents into the pivot table
					if (nodeCapacity >= capacity/2 && nodeRadius <= radius) {
						Util.write(groupFilename, content + "\n");// write the group into file
					}else if(nodeCapacity>0) {					
						GroupedTrajectory.addAll(candidates);// add to another file
					}
				}else {
					queue.add(candidates);// conduct the iteration again
				}	
			}
			if(queue.isEmpty() && --numIteration > 0 && !GroupedTrajectory.isEmpty()) {// conduct another round of pruning.
				firstIteration = 1;
				queue.add(GroupedTrajectory);
				GroupedTrajectory = new HashSet<>();
				System.err.println();
			}
		}
		for(int idx: GroupedTrajectory) {
			Util.write(groupFilename+".single", idx+"\n");
		}
	}
	
	public int [] generateSeeds(int kMax, int traNum) {
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
	 * range query, 
	 */
	public void runRangeClustering() {
		// input a range or a path,
		// get all the candidates in the path
		// run the clustering same to our method without using pivot table
	}
	
	/*
	 * run experiments to compare different optimizations by increasing k and data scale
	 */
	public void runExperiments(String[] args) throws IOException {
		mindex.buildGraph(graphfile, forwardGraph, backwardGraph);// build the graph in the begining
		PrintStream originalOut = System.out;
		for(int datascale = 100; datascale <= Integer.valueOf(args[2]); datascale*=10) {//increase the data scale
			if(args[5].equals("tdrive") && datascale==1000000)
				datascale = 250000;
			if(args[5].equals("porto") && datascale==10000000)
				datascale = 1500000;
			System.setOut(originalOut);
			trajectoryNumber = datascale;
			mtreebuild = false;
			second = true;// enable the several running in the same index
			Initialization();
			String indexfilename = "./logs/efficiency/"+args[5]+"_0_"+datascale+"_index";
			PrintStream indexfileOut = new PrintStream(indexfilename);
			System.setOut(indexfileOut);//record the data loading and index building time.
			
			loadData(datafile, trajectoryNumber, edgefile);	// load the data and create index
			String indexfile = datafile+"_"+Integer.toString(trajectoryNumber)+"_"+Integer.toString(10)+"_"+Integer.toString(20)+"_index.all";
			readIndex(indexfile);// read the file to build index.
			
			for(int expTimes = 0; expTimes< 10; expTimes++) {// run the queries set			
				String groupname = "./seeds/"+expTimes+"_"+datascale;//set the same random seeds first, multiple times
				indexfileOut = new PrintStream(groupname);
				System.setOut(indexfileOut);//record the data loading and index building time.			
				this.centroids = generateSeeds(100, datascale);//generate full seeds and write to the file
				for( k = 10; k < 60; k+=10) {//test different k
					CENTERS = new ArrayList<ClusterPath>();	
					initializeClustersRandom(k);// initialize the center
					interMinimumCentoridDis = new double[k];
					innerCentoridDis = new double[k][];				
					innerCentoridDisGroup = new double[k][];
					assBoundSign = false;
					refinementSign = false;
					invertedIndex = false;
					pivotPruning = false;
					graphPathExtraction = false;
					graphPathExtractionSimpliedRoad = false;
				
					
					// run the llyod's algorithm without any optimizations
					String filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_Base";
					PrintStream fileOut = new PrintStream(filename);
					System.setOut(fileOut);	
					if(datascale<=10000) {// limited by the distance matrix
					runKPath(k, folder, datamap.keySet());// run the baseline
					runrecord.printLog();
					concludeCenter();
					runrecord.clear();
					clearStatistics();
					
					// run the assignment with the lower bounds based on triangle inequality
					CENTERS = new ArrayList<ClusterPath>();	
					initializeClustersRandom(k);// initialize the center, output the centroids into the file
					interMinimumCentoridDis = new double[k];
					innerCentoridDis = new double[k][];
					innerCentoridDisGroup = new double[k][];
					assBoundSign = true;
					filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_ABound";
					fileOut = new PrintStream(filename);
					System.setOut(fileOut);	
					runKPath(k, folder, datamap.keySet());// run the baseline
					runrecord.printLog();
					concludeCenter();
					runrecord.clear();
					clearStatistics();
					}
					
					// run the refinement with the histogram
					CENTERS = new ArrayList<ClusterPath>();	
					initializeClustersRandom(k);// initialize the center
					interMinimumCentoridDis = new double[k];
					innerCentoridDis = new double[k][];
					innerCentoridDisGroup = new double[k][];
					refinementSign = true;
					assBoundSign = true;
					filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_RHis";//use the bound
					fileOut = new PrintStream(filename);
					System.setOut(fileOut);	
					runKPath(k, folder, datamap.keySet());// run the baseline
					runrecord.printLog();
					concludeCenter();
					runrecord.clear();
					clearStatistics();
					
					// run the assignment with the inverted index				
					CENTERS = new ArrayList<ClusterPath>();	
					initializeClustersRandom(k);// initialize the center
					interMinimumCentoridDis = new double[k];
					innerCentoridDis = new double[k][];
					innerCentoridDisGroup = new double[k][];
					invertedIndex = true;
					filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_AII";
					fileOut = new PrintStream(filename);
					System.setOut(fileOut);	
					runKPath(k, folder, datamap.keySet());// run the baseline
					runrecord.printLog();
					concludeCenter();
					runrecord.clear();
					clearStatistics();
					
					// run the assignment with the pivot table								
					CENTERS = new ArrayList<ClusterPath>();
					initializeClustersRandom(k);// initialize the center	
					interMinimumCentoridDis = new double[k];
					innerCentoridDis = new double[k][];
					innerCentoridDisGroup = new double[k][];
					pivotPruning = true;// use the pivot table in the assignment
					filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_APT";
					fileOut = new PrintStream(filename);
					System.setOut(fileOut);
					runKPath(k, folder, datamap.keySet());// run the baseline
					runrecord.printLog();
					concludeCenter();
					runrecord.clear();
					clearStatistics();
					
					
					// run the graph search in the refinement
					for(int iter = 1000; iter <= 5000; iter+=1000) {						
						CENTERS = new ArrayList<ClusterPath>();	
						initializeClustersRandom(k);// initialize the center
						for(int t=0; t<k; t++) {//test different iterations to show the performance
							ClusterPath cl = CENTERS.get(t);
							cl.setIteration(iter);
						}
						interMinimumCentoridDis = new double[k];
						innerCentoridDis = new double[k][];
						innerCentoridDisGroup = new double[k][];
						graphPathExtraction = true;
						filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_RGraph"+"_"+iter;
						fileOut = new PrintStream(filename);
						System.setOut(fileOut);	
						runKPath(k, folder, datamap.keySet());// run the baseline
						runrecord.printLog();
						concludeCenter();
						runrecord.clear();
						clearStatistics();
						//output the result to a new file.
					}
					
					
					/*
					 * run the simplied graph search in the refinement, test different iterations
					 * this method will not be tested 
					 */
				/*	filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_RGraphS";
					fileOut = new PrintStream(filename);
					System.setOut(fileOut);	
					CENTERS = new ArrayList<ClusterPath>();	
					graphPathExtractionSimpliedRoad = true;
					initializeClustersRandom1(k, true, nonImportantRoad);// initialize the center
					interMinimumCentoridDis = new double[k];
					innerCentoridDis = new double[k][];			
					innerCentoridDisGroup = new double[k][];
					runKPath(k, folder, datamap.keySet());// run the baseline
					runrecord.printLog();
					concludeCenter();
					runrecord.clear();
					clearStatistics();*/
				}
			}
		}
	}
	
	/*
	 * run experiments to compare different optimizations by increasing k and data scale
	 */
	public void runExperimentsElbow(String[] args) throws IOException {
		mindex.buildGraph(graphfile, forwardGraph, backwardGraph);// build the graph in the begining
		PrintStream originalOut = System.out;
		for(int datascale = Integer.valueOf(args[2]); datascale <= Integer.valueOf(args[2]); datascale*=10) {//increase the data scale
			if(args[5].equals("tdrive") && datascale==1000000)
				datascale = 250000;
			if(args[5].equals("porto") && datascale==10000000)
				datascale = 1500000;
			System.setOut(originalOut);
			trajectoryNumber = datascale;
			mtreebuild = false;
			second = true;// enable the several running in the same index
			Initialization();
			String indexfilename = "./logs/efficiency/"+args[5]+"_0_"+datascale+"_index";
			PrintStream indexfileOut = new PrintStream(indexfilename);
			System.setOut(indexfileOut);//record the data loading and index building time.
			loadData(datafile, trajectoryNumber, edgefile);	// load the data and create index
			String indexfile = datafile+"_"+Integer.toString(trajectoryNumber)+"_"+Integer.toString(10)+"_"+Integer.toString(20)+"_index.all";			
			for(int expTimes = 0; expTimes< 1; expTimes++) {// run the queries set			
				String groupname = "./seeds/"+expTimes+"_"+datascale;//set the same random seeds first, multiple times
				indexfileOut = new PrintStream(groupname);
				System.setOut(indexfileOut);//record the data loading and index building time.			
				this.centroids = generateSeeds(100, datascale);//generate full seeds and write to the file
				for( k = 2; k < 100; k++) {//test different k
					CENTERS = new ArrayList<ClusterPath>();	
					initializeClustersRandom(k);// initialize the center
					interMinimumCentoridDis = new double[k];
					innerCentoridDis = new double[k][];				
					innerCentoridDisGroup = new double[k][];
					assBoundSign = true;
					refinementSign = true;
					invertedIndex = true;
					String filename = "./logs/efficiency/"+args[5]+"_0_"+k+"_"+datascale+"_"+expTimes+"_AII";
					PrintStream fileOut = new PrintStream(filename);
					System.setOut(fileOut);	
					runKPath(k, folder, datamap.keySet());// run the baseline		
					// compute the cvi, sum distance, and draw the figures.					
					runrecord.printLog();
					concludeCenter();
					runrecord.clear();
					clearStatistics();
				}
			}
		}
	}
	
	
	
	/*
	 * streaming kpaths processing based on a given dataset and inverted index in a sliding window.
	 */
	public int[] streamingClustering(Map<Integer, int[]> trajectories, 
			Map<Integer, Set<Integer>> edgeInvertedList, String args[], int counter) throws IOException{
    //    if(CENTERS==null)
        	CENTERS = new ArrayList<>();// renew the center

		this.datamap = trajectories;
        trajectoryNumber = datamap.size();
        this.edgeIndex = edgeInvertedList;
        traLength = new HashMap<>();
        initializeClustersRandomStreaming(k);// initialize the center
        for(int idx: datamap.keySet())
        	traLength.put(idx, datamap.get(idx).length);
        
        k = Integer.valueOf(args[1]);
        refinementSign = true;
		assBoundSign = true;
		invertedIndex = true;
		pivotPruning = false;
		graphPathExtraction = false;// we use the centroid

		interMinimumCentoridDis = new double[k];
		innerCentoridDis = new double[k][];
		innerCentoridDisGroup = new double[k][];
		runKPath(k, folder, datamap.keySet());// run the baseline
		runrecord.printLog();
		runrecord.clear();
		clearStatistics();
        return null;
	}
	
	/*
	 * read the Pivot table index from the file we have built
	 */
	void readIndex(String filename) {
		GroupedTrajectory = new HashSet<>();
		pivotGroup = new HashMap<>();
		pivotRadius = new HashMap<>();
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(filename)));			
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split(":");
				double radius = Double.valueOf(abc[0]);
				String[] traids = abc[1].split(",");
				int pivot = Integer.valueOf(traids[0]);
				GroupedTrajectory.add(pivot);
				if (traids.length > 1) {
					ArrayList<Integer> arrayList = new ArrayList<>();
					for (int idString=1; idString<traids.length; idString++) {
						arrayList.add(Integer.valueOf(traids[idString]));
						GroupedTrajectory.add(Integer.valueOf(traids[idString]));
					}
					pivotGroup.put(pivot, arrayList);
				} else {
					pivotGroup.put(pivot, null);
				}
				pivotRadius.put(pivot, radius);
			}
			in.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * scan every trajectory to get the maximum distance as the radius.
	 */
	int getRadius(Set<Integer> candidates, int pivot) {
		int max = 0;
		int []centroid = datamap.get(pivot);
		for(int idx: candidates) {
			int[] data = datamap.get(idx);
			int dis = Intersection(centroid, data, centroid.length, data.length);
			if(dis > max) {
				max = dis;
			}			
		}
		return max;
	}
	
	/*
	 * load the dataset for testing on the new distance function EBD
	 */
	public void loadDataOnly(boolean path) throws IOException {
		graphPathExtraction = path;
		mtreebuild = false;
		second = true;		
		Initialization();
		loadData(datafile, trajectoryNumber, edgefile);	// load the data and create index
	}
	
	public void printResults() {
		for(ClusterPath a: CENTERS) {
			a.printClusterTrajectories();
		}
	}
	
	/*
	 * run single test
	 */
	public void staticKpath(boolean path) throws IOException {
		graphPathExtraction = path;
		mtreebuild = false;
		second = true;		
		Initialization();
		invertedIndex = true;
		refinementSign = true;
		assBoundSign = true;
		graphPathExtraction = false;// use the graph
		loadData(datafile, trajectoryNumber, edgefile);	// load the data and create index
		if(second==true)
			createTrajectoryHistogram(datamap, trajectoryNumber);  // build inverted index if there is not index
        String new_folder = mapv_path + folder;
    //    boolean creat_folder = new File(new_folder).mkdirs();
    //  if(!creat_folder)
     //   	return;      
        mindex.buildGraph(graphfile, forwardGraph, backwardGraph);
        initializeClustersRandom(k); 	 //randomly choose k
    //    initializeClusterskMeansPlusPlus(k);// test kmeans++;
        int iterations = runKPath(k, folder, datamap.keySet());
		runrecord.printLog();
		int []result;
		if(graphPathExtraction){
			System.out.println("\nThe final centroid trajectory length are:");
			printCluterTrajectory(k, iterations+1, folder);
		}else {
			System.out.println("\nThe final centroid trajectory ids are:");
			result = printCluterTraID(k, iterations+1, folder);
		}
		int cou= 0;
		for(int edgeid :  edgeHistogram.keySet()) {
			cou +=edgeHistogram.get(edgeid);
		}
		frequencyThreshold = cou/edgeHistogram.size();
		SIGMOD07.SIGMOD07FrequentEdgesMapv(edgeHistogram, edgeInfo, frequencyThreshold, Integer.toString(trajectoryNumber)+".sigmod2018");		
		ArrayList<int []> clusterids1 = new ArrayList<int []>();
		ArrayList<Integer> clusterids = new ArrayList<Integer>();
		for(int t = 0; t < k; t++) {
			int []a = CENTERS.get(t).getTrajectoryData();
			for(int ids :a)
				clusterids.add(ids);
			clusterids1.add(datamapunsorted.get(CENTERS.get(t).getTrajectoryID()));
		}
		mapv.generateHighEdges(edgeInfo, clusterids, Integer.toString(trajectoryNumber)+".kpaths");
		mapv.generateClusterPath1(datamap, edgeInfo, clusterids1, Integer.toString(trajectoryNumber)+".kpaths_color");
		
		groupsBuilt = new HashMap<>();		//store the centroid as group
		for(int i=0; i<k; i++) {
			ArrayList<Integer> arrayList = new ArrayList<Integer>(CENTERS.get(i).getClusterTrajectories());
			groupsBuilt.put(i, arrayList);
		}
	}
	
	/*
	 * conduct the shortest path search
	 */
	public void ShortestPathSearch(Map<Integer, String> edgeInfo, int vertexnumber) {
		// TODO Auto-generated constructor stub
		Map<String, Integer> vertexinfo = new HashMap<>();
		DefaultDirectedGraph<Integer, DefaultEdge> directedGraph =
	            new DefaultDirectedGraph<Integer, DefaultEdge>(DefaultEdge.class);		
		for(int i=0; i<vertexnumber; i++) {//add all the vertex
			directedGraph.addVertex(i);
		}
		
		for(int edge: edgeInfo.keySet()) {//add all the edges
			String vertexs = edgeInfo.get(edge);
			String vertex[] = vertexs.split(",");
			vertexinfo.put(vertexs, edge);
			directedGraph.addEdge(Integer.valueOf(vertex[0]), Integer.valueOf(vertex[1]));
		}
		
		// Prints the shortest path from vertex i to vertex c. This certainly
        // exists for our particular directed graph.
        System.out.println("Shortest path from i to c:");
        int i=0, c= 1;
        DijkstraShortestPath<Integer, DefaultEdge> dijkstraAlg =
            new DijkstraShortestPath<>(directedGraph);
        SingleSourcePaths<Integer, DefaultEdge> iPaths = dijkstraAlg.getPaths(i);
        System.out.println(iPaths.getPath(c) + "\n");
        GraphPath<Integer, DefaultEdge> aDefaultEdge = iPaths.getPath(c);
        //translate to a set of edges, and further compute the distance.
	}
	
	/*
	 * run single test
	 */
	public void staticKpath(boolean path, int num) throws IOException {
		graphPathExtraction = path;
		mtreebuild = false;
		second = true;		
		Initialization();
		invertedIndex = true;
		refinementSign = true;
		assBoundSign = true;
	//	graphPathExtraction = false;// use the graph
		loadData(datafile, trajectoryNumber, edgefile);	// load the data and create index
		
		if(second==true)
			createTrajectoryHistogram(datamap, trajectoryNumber);  // build inverted index if there is not index     
        mindex.buildGraph(graphfile, forwardGraph, backwardGraph);
        initializeClustersRandom(k); 	 //randomly choose k
    //    initializeClusterskMeansPlusPlus(k);// test kmeans++;
        int iterations = runKPath(k, folder, datamap.keySet());
		runrecord.printLog();
		int []result;
		if(graphPathExtraction){
			System.out.println("\nThe final centroid trajectory length are:");
			printCluterTrajectory(k, iterations+1, folder);
		}else {
			System.out.println("\nThe final centroid trajectory ids are:");
			result = printCluterTraID(k, iterations+1, folder);
		}
		int cou= 0;
		for(int edgeid :  edgeHistogram.keySet()) {
			cou +=edgeHistogram.get(edgeid);
		}
		frequencyThreshold = cou/edgeHistogram.size();
	//	SIGMOD07.SIGMOD07FrequentEdgesMapv(edgeHistogram, edgeInfo, frequencyThreshold, Integer.toString(trajectoryNumber)+".sigmod2018");		
	//	SIGMOD07.topkfrequentedegs(edgeHistogram, edgeInfo, 10, Integer.toString(trajectoryNumber)+".10frequentEdges");
		SIGMOD07.topkTrajectoryWithHighestFrequencySum(edgeHistogram, edgeInfo, 10, Integer.toString(trajectoryNumber), datamapunsorted);
		ArrayList<int []> clusterids1 = new ArrayList<int []>();
		ArrayList<Integer> clusterids = new ArrayList<Integer>();
		for(int t = 0; t < k; t++) {
			int []a = CENTERS.get(t).getTrajectoryData();
			for(int ids :a)
				clusterids.add(ids);
			if(!graphPathExtraction)
				clusterids1.add(datamapunsorted.get(CENTERS.get(t).getTrajectoryID()));
			else {
				clusterids1.add(CENTERS.get(t).getTrajectoryDataUnsorted());//the graph-based
			//	System.out.println(Arrays.toString(CENTERS.get(t).getTrajectoryDataUnsorted()));
			}
		}
	//	mapv.generateHighEdges(edgeInfo, clusterids, Integer.toString(trajectoryNumber)+".kpaths");
		if(!graphPathExtraction)
			mapv.generateClusterPath1(datamap, edgeInfo, clusterids1, Integer.toString(trajectoryNumber)+".kpaths_color"+num);
		else
			mapv.generateClusterPath1(datamap, edgeInfo, clusterids1, Integer.toString(trajectoryNumber)+".kpaths_color"+num+"_graph");
	/*	groupsBuilt = new HashMap<>();		//store the centroid as group
		for(int i=0; i<k; i++) {
			ArrayList<Integer> arrayList = new ArrayList<Integer>(CENTERS.get(i).getClusterTrajectories());
			groupsBuilt.put(i, arrayList);
		}*/
	}
}
