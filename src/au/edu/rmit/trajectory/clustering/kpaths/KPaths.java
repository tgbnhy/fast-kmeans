package au.edu.rmit.trajectory.clustering.kpaths;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import me.lemire.integercompression.*;
import me.lemire.integercompression.differential.*;
import com.kamikaze.pfordelta.*;
import au.edu.rmit.trajectory.clustering.TrajectoryMtree;
import scala.Int;
import scala.annotation.meta.getter;

/*
 * static kpath is used for clustering all the taxi trips, so we can find the k-paths for designing bus routes, this is not a real time problem, so we can 
 */
public class KPaths extends Thread {
	// stores the clusters
	protected ArrayList<ClusterPath> CENTERS; // it stores the k clusters
	ArrayList<ClusterPath> PRE_CENS; // it stores the previous k clusters
	
	// the parameters
	protected int TRY_TIMES = Integer.valueOf(LoadProperties.load("try_times"));//iteration times
	protected String mapv_path = LoadProperties.load("vis_path");
	String mapv_path_traclu_sigmod07 = LoadProperties.load("TraClus");
	int frequencyThreshold = Integer.valueOf(LoadProperties.load("frequencyThreshold"));
	int streamingDuration = Integer.valueOf(LoadProperties.load("streamingDuration"));
	int streamEdges = Integer.valueOf(LoadProperties.load("streamEdges"));
	
	protected RunLog runrecord = new RunLog(); 
	ArrayList<int[]> centroids;
	protected int trajectoryNumber;// the number of trajectories in the dataset
	protected String folder;
	protected int k;
	
	boolean dataEnough;
	boolean dataOut;
	boolean iterationStops;
	boolean readingdata;
	
	protected String datafile;
	int traNumber;
	String edgefile;
	protected String graphfile;
	int slidingwindow = 10;//a window to control the inverted index, 
	
	protected int distanceModel;//0 for EBD, 1 for EDR, 2 for DTW, 3 For Frechet, 4 for haursdoff, 5 for ERP
	//for Yinyang and bound computation
	protected  Map<Integer, double[]> trajectoryBounds;// build a lower bound list with all groups and upper bound for each trajectory, 0: upper bound, 1: global lower bound, 2~: lower bound with all different group	
	protected Map<Integer, ArrayList<Integer>> group;// group id, centers id belong to this group
	protected Map<Integer, Integer> centerGroup;//center id, group id	
	protected double[][] innerCentoridDis;//stores the distance between every two centorids
	protected double[][] innerCentoridDisGroup;//stores the distance between every two centorids
	protected double[] interMinimumCentoridDis;//store the distance to nearest neighbor of each centorid
	
	//for storage
	protected Map<Integer, int[]> datamap; // the trajectory dataset
	protected Map<Integer, int[]> datamapunsorted; // the trajectory dataset

	protected Map<Integer, Integer> traLength; // the trajectory dataset
	protected Map<Integer, Integer> trajectoryHistogram;//the histogram of each trajectory
	protected Map<Integer, Set<Integer>> edgeIndex;// the index used for similarity search
	protected Map<Integer, Integer> edgeHistogram;
	protected Map<Integer, String> edgeInfo;// the edge information
	Map<Integer, Integer> edgeType;
	protected Map<Integer, double[]> idPoints;
	
	//for graph
	protected HashMap<Integer, ArrayList<Integer>> forwardGraph = new HashMap<Integer, ArrayList<Integer>>();//the linked edge whose start is the end of start
	protected HashMap<Integer, ArrayList<Integer>> backwardGraph = new HashMap<Integer, ArrayList<Integer>>();//the linked edge whose start is the end of start
	protected ArrayList<int[]> centoridData = new ArrayList<>();//initialize the centroid
	HashMap<String, Integer> road_types;
	ArrayList<Integer> nonImportantRoad;//for the road network simplification
	
	//Mtree index 
	protected TrajectoryMtree mindex = new TrajectoryMtree();
	protected HashMap<Integer, Double> pivotRadius;// the radius
	protected HashMap<Integer, ArrayList<Integer>> pivotGroup;// the group for pivot index
	protected Set<Integer> GroupedTrajectory;
	
	protected boolean assBoundSign;// a sign used to set whether we use the bounds
	protected boolean refinementSign;// a sign used to set whether we use O(n) refinement
	protected boolean invertedIndex;// a sign used to set whether we use the inverted index
	protected boolean mtreebuild;// a sign used to set whether we use the pivot-table
	protected boolean pivotPruning;// a sign used to set whether we use the pivot-table
	protected boolean graphPathExtraction;// a sign used to set whether we use the optimization
	protected boolean graphPathExtractionSimpliedRoad;// a sign used to set whether we use the graph with a simple road network
	protected boolean Euclidean;
	protected boolean second = true; //
	
	int numFilCenter=0;
	protected int numFillocal=0;
	protected int numFilGroup=0;
	protected int numFilGlobal=0;
	protected int numeMovedTrajectories=0;
	protected int indexFil=0;
	protected int numCompute=0;
	
	double maxDistance[];
	
	public void clearStatistics() {
		numFilCenter=0;
		numFillocal=0;
		numFilGroup=0;
		numFilGlobal=0;
		numeMovedTrajectories=0;
		indexFil=0;
		numCompute=0;
	}
	public KPaths(String datapath) {
		trajectoryNumber=0;
	}
	
	public double[] getmaxDistance() {
		return maxDistance;
	}
	
	public KPaths(String []args) {
		datafile = args[0];
		k = Integer.valueOf(args[1]);
		trajectoryNumber = Integer.valueOf(args[2]);
		edgefile = args[3];
		graphfile = args[4];
		maxDistance = new double[6];
	}
	
	public ArrayList<ClusterPath> getCENTERS() {
		return CENTERS;
	}
	// stop the iteration when the clusters do not change compared with last time
	protected boolean timeToEnd() {
		for (ClusterPath cc : CENTERS) {
			if(cc.getCenterChanged()==true) {
				return false;
			}
		}
		return true;
	}

	/*
	 * load the data from file
	 */
	public void loadData(String path, int number, String edgePath) throws IOException{
		int idx=0;
		int gap = number/k;
		Random rand = new Random();
		int counter = 0;
		edgeIndex = new HashMap<>();
		edgeHistogram = new HashMap<>();
		traLength = new HashMap<>();
		datamap = new HashMap<>();
		datamapunsorted = new HashMap<>();
		int dataFootprint = 0;
		IntegratedIntCompressor iic = new IntegratedIntCompressor();
		int countlength = 0;
		if(second == true)//conduct the second kpath if it is false
			readRoadNetwork(edgePath);
		try {
			int[] roadFrequency = new int[road_types.size()];
			Scanner in = new Scanner(new BufferedReader(new FileReader(path)));			
			long time = 0;
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split("\t");
				String[] vertexSeries = abc[1].split(",");
				int[] vertexes = new int[vertexSeries.length];
				int[] vertexes1 = new int[vertexSeries.length];
				long timeStart = System.nanoTime();
				for(int t=0; t < vertexSeries.length; t++) {					
					vertexes[t] = Integer.valueOf(vertexSeries[t]);
					vertexes1[t] = Integer.valueOf(vertexSeries[t]);
					int edgeID = vertexes[t];
					if(second == true) {
						if(edgeIndex.containsKey(edgeID)) {
							Set<Integer> lists = edgeIndex.get(edgeID);
							lists.add(idx);	//enlarge the lists
							edgeIndex.put(edgeID, lists);
						}else {
							Set<Integer> lists = new HashSet<Integer>();
							lists.add(idx);	
							edgeIndex.put(edgeID, lists);
						}
					}
					/*	if(edgeType.containsKey(edgeID)) {
						int roadType = edgeType.get(edgeID);
						roadFrequency[roadType]++;
					}*/
					if(edgeHistogram.containsKey(edgeID)) {
						edgeHistogram.put(edgeID, edgeHistogram.get(edgeID)+1);
					}else {
						edgeHistogram.put(edgeID, 1);
					}
				}
			//	System.out.println(Arrays.toString(vertexes));
				datamapunsorted.put(idx, vertexes1);
				if(distanceModel==0) {//only sort for fast distance computation
					Arrays.sort(vertexes);// this sort the array
				//	int[] compressed = iic.compress(vertexes);
				//	dataFootprint += compressed.length;
				}
				mtreebuild = false;
				if(mtreebuild) {//build the mtree, will not be used anymore
					if(second == false) {
						mindex.buildMtree(vertexes, idx);//create the M-tree
						System.out.println(idx);
					}else if(idx==counter && centoridData.size()<k) {// initialize the centroid in the beginning
						centoridData.add(vertexes);
						counter += rand.nextInt(gap);
						ClusterPath cl = new ClusterPath(vertexes, 0);
						CENTERS.add(cl);
					}
					idx++;
				}else {
				//	System.err.println(idx);
					countlength+=vertexSeries.length;
					traLength.put(idx, vertexSeries.length);
					datamap.put(idx++, vertexes);
				}
				long timeEnd = System.nanoTime();
				time += timeEnd-timeStart;
				if(idx>number)
					break;
			}
			System.out.println("Data loading and inverted index building time:"+ time/1000000000.0 + "s");
			/*	System.out.println(Arrays.toString(roadFrequency));
			System.out.println(road_types.keySet().toString());
				TreeMap<Integer, String> freStreet = new TreeMap<>();
				int count=0;
			for(String road : road_types.keySet()) {
				freStreet.put(roadFrequency[count++], road);
			}
			System.out.println(freStreet.toString());*/
			in.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	//	System.err.println("the whole length is"+countlength);
		if(mtreebuild) {
			if(second == false) {
				System.out.println("the trajectory dataset is loaded");			
				System.out.println("the frequency histogram of edge is built");
				System.out.println("the inverted index of edge is built");
				System.out.println("the M-tree is built");
				mindex.buildHistogram();//build the histogram
				mindex.writeMtree(mapv_path);//write to the disk
			}
		}
	//	System.err.println("compress data is: "+dataFootprint);
	//	invertedIndex();
	}
	
	public void invertedIndex() {		
		int original = 0;
		int compress = 0;
		IntegratedIntCompressor iic = new IntegratedIntCompressor();
		for(int edgeid: edgeIndex.keySet()) {
			Set<Integer> arrayList = edgeIndex.get(edgeid);
			TreeSet<Integer> aIntegers = new TreeSet<>(arrayList);
			original += arrayList.size();
			int[] data = new int[arrayList.size()];
			int i=0;
			for(int traid: aIntegers) {
				data[i++] = traid;
			}
	        int[] compressed = iic.compress(data); // compressed array
	        compress += compressed.length;
		}
	//	System.err.println("compress index is: "+original+" "+compress);
	}
	public Map<Integer, int[]> getdatamap(){
		return datamap;
	}
	
	public Map<Integer, double[]> getPointcoordiate(){
		return idPoints;
	}
	
	public int getTraNum(){
		return trajectoryNumber;
	} 
	/*
	 * the edge information
	 */
	void readRoadNetwork(String edgePath) throws FileNotFoundException {
		road_types = new HashMap<>();
		edgeType = new HashMap<>();
		idPoints = new HashMap<>();
		nonImportantRoad = new ArrayList<>();
		int type=0;
		String filename = "./logs/edgelength";
		PrintStream fileOut = new PrintStream(filename);
		System.setErr(fileOut);
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(edgePath)));
			while (in.hasNextLine()) {	// load the geo-information of all the edges in the graph
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split(";");
				edgeInfo.put(Integer.valueOf(abc[0]), abc[1]+","+abc[2]);
				String xString = abc[1].split(",")[0];
				String yString = abc[2].split(",")[0];
				double xy[] = new double[2];
				xy[0] = Double.valueOf(xString);
				xy[1] = Double.valueOf(yString);
				idPoints.put(Integer.valueOf(abc[0]), xy);//this is for DTW, Frechet, Hausdorff, ERP
				
				if(abc.length>7) {
					 int roadType = 0;
					 if(!road_types.containsKey(abc[6])) {
						 road_types.put(abc[6], type);//we build the edge histogram
						 roadType = type++;
					 }
					 else{
						 roadType = road_types.get(abc[6]);//"primary road"
					 }
					 edgeType.put(Integer.valueOf(abc[0]), roadType);
				//	 if(checkNonImportantRoad(abc[6]))// add the non important road into arraylist
					//	 nonImportantRoad.add(Integer.valueOf(abc[0]));
				}
				System.err.println(abc[3]);
			}
			in.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	boolean checkNonImportantRoad(String edgeType) {
		if(edgeType == "" || edgeType == "unclassified" || edgeType == "tertiary_link" || edgeType == "primary_link" || 
				edgeType == "living_street" || edgeType == "trunk" || edgeType == "motorway_link" || edgeType == "motorway" ||
				edgeType == "secondary" || edgeType == "residential" || edgeType == "service" || edgeType == "pedestrian" || 
				edgeType == "trunk_link" || edgeType == "primary") {
			return false;
		}
		else
			return true;
	}
	
	/*
	 * read the types
	 */
	void concludeWholeDatasetRoadTypes() {
		
	}
	
	/*
	 * conclude the edge information of the final centroids, for concluding the road
	 */
	void concludeCenter() {
		int a[]  = new int[road_types.size()+1];
		for(int t = 0; t < k; t++) {
			int []trajectory = CENTERS.get(t).getTrajectoryData();
			for(int i=0; i<trajectory.length; i++) {
				if(edgeType.containsKey(trajectory[i])) {
					int edgeTypes = edgeType.get(trajectory[i]);
					a[edgeTypes]++;
				}else {
					a[road_types.size()]++;// the unknown one
				}
			}
		}
		System.out.println("there are "+(road_types.size()+1)+" types of roads");
		for(String id: road_types.keySet()) {
			if(a[road_types.get(id)]>0) {
				System.out.println(id+" "+road_types.get(id)+" "+a[road_types.get(id)]);
			}
		}
	}
	
	/*
	 * build inverted index on the edge
	 */
	public void createTrajectoryHistogram(Map<Integer, int[]> datamap, int trajectoryNumber) {	
	/*	trajectoryHistogram = new HashMap<>();
		for(int idx:datamap.keySet()) {	//scan each trajectory, compute the frequency for each trajectory
			int[] tra = datamap.get(idx);
			int tra_fre =0;
			for(int t=0; t<tra.length; t++) {	//scan each edge
				tra_fre += edgeHistogram.get(tra[t]); //the frequency is the sum of edge frequency in each trajectory.
			}
			trajectoryHistogram.put(idx, tra_fre);
		}	*/	
		System.out.println("==============================================================\n");
	}
	
	/*
	 * initialize the k clusters by randomly choosing from existing trajectories
	 */
	void initializeClustersRandom(int k) {
		Random rand = new Random();
		for(int t=0; t<k; t++) {
			int  n = rand.nextInt(trajectoryNumber) + 1;
			int[] cluster = datamap.get(n);
			ClusterPath cl = new ClusterPath(cluster, n);
			CENTERS.add(cl);		
		}
	}
	
	/*
	 * initialize the k clusters by randomly choosing from existing trajectories
	 */
	void initializeClustersRandom(int k, Set<Integer> data) {
		Random rand = new Random();
		ArrayList<Integer> arrayList = new ArrayList<>(data);
		for(int t=0; t<k; t++) {
			int  n = rand.nextInt(data.size());
			int idx = arrayList.get(n);
			int[] cluster = datamap.get(idx);
			ClusterPath cl = new ClusterPath(cluster, n);
			CENTERS.add(cl);		
		}
	}
	
	/*
	 * initialize the k clusters by choosing from existing trajectories incrementally
	 */
	protected void initializeClustersIncrease(int k, int delta) {
		int n = 0;
		for(int t=0; t<k; t++) {
			n += delta;
			int[] cluster = datamap.get(n);
			System.out.print(cluster.length+",");
			ClusterPath cl = new ClusterPath(cluster, n);
			CENTERS.add(cl);		
		}
		System.out.println();
	}
	
	/*
	 * initialize the k clusters by choosing from existing trajectories which have high frequency and do not intersect with each other
	 */
	void initializeClustersHighFrequency(int k, int range) {
		Random rand = new Random();
		//sort the trajectoryHistogram by value decreasingly, choose the top 1000 or more randomly.
		Map<Integer, Integer> sortedMap = trajectoryHistogram.entrySet().stream()
			    .sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
			    .collect(Collectors.toMap(Entry::getKey, Entry::getValue,
			                              (e1, e2) -> e1, LinkedHashMap::new));		
		ArrayList<Integer> keyset = new ArrayList<Integer>(sortedMap.keySet());
		ArrayList<Integer> allEdges = new ArrayList<Integer>();
		for(int t=0; t<k; t++) {
			int n = rand.nextInt(range) + 1;
			int idx = keyset.get(n);
			int[] cluster = datamap.get(idx);
			for(int i=0; i<cluster.length; i++) {
				if(allEdges.contains(i)) {
					t--;
					continue;
				}
			}
			Collections.addAll(allEdges, Arrays.stream(cluster).boxed().toArray(Integer[]::new));
			ClusterPath cl = new ClusterPath(cluster, idx);
			CENTERS.add(cl);
		}
	}
	
	/*
	 *  print the cluster ids and generate the trajectory into file for mapv visualization
	 */
	public int[] printCluterTraID(int k, int iteration, String folder) {		
		int []clusterid = new int[k];
		for(int t = 0; t < k; t++) {
			System.out.print(CENTERS.get(t).getTrajectoryID()+",");
			clusterid[t] = CENTERS.get(t).getTrajectoryID();
		}
		System.out.println();
		return clusterid;
	}
	
	/*
	 *  print the cluster ids and generate the trajectory into file for mapv visualization
	 */
	public void printCluterTrajectory(int k, int iteration, String folder) {		
		String output = mapv_path+folder+"\\"+Integer.toString(iteration);
		centroids = new ArrayList<>();
		for(int t = 0; t < k; t++) {
			System.out.print(CENTERS.get(t).getTrajectoryData().length+",");
			centroids.add(CENTERS.get(t).getTrajectoryData());
		}
		System.out.println();
	//	mapv.generateClusterPath1(datamap, edgeInfo, centroids, output);
	}

	/*
	 * assign by building the histogram again
	 */
	public ArrayList<ClusterPath> assignRebuildInvertedindex(int k, ArrayList<ClusterPath> new_CENTERS, boolean yinyang, 
			int groupnumber, Set<Integer> candidateset) {
		Set<Integer> candidateofAllclusters = new HashSet<Integer>();
		Map<Integer, int[]> clustData = new HashMap<Integer, int[]>();
		int minlength = Integer.MAX_VALUE;
		int min_length_id=0;
		double sumDist = 0;
		for (int j = 0; j < k; j++) {
			int[] clustra = CENTERS.get(j).getTrajectoryData();			
			clustData.put(j, clustra);
			if(invertedIndex) {
				Set<Integer> candilist = CENTERS.get(j).creatCandidateList(edgeIndex, datamap);//generate the candidate list
				candidateofAllclusters.addAll(candilist);// merge it to a single list
				if(clustra.length < minlength) {// get the minimum length
					minlength = clustra.length;
					min_length_id = j;
				}
			}
		}
	
		for (int idx:candidateset) {// first assignment, check whether we can assign in batch
			long Time1 = System.nanoTime();
			int[] tra = datamap.get(idx);//the trajectory data is read	
			long Time2 = System.nanoTime();
			runrecord.addIOTime((Time2-Time1)/1000000000.0);			
			double min_dist = Double.MAX_VALUE;
			int min_id = 0;		
			double [] bounds = null;
			if(yinyang) {//create the bounds for pruning in the first iterations.
				bounds = new double[groupnumber+2];
				Arrays.fill(bounds, Double.MAX_VALUE);// only check for the first rounds
			}			
			if(invertedIndex && !candidateofAllclusters.contains(idx)) {//if it is never contained by any list, we can assign it to the cluster with minimum length
				min_dist = Math.max(tra.length, minlength);
				min_id = min_length_id;
				if(yinyang) {// initialize the lower bound
					for (int j = 0; j < k; j++) {
						int length = clustData.get(j).length;
						double dist =  Math.max(tra.length, length);
						int groupNumber=centerGroup.get(j);
						if(j==min_length_id) {
							j=k;
						}
						if(dist<bounds[groupNumber+2]) {
							bounds[groupNumber+2] = dist;
						}
					}
				}
				indexFil+=k;
			}else {
				for (int j = 0; j < k; j++) {
					int centroid = CENTERS.get(j).getTrajectoryID();
					if(idx == centroid) {// the centroid itself.
						min_dist = 0;
						min_id = j;
					}
					Set<Integer> canlist = CENTERS.get(j).getCandidateList();// get the candidate list of each cluster
					double dist = 0;
					int[] clustra = clustData.get(j);
					if(invertedIndex) {
						if(!canlist.contains(idx)) {						// it is not contained
							indexFil++;
							dist = Math.max(tra.length, clustra.length);
						}else {			
							numCompute++;
							dist = (double)Intersection(tra, clustra, tra.length, clustra.length);
						}
					}else {
						numCompute++;
						if(distanceModel==0)
							dist = (double)Intersection(tra, clustra, tra.length, clustra.length);
						else if(distanceModel==1)
							dist = Util.EDR(tra, clustra);
						else if(distanceModel==2)
							dist = Util.DTW(tra, clustra, idPoints);
						else if(distanceModel==3)
							dist = Util.Frechet(tra, clustra, idPoints);
						else if (distanceModel == 4){
							dist = Util.Hausdorff(tra, clustra, idPoints);
						}else if (distanceModel == 5){
							dist = Util.ERP(tra, clustra, idPoints);
						}
					}
					if (min_dist > dist) {
						min_dist = dist; // maintain the one with min distance
						min_id = j;
					}					
					if(yinyang) {// initialize the lower bound
						int groupid=centerGroup.get(j);// get the group numbers
						if(dist < bounds[groupid+2]) {
							bounds[groupid+2] = dist;
						}
					}
				}
			}
			sumDist += min_dist;
			if(yinyang) {// for initialize the bounds
				int groupid=centerGroup.get(min_id);
				bounds[groupid+2] = Double.MAX_VALUE;// set the optimal group as max distance as we do not need to need to consider this group
				bounds[0] = min_dist;// initialize the upper bound
				trajectoryBounds.put(idx, bounds);//the initial bound
			}
			ClusterPath newCluster = new_CENTERS.get(min_id);
			Time1 = System.nanoTime();
			if(distanceModel==0)
				newCluster.updateHistorgramGuava(tra, idx); //update the edge histogram using every new trajectory			
			Time2 = System.nanoTime();						
			newCluster.addTrajectoryToCluster(idx);	// update the new trajectory to this cluster.
			runrecord.addHistorgramTime((Time2-Time1)/1000000000.0);
		}
		System.out.println("the sum distance after assignment is: "+sumDist);
		return new_CENTERS;
	}
	
	/*
	 * given a pivot and computed distance to the centroid, 
	 * judge whether we can assign all the trajectories to this centroid once based on the score
	 */
	public double batchAssign(int pivot, double mindistance, int minid,
			double []bounds, ClusterPath newCluster, int groupnumber, Map<Integer, int[]> clustData) {
		double sum = 0;
		double radius = pivotRadius.get(pivot);
		for(int idx: pivotGroup.get(pivot)) {
			double [] newbounds = new double[groupnumber+2];
			for(int i=0; i<groupnumber+2; i++) {
				newbounds[i] = bounds[i]-radius;//the bound minus with the radius
			}
			trajectoryBounds.put(idx, newbounds);
			int []tra = datamap.get(idx);
			sum += Util.Intersection(tra, clustData.get(minid));// update the sum distance
			numCompute++;
			newCluster.updateHistorgramGuava(tra, idx); //update the edge histogram using every new trajectory		
			newCluster.addTrajectoryToCluster(idx);
		}
		return sum;
	}
	
	//use a queue to conduct the assignment based on pivot-table index, special for EBD
	public ArrayList<ClusterPath> assignPivotTable(int k, ArrayList<ClusterPath> new_CENTERS, boolean yinyang, 
			int groupnumber, Set<Integer> candidateset) {
		Set<Integer> candidateofAllclusters = new HashSet<Integer>();
		Map<Integer, int[]> clustData = new HashMap<Integer, int[]>();
		int minlength = Integer.MAX_VALUE;
		int min_length_id=0;
		double sumDist = 0;
		for (int j = 0; j < k; j++) {
			int[] clustra = CENTERS.get(j).getTrajectoryData();			
			clustData.put(j, clustra);
			if(invertedIndex) {
				Set<Integer> candilist = CENTERS.get(j).creatCandidateList(edgeIndex, datamap);//generate the candidate list
				candidateofAllclusters.addAll(candilist);// merge it to a single list
				if(clustra.length < minlength) {// get the minimum length
					minlength = clustra.length;
					min_length_id = j;
				}
			}
		}
		PriorityQueue<Integer> queueB = new PriorityQueue<Integer>(pivotGroup.keySet());//insert all the pivots
		int counter = 0;//to count the scanned pivot
		while (!queueB.isEmpty()) {// first assignment, check whether we can assign in batch
			int idx = queueB.poll();//pull the candidates
			long Time1 = System.nanoTime();
			int[] tra = datamap.get(idx);//the trajectory data is read
			long Time2 = System.nanoTime();
			runrecord.addIOTime((Time2-Time1)/1000000000.0);			
			double min_dist = Double.MAX_VALUE;
			int min_id = 0;		
			double [] bounds = new double[groupnumber+2];
			Arrays.fill(bounds, Double.MAX_VALUE);// only check for the first rounds
			if(invertedIndex && !candidateofAllclusters.contains(idx)) {//if it is never contained by any list, we can assign it to the cluster with minimum length
				min_dist = Math.max(tra.length, minlength);
				min_id = min_length_id;
				if(yinyang) {// initialize the lower bound
					for (int j = 0; j < k; j++) {
						int length = clustData.get(j).length;
						double dist =  Math.max(tra.length, length);
						int groupNumber=centerGroup.get(j);
						if(j==min_length_id) {
							j=k;
						}
						if(dist<bounds[groupNumber+2]) {
							bounds[groupNumber+2] = dist;
						}
					}
				}
				indexFil+=k;
			}else {
				for (int j = 0; j < k; j++) {
					int centroid = CENTERS.get(j).getTrajectoryID();
					if(idx == centroid) {// the centroid itself.
						min_dist = 0;
						min_id = j;
					}
					Set<Integer> canlist = CENTERS.get(j).getCandidateList();// get the candidate list of each cluster
					double dist = 0;
					int[] clustra = clustData.get(j);
					if(invertedIndex) {
						if(!canlist.contains(idx)) {						// it is not contained
							dist = Math.max(tra.length, clustra.length);
							indexFil++;
						}
						else {			
							dist = (double)Intersection(tra, clustra, tra.length, clustra.length);
							numCompute++;
						}
					}else {
						numCompute++;
						dist = (double)Intersection(tra, clustra, tra.length, clustra.length);
					}
					if (min_dist > dist) {
						min_dist = dist; // maintain the one with min distance
						min_id = j;
					}					
					if(yinyang) {// initialize the lower bound
						int groupid=centerGroup.get(j);// get the group numbers
						if(dist < bounds[groupid+2]) {
							bounds[groupid+2] = dist;
						}
					}
				}
			}
			sumDist += min_dist;
			if(yinyang) {// for initialize the bounds
				int groupid=centerGroup.get(min_id);
				bounds[groupid+2] = Double.MAX_VALUE;// set the optimal group as max distance as we do not need to need to consider this group
				bounds[0] = min_dist;// initialize the upper bound
				trajectoryBounds.put(idx, bounds);//the initial bound
			}
			ClusterPath newCluster = new_CENTERS.get(min_id);
			Time1 = System.nanoTime();
			newCluster.updateHistorgramGuava(tra, idx); //update the edge histogram using every new trajectory			
			Time2 = System.nanoTime();						
			newCluster.addTrajectoryToCluster(idx);	// update the new trajectory to this cluster.
			if(pivotGroup.containsKey(idx) && pivotGroup.get(idx)!=null) {// call the batch assigning if there is an index
				if(min_dist < 2*pivotRadius.get(idx)) {
					queueB.addAll(pivotGroup.get(idx));//cannot pruned
				}else {
					int groupsize = pivotGroup.get(idx).size();
					counter += groupsize;
					sumDist += batchAssign(idx, min_dist, min_id, bounds, newCluster, groupnumber, clustData);
					numFilGlobal+= k*groupsize;
				}
			}
			runrecord.addHistorgramTime((Time2-Time1)/1000000000.0);
		}
		System.out.println("the sum distance after assignment is: "+sumDist+"\n#pruned trajectories by pivot table: "+counter);
		return new_CENTERS;
	}
	
	/*
	 * the data needs to be sorted before the intersection, the edge-based distance (EBD)
	 */
	public int Intersection(int arr1[], int arr2[], int m, int n) {
		int i = 0, j = 0;
		int dist = 0;
		while (i < m && j < n) {
			if (arr1[i] < arr2[j])
				i++;
			else if (arr2[j] < arr1[i])
				j++;
			else
			{
				dist++;
				i++;
				j++;
			}
		}
		return Math.max(m, n)-dist;
	}
	
	/*
	 * single k-path operation
	 */
	public double singleKpath(int k, double overallDis, boolean yinyang, int groupnumber, String folder, 
			Set<Integer> candidateset) {
		if(yinyang)
			printCluterTraID(k, 1, folder);
		PRE_CENS = new ArrayList<ClusterPath>(CENTERS);		//maintain current centers for judging convergence						
		ArrayList<ClusterPath> new_CENTERS = new ArrayList<ClusterPath>(); // it stores the k clusters
		for(int i=0; i<k; i++) {
			if(CENTERS.get(i).getgraphPathExtractionSimpliedRoad()==false) {
				ClusterPath newCluster = new ClusterPath(CENTERS.get(i).getClusterPath().getVIseries(), CENTERS.get(i).getTrajectoryID());
				new_CENTERS.add(newCluster);
			}else {
				ClusterPath newCluster = new ClusterPath(CENTERS.get(i).getClusterPath().getVIseries(), CENTERS.get(i).getTrajectoryID(),
						CENTERS.get(i).getgraphPathExtractionSimpliedRoad(), CENTERS.get(i).getnonImportantRoad());
				new_CENTERS.add(newCluster);
			}
		}
		long startTime1 = System.nanoTime();
		if(pivotPruning) {//
			new_CENTERS = assignPivotTable(k, new_CENTERS, yinyang, groupnumber, candidateset);
		}else
			new_CENTERS = assignRebuildInvertedindex(k, new_CENTERS, yinyang, groupnumber, candidateset);	//update the CENTERS
		long endtime = System.nanoTime();
		runrecord.addAssignmentTime((endtime-startTime1)/1000000000.0);
		System.out.print("assign time cost: ");
		System.out.printf("%.3f", (endtime-startTime1)/1000000000.0);
		System.out.println("s");
		long startTime = System.nanoTime();
		for(int i=0; i<k; i++) {// generate the new centroid for each cluster
			if(distanceModel == 0 && refinementSign) {
				if(graphPathExtraction)// using the graph information
					new_CENTERS.get(i).extractNewPathCPEP(forwardGraph, backwardGraph, i); // test the optimal		
				else {
					new_CENTERS.get(i).extractNewPathGuava(datamap, runrecord, traLength, trajectoryHistogram); 
				}
			}else {//using the distance matrix based k-medoids 
				new_CENTERS.get(i).extractNewPathStupid(datamap, runrecord, distanceModel, idPoints);
			}
			overallDis += new_CENTERS.get(i).getSumDistance();
		}
		endtime = System.nanoTime();
		runrecord.addRefinementTime((endtime-startTime)/1000000000.0);
		CENTERS = new ArrayList<ClusterPath>(new_CENTERS);
		if(yinyang)
			System.out.println("iteration 1, the sum distance is "+overallDis+", time cost: "+(endtime-startTime1)/1000000000.0+"s\n");		
		return overallDis;
	}
	
	/*
	 * conduct the clustering, we are using the Lloyd's algorithm
	 */
	public int kPath(int k, String folder) {
		int t = 0;
		for(; t < TRY_TIMES; t++){
			printCluterTraID(k, t, folder);
			double overallDis = 0;
			overallDis = singleKpath(k, overallDis, false, 0, folder, datamap.keySet());
			System.out.println("iteration "+(t+1)+", the sum distance is "+overallDis);
			if(timeToEnd()) {
				System.out.println("\nIteration stops now");
				runrecord.setIterationtimes(t+1);
				break;//convergence
			}
		}
		return t;
	}
	
	public void Initialization() throws IOException {
		CENTERS = new ArrayList<ClusterPath>();	
		interMinimumCentoridDis = new double[k];
		innerCentoridDis = new double[k][];
		innerCentoridDisGroup = new double[k][];
		edgeIndex = new HashMap<>();
		edgeHistogram = new HashMap<>();
		trajectoryHistogram = new HashMap<>();
		edgeInfo = new HashMap<>();	
		datamap = new HashMap<>();
		traLength = new HashMap<>();
		Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat("DDHHmmss");
        folder = sdf.format(cal.getTime());
	}
	

	/*
	 * this methond will finish the baseline and pivot-based index.
	 */
	public void staticKpath(String[] args, boolean indexbased) throws IOException {		
		mtreebuild = indexbased;
		if(mtreebuild) {
			second = false;
			graphPathExtraction = true;
		}else {
			second = true;
			graphPathExtraction = false;
		}		
		int iterations=0;
		Initialization();
		loadData(datafile, trajectoryNumber, edgefile);	// load the data and create index
		if(second==true)
			createTrajectoryHistogram(datamap, trajectoryNumber);  // build inverted index if there is not index
        String new_folder = mapv_path+folder;
        boolean creat_folder = new File(new_folder).mkdirs();
        if(!creat_folder)
        	return;
       
        if(mtreebuild) {
        	mindex.buildGraph(graphfile, forwardGraph, backwardGraph);
        	mindex.runkpath(folder,args);
        }else {
        	initializeClustersIncrease(k, 100);
            //	initializeClustersRandom(k); 		// initialize the kmeans in the first iteration
        	// initializeClustersHighFrequency(k, 1000);
        	iterations = kPath(k, new_folder);        
        }
	//	SIGMOD07.SIGMOD07FrequentEdgesMapv(edgeHistogram, edgeInfo, frequencyThreshold, mapv_path_traclu_sigmod07+Integer.toString(frequencyThreshold));		
		if(graphPathExtraction == false){
			runrecord.printLog();
			System.out.println("\nThe final centroid trajectory ids are:");
			int []result = printCluterTraID(k, iterations+1, folder);
		}
	}
	
	public static void main(String[] args) throws IOException, SQLException, InterruptedException {
	//	Process run = new Process(args);
	//	DataReading.convertToEdges("E:\\dataset\\nantong\\nantong_simple_date.csv", "E:\\dataset\\nantong\\new_nantong_simple_date.csv", "E:\\dataset\\nantong\\edge_mapping", "E:\\dataset\\nantong\\car_mapping");
	//	DataReading.convertToEdges("E:\\dataset\\nantong\\nantong_7_1.txt", "E:\\dataset\\nantong\\test.csv", "E:\\dataset\\nantong\\edge_mapping1", "E:\\dataset\\nantong\\car_mapping1");
	//	long s = System.nanoTime();
	//	StreamKpath.readDataFromFile("E:\\dataset\\nantong\\new_nantong_simple_date.csv");
	//	DataReading.convertStandardFormat("E:\\dataset\\nantong\\nantong.csv", "E:\\dataset\\nantong\\nantong_traffic_record");		
	//	System.out.println((System.nanoTime()-s)/1000000000.0);
	//	testStreamkPath();
	//	run.staticKpath(args);
	}
}
