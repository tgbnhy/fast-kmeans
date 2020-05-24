package au.edu.rmit.trajectory.clustering.kpaths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

/*
 * this store the information of each cluster
 */
public class ClusterPath {
	protected Set<Integer> clusterTrajectory; // it stores all the idxs of trajectories
	protected Set<Integer> clusterPivot; // it stores all the pivot trajectories assigned to this clusters
	protected VIseries centroid; // store the centorid path
	int iteration = 5000; // the maximum iteration times in the CPEP, we need to make it changeable.
	
	// build histogram for edges and length using hashmap and Guava
	protected Multiset<Integer> edgeOcc; 
	protected Multiset<Integer> lengthOcc;
	protected boolean centerChanged;
	int minlength, maxlength;
	int []lengthAccu;//for the score computation
	Map<String, Double> checkedList;//store the start, end edges, and the score
    PriorityQueue<Path> queue;
    ArrayList<Integer> sortedFrequency;
    int sumEdgeOcc;
    int pathMinlength;
	
	protected int idx;// it stores the idx of the trajectory which is the centroid.
	protected double sumdistance=0; // the sum distance in this cluster.
	protected Set<Integer> candidateList = new HashSet<>();//built from the inverted index
	protected int []finalPath;
	
	boolean graphPathExtractionSimpliedRoad;
	ArrayList<Integer> nonImportantRoad;
	
	public ClusterPath(int[] cl, int idx1) {		
		clusterTrajectory = new HashSet<Integer>();
		clusterPivot = new HashSet<>();
		centroid = new VIseries();
		centroid.setVIseries(cl);
		if(cl!=null)
			centroid.length = cl.length;
		else
			centroid.length = 0;
		centroid.idx = idx1;
		this.idx = idx1;
		this.finalPath = cl;
		this.edgeOcc = HashMultiset.create();
		this.lengthOcc = HashMultiset.create();
		graphPathExtractionSimpliedRoad = false;
		sumEdgeOcc = 0;
		centerChanged = true;
		finalPath = new int[cl.length];
		int i=0;
		for(int edge:cl) {
			finalPath[i++] = edge;
		}
	}
	
	public int getEdgeHistogramSize() {
		
		return edgeOcc.size();
	}
	
	public ClusterPath(Set<Integer> cl, int idx1) {		
		clusterTrajectory = cl;
		this.idx = idx1;
	}
	
	public ClusterPath(int[] cl, int idx1, boolean graphPathExtractionSimpliedRoad1, ArrayList<Integer> nonImportantRoad) {		
		clusterTrajectory = new HashSet<Integer>();
		centroid = new VIseries();
		centroid.setVIseries(cl);
		if(cl!=null)
			centroid.length = cl.length;
		else
			centroid.length = 0;
		centroid.idx = idx1;
		this.idx = idx1;
		this.finalPath = cl;
		this.edgeOcc = HashMultiset.create();
		this.lengthOcc = HashMultiset.create();
		sumEdgeOcc = 0;
		centerChanged = true;
		graphPathExtractionSimpliedRoad = true;
		this.nonImportantRoad = nonImportantRoad;
	}
	
	public void setIteration(int iteration) {
		this.iteration = iteration;
	}

	public boolean getgraphPathExtractionSimpliedRoad() {
		return graphPathExtractionSimpliedRoad;
	}
	
	public ArrayList<Integer> getnonImportantRoad() {
		return nonImportantRoad;
	}
	/*
	 * generate the list id of trajectory that share edge with cluster
	 */
	public Set<Integer> creatCandidateList(Map<Integer, Set<Integer>> edgeIndex, Map<Integer, int[]> datamap) {
		if (centerChanged) {//create a new list if changed, otherwise return previous to avoid rebuild
			int[] trajectory = datamap.get(idx);// if the centroid is selected from
			if ( finalPath != null ) {
				trajectory = finalPath;
			}
			candidateList = new HashSet<Integer>();
			for (int edgeid : trajectory) {
				Set<Integer> traidlist = edgeIndex.get(edgeid);
				Collections.addAll(candidateList, traidlist.toArray(new Integer[0]));
			}
		}
		return candidateList;
	}
	
	/*
	 * generate the list id of trajectory that share edge with cluster, or we can ignore this
	 */
	public Set<Integer> creatCandidateListNoDatamap(Map<Integer, Set<Integer>> edgeIndex, int[] trajectory) {
		if (centerChanged) {//create a new list if changed, otherwise return previous to avoid rebuild
			candidateList = new HashSet<Integer>();
			for (int edgeid : trajectory) {
				Set<Integer> traidlist = edgeIndex.get(edgeid);
				Collections.addAll(candidateList, traidlist.toArray(new Integer[0]));
			}			
		}
		return candidateList;
	}
	
	public Set<Integer> getCandidateList() {
		return candidateList;
	}
	
	/*
	 * add the trajectory into the clusters
	 */
	void mergeTrajectoryToCluster(ArrayList<Integer> index){
		clusterTrajectory.addAll(index);
	}
	
	/*
	 * remove the trajectory into the clusters
	 */
	void removeTrajectoryToCluster(ArrayList<Integer> index){
		clusterTrajectory.removeAll(index);
	}
	
	/*
	 * add the trajectory into the clusters
	 */
	void addTrajectoryToCluster(int index){
		clusterTrajectory.add(index);
	}
	
	/*
	 * remove the trajectory into the clusters
	 */
	void removeTrajectoryToCluster(int value){
		if(clusterTrajectory.contains(value)) {
			clusterTrajectory.remove(new Integer(value));
		}
	}
		
	public void updateHistorgramGuava(int[] tra, int idx) {
		lengthOcc.add(tra.length);
		for(int edge: tra) {
			edgeOcc.add(edge);
		}
	}
	
	void removeHistorgramGuava(int[] tra, int idx) {
		lengthOcc.remove(tra.length, 1);
		for(int edge: tra) {
			edgeOcc.remove(edge, 1);
		}
	}
	
	public void updateHistorgramGuava(Multiset<Integer> edgeH, Multiset<Integer> lengthH) {
		edgeOcc.addAll(edgeH);
		lengthOcc.addAll(lengthH);
	}
	
	public void removeHistorgramGuava(Multiset<Integer> edgeH, Multiset<Integer> lengthH) {
		edgeOcc.removeAll(edgeH);
		lengthOcc.removeAll(lengthH);
	}
	
	/*
	 * return the sum distance in the cluster
	 */
	public double getSumDistance() {
		return sumdistance;
	}
	
	/*
	 * return the clusters
	 */
	VIseries getClusterPath() {
		return centroid;
	}
	
	/*
	 * return the trajectory array of this cluster
	 */
	public Set<Integer> getClusterTrajectories() {
		return clusterTrajectory;
	}
	
	/*
	 * return the trajectory array of this cluster
	 */
	public Set<Integer> printClusterTrajectories() {
		System.out.print(idx+":");
		for(int traid:clusterTrajectory)
			System.out.print(traid+",");
		System.out.println();
		return clusterTrajectory;
	}
	
	public int getTrajectoryID() {
		return idx;
	}
	
	public int[] getTrajectoryData() {
		return finalPath;
	}
	
	public int[] getTrajectoryDataUnsorted() {
		return centroid.getVIseries();
	}
	
	public boolean getCenterChanged() {
		return centerChanged;
	}
	
	//for bound computing
	int getFirstNSum(ArrayList<Integer> sortedFrequency, int number) {
		int sum = 0;
		for(int i=0; i<number && i< sortedFrequency.size(); i++) {
			sum+=sortedFrequency.get(i);
		}
		return sum;
	}
	
	/*
	 * there is no trajectory in the center
	 */
	public void accumulateLenghtOcc(){				
		ArrayList<Integer> sortedTralen = new ArrayList<Integer>(lengthOcc.elementSet());
		Collections.sort(sortedTralen);		//sorted length increasingly
		if (!lengthOcc.isEmpty()) {
			minlength = sortedTralen.get(0);
			maxlength = sortedTralen.get(sortedTralen.size() - 1);
			lengthAccu = new int[maxlength - minlength + 1];
			sumEdgeOcc = 0;
			sortedFrequency = new ArrayList<Integer>();
			for (int edge : edgeOcc.elementSet()) {// reverse way
				sortedFrequency.add(edgeOcc.count(edge));
			}
			Collections.sort(sortedFrequency, Collections.reverseOrder());
			for (int length = minlength; length <= maxlength; length++) {// this is not right
				int occ = 0;
				for (int i = length; i > minlength; i--) {
					if (lengthOcc.contains(i))
						occ += lengthOcc.count(i)*(i-minlength);// add the extra length
				}
				lengthAccu[length - minlength] = occ;
				sumEdgeOcc += lengthOcc.count(length)*length;
			}
		}
	}
	

	/*
	 * We use Multiset in the google library to update the histogram in a faster way by scanning each trajectory
	 */
	double extractNewPathGuava(Map<Integer, int[]> datamap, RunLog runrecord, Map<Integer, Integer> traLengthmap, 
			Map<Integer, Integer> trajectoryHistogram) {
		int min = Integer.MAX_VALUE;
		int traid = 0 ;
		centerChanged = true;
		accumulateLenghtOcc();
		for(int traidx: clusterTrajectory) {		//read each trajectory in this cluster, it is slow as we need to read every trajectory, we can construct a weighted graph and choose the route			
			int traLength = traLengthmap.get(traidx);
		//	int bound = Math.min(trajectoryHistogram.get(traidx), getFirstNSum(sortedFrequency, traLength));
			int bound = getFirstNSum(sortedFrequency, traLength);
			int sumLength = sumEdgeOcc;// the sum of lengths
			sumLength += lengthAccu[traLength - minlength];			
			if((sumLength - bound) >= min) {
				continue;// skip reading the trajectory data
			}
			int[] tra = datamap.get(traidx);
			int sumFrequency = 0;
			int previous = -1;
			for(int edge: tra) {		//compute the vertex frequency of each trajectory
				if(previous == edge) {
					continue;// avoid adding many times
				}
				sumFrequency += edgeOcc.count(edge);//
				previous = edge;
			}
			int sumDis = sumLength - sumFrequency;
			if(min > sumDis) {//choose the one which can minimize the sum of distance
				min = sumDis;
				traid = traidx;
			}
		}
		if(idx==traid)
			centerChanged = false;//center does not change;
		else
			centerChanged = true;
		sumdistance = min;
		int [] a = datamap.get(traid);
		double drift = 0;
		if(centerChanged && a!=null)
			drift = Intersection(finalPath, a, finalPath.length, a.length);
		centroid.setVIseries(a);
		centroid.idx = traid;
		idx = traid;
		finalPath = a;
		return drift;
	}
	
	/*
	 * this method will build a matrix first, then choose the one that has the minimum sum.
	 */
	void extractNewPathStupid(Map<Integer, int[]> datamap, RunLog runrecord, int disModel, Map<Integer, double[]> idPoints) {
		ArrayList<Integer> arrayList = new ArrayList<>(clusterTrajectory);
		if(clusterTrajectory.isEmpty())
			return;
		int size = clusterTrajectory.size();
		double [][]matrix = new double[size][size];
		centerChanged = true;
		for(int i=0; i<size; i++) {
			for(int j=0; j<size; j++) {
				if(i==j) {
					matrix[i][j] = 0;
					continue;
				}
				if(matrix[i][j]==0) {//compute the matrix
					if(disModel==0)
						matrix[i][j] = Util.Intersection(datamap.get(arrayList.get(i)), datamap.get(arrayList.get(j)));
					else if(disModel==1)
						matrix[i][j] = Util.EDR(datamap.get(arrayList.get(i)), datamap.get(arrayList.get(j)));
					else if(disModel==2)
						matrix[i][j] = Util.DTW(datamap.get(arrayList.get(i)), datamap.get(arrayList.get(j)), idPoints);
					else if(disModel==3)
						matrix[i][j] = Util.Frechet(datamap.get(arrayList.get(i)), datamap.get(arrayList.get(j)), idPoints);
					else if(disModel==4)
						matrix[i][j] = Util.Hausdorff(datamap.get(arrayList.get(i)), datamap.get(arrayList.get(j)), idPoints);
					else if(disModel==5)
						matrix[i][j] = Util.ERP(datamap.get(arrayList.get(i)), datamap.get(arrayList.get(j)), idPoints);
				}else if(matrix[j][i] == 0){
					matrix[j][i] = matrix[i][j];
				}
			}
		}
		double min = Double.MAX_VALUE;
		int traid =0;
		for(int i=0; i<size; i++) {
			double sum = 0;
			for(int j=0; j<size; j++) {
				sum += matrix[i][j];
			}
			if(sum<min) {
				min = sum;
				traid = arrayList.get(i);
			}
		}
		if(idx==traid)
			centerChanged = false;//center does not change;
		else
			centerChanged = true;
		sumdistance = min;
		int [] a = datamap.get(traid);
		centroid.setVIseries(a);
		centroid.idx = traid;		
		idx = traid;
		finalPath = a;
	}
	
	/* A greedy algorithm when the optimal result does not change any more and better than last result, 
	 * we will stop exploring when the length is close to the maximum length it can be
	 */
	public double extractNewPathCPEP(
			HashMap<Integer, ArrayList<Integer>> forwardGraph, 
			HashMap<Integer, ArrayList<Integer>> backwardGraph, int clusterid) {
		checkedList = new HashMap<>();
		queue = new PriorityQueue<Path>();
		ArrayList<Integer> optimal = new ArrayList<>();
		
		accumulateLenghtOcc();
		double optimalScore = Double.MAX_VALUE;
/*		if(finalPath!=null) {
			centerChanged = false;
			optimalScore = computeScore(finalPath);//initialized to min score using last centroid
			for(int a:centroid.getVIseries())//using the unsorted
				optimal.add(a);
			double lowerbound = estimateLowerbound(optimalScore, optimal.size());
			if(optimalScore > lowerbound) {
				Path aPath = new Path(optimal, optimalScore, lowerbound);
				queue.add(aPath);
			}
		}*/
		for(int edge: edgeOcc.elementSet()) {//initialize the edges as paths
			if(graphPathExtractionSimpliedRoad && nonImportantRoad.contains(edge)) {
				continue;//remove the road will never shown in the centroids
			}
			ArrayList<Integer> arrayList = new ArrayList<>();
			arrayList.add(edge);
			double score = computeScore(arrayList);
			double lowerbound = estimateLowerbound(score, 1);
			if(optimalScore > lowerbound) {
				Path aPath = new Path(arrayList, score, lowerbound);
				queue.add(aPath);
			}
		}
		
		int cou= 0;
	//	System.out.println("\n");//termination times
		int terminationTime = 0;
		while(!queue.isEmpty()) {
			cou++;
			Path candidate = queue.poll();
			double lowerbound = candidate.getLowerbound();// compute the score	
			double score = candidate.getScore();
		//	System.out.println(cou+ "\t" +optimalScore);
			if(lowerbound > optimalScore || cou>=iteration) {//termination as all possible has been checked
				break;
			}			
			ArrayList<Integer> Can = candidate.getPath();
			int start = Can.get(0);
			int end = Can.get(Can.size() - 1);// the last edge			
			if(backwardGraph.containsKey(start)) {
			ArrayList<Integer> backAppend = backwardGraph.get(start);			
			if(backAppend!=null)
			for (int ids : backAppend) {
				if(graphPathExtractionSimpliedRoad && nonImportantRoad.contains(ids)) {
					continue;//remove the road will never shown in the centroids
				}
				if (edgeOcc.contains(ids) && !Can.contains(ids)) {
					ArrayList<Integer> newCan = new ArrayList<>(Can);
					newCan.add(0, ids);	//insert to the beginning					
					if(forwardGraph.containsKey(end)) {
						ArrayList<Integer> startAppend1 = forwardGraph.get(end);// no circle
						if(startAppend1!=null && startAppend1.contains(ids)) {
							continue;
						}
					}		
					double newscore = computeScoreWithPrevious(ids, score, newCan.size());
					if (newscore < optimalScore) {
						terminationTime = cou;
						centerChanged = true; // the center has changed
						optimalScore = newscore;
						optimal = newCan;
					}
					lowerbound = estimateLowerbound(newscore, newCan.size());
					String signature = signaturePath(ids, end);
					if (neverChecked(signature, lowerbound)) {
						continue;
					}
					if (lowerbound < optimalScore) {
						Path newpath = new Path(newCan, newscore, lowerbound);
						queue.add(newpath);
					}
				}
			}
			}
			
			if(forwardGraph.containsKey(end)) {
			ArrayList<Integer> startAppend = forwardGraph.get(end);
			
			if(startAppend!=null)
			for (int ids : startAppend) {
				if(graphPathExtractionSimpliedRoad && nonImportantRoad.contains(ids)) {
					continue;//remove the road will never shown in the centroids
				}
				if (edgeOcc.contains(ids) && !Can.contains(ids)) {//connected and no repetitive edges
					ArrayList<Integer> newCan = new ArrayList<>(Can);
					newCan.add(ids);				
					if(forwardGraph.containsKey(ids)) {
						ArrayList<Integer> startAppend1 = forwardGraph.get(ids);
						if(startAppend1!=null && startAppend1.contains(start)) {
							continue;
						}
					}				
					double newscore = computeScoreWithPrevious(ids, score, newCan.size());
					if (newscore < optimalScore) {
						terminationTime = cou;
						centerChanged = true;
						optimalScore = newscore;
						optimal = newCan;
					}
					lowerbound = estimateLowerbound(newscore, newCan.size());
					String signature = signaturePath(start, ids);
					if (neverChecked(signature, lowerbound))	//repetitive path
						continue;
					if (lowerbound < optimalScore) {//if this path is promising
						Path newpath = new Path(newCan, newscore, lowerbound);
						queue.add(newpath);
					}
				}
			}
			}
		}
	//	System.out.println("terminate at\t"+terminationTime);
		sumdistance = optimalScore;
		ArrayList<Integer> graphPath = new ArrayList<>(optimal);
		int[] centoridDataGraph = graphPath.stream().mapToInt(i -> i).toArray();
		centroid.setVIseries(centoridDataGraph);
	//	System.out.println(graphPath.toString());
	//	System.out.println(optimalScore+" "+optimal.size());
		Collections.sort(optimal);
		int[] centoridData = optimal.stream().mapToInt(i -> i).toArray();
		double drift=0;
		if(finalPath!=null && centerChanged)
			drift = Intersection(finalPath, centoridData, finalPath.length, centoridData.length);

		finalPath = centoridData;
		return drift;
	}
	
	/*
	 * the data needs to be sorted before the intersection
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
	 * check whether it exists in the checked list
	 */
	public boolean neverChecked(String aString, double lowerbound) {
		if(checkedList.containsKey(aString)) {
			double previous = checkedList.get(aString);
			if(lowerbound >= previous)
				return false;
			else {
				checkedList.put(aString, lowerbound);// add to the list
				return true;
			}
		}else {
			checkedList.put(aString, lowerbound);// add to the list
			return true;
		}
	}
	
	/*
	 * compute the objective score based on the edge histogram and length histogram,
	 */
	public double computeScore(int[] path) {
		double weight = sumEdgeOcc;
		for(int i=0; i<path.length; i++) {
			weight -= edgeOcc.count(path[i]);
		}
		if(path.length>minlength && path.length<=maxlength) {
			weight += lengthAccu[path.length-minlength];
		}else if(path.length>maxlength) {
			weight -= sumEdgeOcc;
			weight += path.length*clusterTrajectory.size();
		}
		return weight;
	}
	
	/*
	 * compute the objective score based on the edge histogram and length histogram
	 */
	public double computeScore(ArrayList<Integer> path) {
		int []patharray =  path.stream().mapToInt(i -> i).toArray();
		return computeScore(patharray);
	}
	
	/*
	 * compute the objective score based on the edge histogram and length histogram
	 */
	public double computeScoreWithPrevious(int newedge, double weight, int length) {
		weight -= edgeOcc.count(newedge);
		if(length > minlength && length<=maxlength)
			weight += lengthAccu[length - minlength] - lengthAccu[length - minlength-1];
		return weight;
	}
	
	/*
	 * convert to a short unique string using the start and end edge, and the length, this can be used as a dominance table, 
	 */
	public String signaturePath(int a, int b) {
		String key = a+"_"+b;
		return key;
	}
	
	/*
	 * add the rest heavy edges to compute the lower bound of the score
	 */
	public double estimateLowerbound(double score, int length) {
		double lowerbound = score;
		int i = length;
		while (i <= minlength && (i - length) < sortedFrequency.size()) {
			lowerbound -= sortedFrequency.get(i - length);
			i++;
		}
		while (i <= maxlength && ((i - length) < sortedFrequency.size()) 
				&& (lengthAccu[i - minlength] -lengthAccu[i - minlength-1]) < sortedFrequency.get(i - length)) {
			lowerbound += lengthAccu[i - minlength] - lengthAccu[i - minlength-1] - sortedFrequency.get(i - length);
			i++;
		}
		return lowerbound;
	}
	
	/*
	 * the maximum length of the final output path
	 */
	public int estimateMax() {
		int i=minlength+1;
		while((lengthAccu[i - minlength] -lengthAccu[i - minlength-1]) < sortedFrequency.get(i - minlength)) {
			i++;
		}
		return i;
	}
}
