package au.edu.rmit.trajectory.clustering.kpaths;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Comparator;
import java.util.HashMap;

import au.edu.rmit.trajectory.clustering.visualization.mapv;
import java_cup.internal_error;

import java.util.stream.Collectors;

public class SIGMOD07 {

	public SIGMOD07() {
		// TODO Auto-generated constructor stub
	}

	/*
	 * this will produce the edges that have a higher frequency than the threshold, SIGMOD07 Trajectory clustering: a partition-and-group framework
	 * into the file-output,
	 */
	public static void SIGMOD07FrequentEdgesMapv(Map<Integer, Integer> edgeHistogram, Map<Integer, String> edgeInfo,
			int frequencyThreshold, String output) {
		ArrayList<Integer> highFreEdge = new ArrayList<Integer>();
		for (int edgeid : edgeHistogram.keySet()) {
			if (edgeHistogram.get(edgeid) > frequencyThreshold) {
				highFreEdge.add(edgeid);
			}
		}
		mapv.generateHighEdges(edgeInfo, highFreEdge, output);
	}
	
	public static void topkfrequentedegs(Map<Integer, Integer> edgeHistogram, Map<Integer, String> edgeInfo,
			int k, String output) {
		ArrayList<Integer> highFreEdge = new ArrayList<Integer>();
		Map<Integer, Integer> sorted = edgeHistogram
		        .entrySet()
		        .stream()
		        .sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
		        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2,
		                LinkedHashMap::new));
		int i =0;
		for (int edgeid : sorted.keySet()) {
			System.out.println(sorted.get(edgeid));
			if (i++<k) {
				highFreEdge.add(edgeid);
			}else {
				break;
			}
		}
		mapv.generateHighEdges(edgeInfo, highFreEdge, output);
	}
	
	/*
	 * this functions selects the trajectories with top k sum frequency.
	 */
	public static void topkTrajectoryWithHighestFrequencySum(Map<Integer, Integer> edgeHistogram, Map<Integer, String> edgeInfo,
			int k, String output, Map<Integer, int[]> datamapunsorted) {
		ArrayList<int[]> highFreEdge = new ArrayList<int[]>();
		Map<Integer, Integer> trajectoryFrequency = new HashMap<Integer, Integer>();
		
		
		
		int i =0;
		for(int trajectoryID: datamapunsorted.keySet()) {
			int[] sequnce = datamapunsorted.get(trajectoryID);
			int frequencySum = 0;
			for(int edgeID: sequnce) {
				frequencySum += edgeHistogram.get(edgeID);
			}
			int average = frequencySum/sequnce.length;
		//	trajectoryFrequency.put(trajectoryID, frequencySum);
			trajectoryFrequency.put(trajectoryID, average);
		}
		
		Map<Integer, Integer> sorted = trajectoryFrequency
		        .entrySet()
		        .stream()
		        .sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
		        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2,
		                LinkedHashMap::new));
		
		
		// compute the sum frequency of each trajectory, and choose the highest as
		for (int trajID : sorted.keySet()) {
		//	System.out.println(sorted.get(trajID));
			if (i++<k) {
				highFreEdge.add(datamapunsorted.get(trajID));//add the trajectory
			}else {
				break;
			}
		}
		mapv.generateClusterPath1(datamapunsorted, edgeInfo, highFreEdge, output+".10frequentTrajectoryAverage");
	}
}
