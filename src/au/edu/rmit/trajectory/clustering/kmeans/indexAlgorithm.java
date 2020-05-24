package au.edu.rmit.trajectory.clustering.kmeans;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import org.netlib.util.doubleW;

import au.edu.rmit.trajectory.clustering.kpaths.Util;
import edu.wlu.cs.levy.cg.KDTree;
import edu.wlu.cs.levy.cg.KeyDuplicateException;
import edu.wlu.cs.levy.cg.KeySizeException;
import es.saulvargas.balltrees.BallTreeMatrix;
import java_cup.internal_error;
import scala.annotation.implicitNotFound;
import skyline0623.balltree.BallTree;
import skyline0623.balltree.Hypersphere;
import skyline0623.balltree.Point;
import skyline0623.balltree.Process;
import covertree.*;
public class indexAlgorithm<E> {

	int distanceCompute = 0;
	int NodeAccess = 0;
	int dataAccess = 0;
	
	public indexAlgorithm() {
		// TODO Auto-generated constructor stub
	}
	
	public void buildKDtree(int dims, double[][] itemMatrix) throws KeySizeException, KeyDuplicateException {
		KDTree<Long> kt = new KDTree<Long>(dims);
		long idx = 1;
		for(double[] point: itemMatrix) {
			kt.insert(point, idx++);//fast construction: point, value
		//	System.out.println("inset kd"+idx);
		}		
		indexNode rootKmeans = new indexNode(dims);
	//	kt.traverseConvert(rootKmeans, kt.getroot(), dims);	// traversing the index is hard for kd-tree
	}
	
	/*
	 * get the count of the tree.
	 */
	public double[] updateSum(indexNode root, int dimension) {
		if(root.isLeaf()) {	
			return root.getSum();
		}
		else {
			Set<indexNode> listnode = root.getNodelist();
		//	System.out.println(listnode.size());
			double []sum = new double[dimension];
			for(indexNode aIndexNode: listnode) {
				double []sumd = updateSum(aIndexNode, dimension);
				for(int i=0; i<dimension; i++)
					sum[i] += sumd[i];
			}
			root.setSum(sum);
			return sum;
		}
	}
	
	public indexNode buildMtree(double[][] itemMatrix, int dimension, int capacity) {// too slow
		System.out.println("Building M-tree...");
		long startTime1 = System.nanoTime();
		PointMtree mindex = new PointMtree(capacity);		//capacity
		int idx = 1;
		for(double[] point: itemMatrix) {
			mindex.buildMtree(point, idx++);//create the M-tree
		//	System.out.println("inset M-tree"+idx);
		}
		indexNode rootKmeans = new indexNode(dimension);
		mindex.traverseConvert(rootKmeans, mindex.getRoot(), dimension);		// traversing the index
		updateSum(rootKmeans, dimension);
		long endtime = System.nanoTime();
		System.out.println((endtime-startTime1)/1000000000.0);
	//	System.out.println("the count of M-tree is " + rootKmeans.getTotalCoveredPoints());
		return rootKmeans;
	}
	
	
	public indexNode buildCoverTree(int dims, double[][] datamapEuc) {
		//build the tree based on 
		CoverTree<Integer> cTree = new CoverTree<Integer>();
		int idx=0;
		for(double[] point:  datamapEuc) {
			cTree.insert(idx, point);//fast construction: point, value
		//	System.out.println("inset cover tree"+idx);
			idx++;
		}
		indexNode rootKmeans = new indexNode(dims);
		cTree.traverseConvert(rootKmeans, dims);
		updateSum(rootKmeans, dims);
		return rootKmeans;
		
	}
	
	public indexNode buildBalltree(Map<Integer, double[]> datamapEuc, int dimension, int capacity) {// too slow	
	//	System.out.println("Building Ball-tree...");
		long startTime1 = System.nanoTime();	
		for(int idx: datamapEuc.keySet()) {
			double[] point = datamapEuc.get(idx);
			Process.DIMENSION = dimension;
			Process.INSTANCES.add(new Point(point));
			Process.MAX_INSTANCE_NUM_NOT_SPLIT = capacity;
		//	System.out.println("inset Ball-tree"+idx);
		}
		Hypersphere BALL_TREE = BallTree.buildAnInstance(null);		
		indexNode rootKmeans = new indexNode(dimension);
		BALL_TREE.traverseConvert(rootKmeans, dimension);
	//	computeFartherToChild(rootKmeans);
		long endtime = System.nanoTime();
		System.out.print((endtime-startTime1)/1000000000.0+",");
		System.out.println("the count of Ball-tree is " + 
				rootKmeans.getTotalCoveredPoints()+", the radius is "+rootKmeans.getRadius());
	//	System.out.println("the count of Ball-tree is " + rootKmeans.getTotalCoveredPoints());
		return rootKmeans;
	}
	
	public indexNode buildBalltree2(double[][] itemMatrix, int dimension, int capacity) {// too slow	
		System.out.println("Building Ball-tree using Matrix...");
		long startTime1 = System.nanoTime();
		int deepth = (int) (Math.log(itemMatrix.length)/Math.log(2));//the deepth is computed based on binary tree
		indexNode rootKmeans = BallTreeMatrix.create(itemMatrix, capacity, deepth);//we should not set the deepth too deep
		updateSum(rootKmeans, dimension);
		long endtime = System.nanoTime();
		System.out.println("index time cost: "+(endtime-startTime1)/1000000000.0);
		System.out.println("the count of Ball-tree using Matrix is " + 
			rootKmeans.getTotalCoveredPoints()+", the radius is "+rootKmeans.getRadius());
		return rootKmeans;
	}
	
	
	/*
	 * get all the points under a node
	 */
	public ArrayList<Integer> getPointsIDlist(indexNode root) {
		ArrayList<Integer> result = new ArrayList<>();
		if(root == null)
			return null;
		if(root.isLeaf()) {	
			result.addAll(root.getpointIdList());
		}
		else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				result.addAll(getPointsIDlist(aIndexNode));
			}
		}
		return result;
	}
	
	public void setdistanceCompute(int distanceCompute) {
		this.distanceCompute = distanceCompute;
	}
	
	public int getdistanceCompute() {
		return distanceCompute;
	}
	
	public void setNodeAccess(int nodeAccess) {
		this.NodeAccess = nodeAccess;
	}
	
	public int getNodeAccess() {
		return NodeAccess;
	}
	
	
	public void setdataAccess(int dataAccess) {
		this.dataAccess = dataAccess;
	}
	
	public int getdataAccess() {
		return dataAccess;
	}
	
	/*
	 * search all the points within a distance radius, known as similarity search, using M-tree may be better.
	 */
	public ArrayList<Integer> SimilaritySearchBall(double radius, double point[], indexNode root, int dimension, 
			double[][] itemMatrix) { 
		ArrayList<Integer> result = new ArrayList<>();
		if (root.isLeaf()) {
			for (int id : root.getpointIdList()) {
				dataAccess++;
				double distance = Util.EuclideanDis(itemMatrix[id - 1], point, dimension);
				distanceCompute++;
				if (distance < radius)
					result.add(id);
			}
		}else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				NodeAccess++;
				double distance = Util.EuclideanDis(aIndexNode.getPivot(), point, dimension);
				distanceCompute++;
				if(radius - distance >= aIndexNode.getRadius()) {
				//	System.out.println("aaaaaaaaaaaaaaaaaa");
					result.addAll(getPointsIDlist(aIndexNode));
				}else {// cannot be pruned
					result.addAll(SimilaritySearchBall(radius, point, aIndexNode, dimension, itemMatrix));
				}
			}
		}
		return result;
	}
	
	/*
	 * search the nearest neighbor
	 */
	public void NearestNeighborSearchBall(double point[], indexNode root, int dimension, 
			double[][] itemMatrix, double []minDistnearestID) { 
		if (root.isLeaf()) {
			for (int id : root.getpointIdList()) {
				double distance = Util.EuclideanDis(itemMatrix[id - 1], point, dimension);
				distanceCompute++;
				if (distance < minDistnearestID[0]) {
					minDistnearestID[0] = distance;
					minDistnearestID[1] = (double)id;
				}
			}
		}else {
			Set<indexNode> listnode = root.getNodelist();
			for(indexNode aIndexNode: listnode) {
				double distance = Util.EuclideanDis(aIndexNode.getPivot(), point, dimension);
				distanceCompute++;
				if(minDistnearestID[0] > distance - root.getRadius()) {
					NearestNeighborSearchBall(point, aIndexNode, dimension, itemMatrix, minDistnearestID);
				}
			}
		}
	}
	
	/*
	 * get the count of points in the tree.
	 */
	public int getcount(indexNode root) {
		if(root == null)
			return 0;
		if(root.isLeaf()) {	
			return root.getpointIdList().size();
		}
		else {
			Set<indexNode> listnode = root.getNodelist();
			int max = 0;
			for(indexNode aIndexNode: listnode) {
				max += getcount(aIndexNode);
			}
			return max;
		}
	}
	
	/*
	 * get the highest weight of the tree
	 */
	public int getHeight(indexNode root) {
		if(root.isLeaf())
			return 0;
		else {
			Set<indexNode> listnode = root.getNodelist();
			int max = 0;
			for(indexNode aIndexNode: listnode) {
				if(getHeight(aIndexNode)>max) {
					max = getHeight(aIndexNode);
				}
			}
			return max+1;
		}
	}
	
	
	/*
	 * get the radius of the leaf node tree and leaf number, divided by the maximum radius, used for dataset search
	 */
	public int getLeafRadius(indexNode root, ArrayList<Double> radius, double maximum, ArrayList<Double> distanceToFather, ArrayList<Double> numPoint, ArrayList<Double> nodedept, int depth) {
		if(root.isLeaf()) {
			radius.add(root.getRadius()/maximum);
			distanceToFather.add(root.getDisFather()/maximum);
			numPoint.add((double)(root.getTotalCoveredPoints()));
			nodedept.add((double)depth);
			return 1;
		}else {
			Set<indexNode> listnode = root.getNodelist();
			int max = 0;
			for(indexNode aIndexNode: listnode) {
				max += getLeafRadius(aIndexNode, radius, maximum,distanceToFather,numPoint, nodedept, depth+1);
			}
			return max;
		}
	}
	
}
