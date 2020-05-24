package au.edu.rmit.trajectory.clustering;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import au.edu.rmit.mtree.*;
import au.edu.rmit.mtree.tests.Data;
import au.edu.rmit.mtree.tests.DataDouble;
import au.edu.rmit.mtree.utils.Pair;
import au.edu.rmit.mtree.utils.Utils;
public class TrajectoryMtree extends MTree<Data>  {
	private static final PromotionFunction<Data> nonRandomPromotion1 = new PromotionFunction<Data>() {
		@Override
		public Pair<Data> process(Set<Data> dataSet, DistanceFunction<? super Data> distanceFunction) {
			return Utils.minMax(dataSet);
		}
	};	
	
	/*
	 * initilize our trajectory functions, we can specify the capacity and distance function
	 */
	public TrajectoryMtree() {//DEFAULT_MIN_NODE_CAPACITY
		super(20, DistanceFunctions.EBD, 
			new ComposedSplitFunction<Data>(
				nonRandomPromotion1,
				new PartitionFunctions.BalancedPartition<Data>()
			)
		);
	}
	
	/*
	 * initilize our trajectory functions, we can specify the capacity and distance function
	 */
	public TrajectoryMtree(int capacity, int k) {//DEFAULT_MIN_NODE_CAPACITY
		super(capacity, DistanceFunctions.EUCLIDEAN, 
			new ComposedSplitFunction<Data>(
				nonRandomPromotion1,
				new PartitionFunctions.BalancedPartition<Data>()
			)
		);
	}
	
	/*
	 * initilize our trajectory functions, we can specify the capacity and distance function
	 */
	public TrajectoryMtree(int capacity) {//DEFAULT_MIN_NODE_CAPACITY
		super(capacity, DistanceFunctions.EBD, 
			new ComposedSplitFunction<Data>(
				nonRandomPromotion1,
				new PartitionFunctions.BalancedPartition<Data>()
			)
		);
	}

	public void add(Data data) {
		super.add(data);
		_check();
	}

	public boolean remove(Data data) {
		boolean result = super.remove(data);
		_check();
		return result;
	}
	
	DistanceFunction<? super Data> getDistanceFunction() {
		return distanceFunction;
	}
	
	public void buildMtree(int [] trajectory, int traid) {
		Data data = new Data(trajectory, traid);
		add(data);
	}
	
	public void buildHistogram() {
		buildHistogram(this.root);
	}
	
	public void writeMtree(String folder) {
		writeMtreetoFile(this.root, 1, this.distanceFunction.getID(this.root.getData()), folder);
	}
	
	/*
	 * grouping trajectories into groups and write into files
	 */
	public void groupTrajectoryPivot(Map<Integer, int[]> datamap, Map<Integer, ArrayList<Integer>> edgeIndex,
			Map<Integer, Integer> edgeHistogram) {
		
	}
}
