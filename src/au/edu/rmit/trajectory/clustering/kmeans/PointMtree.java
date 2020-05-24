package au.edu.rmit.trajectory.clustering.kmeans;

import java.util.Set;
import au.edu.rmit.mtree.*;
import au.edu.rmit.mtree.tests.DataDouble;
import au.edu.rmit.mtree.utils.Pair;
import au.edu.rmit.mtree.utils.Utils;

public class PointMtree extends MTree<DataDouble>  {
	private static final PromotionFunction<DataDouble> nonRandomPromotion1 = new PromotionFunction<DataDouble>() {
		@Override
		public Pair<DataDouble> process(Set<DataDouble> dataSet, DistanceFunction<? super DataDouble> distanceFunction) {
			return Utils.minMax(dataSet);
		}
	};	
	
	/*
	 * initilize our trajectory functions, we can specify the capacity and distance function
	 */
	public PointMtree() {//DEFAULT_MIN_NODE_CAPACITY
		super(20, DistanceFunctions.EBD, 
			new ComposedSplitFunction<DataDouble>(
				nonRandomPromotion1,
				new PartitionFunctions.BalancedPartition<DataDouble>()
			)
		);
	}
	
	/*
	 * initilize our trajectory functions, we can specify the capacity and distance function
	 */
	public PointMtree(int capacity) {//DEFAULT_MIN_NODE_CAPACITY
		super(capacity, DistanceFunctions.EUCLIDEAN, 
			new ComposedSplitFunction<DataDouble>(
				nonRandomPromotion1,
				new PartitionFunctions.BalancedPartition<DataDouble>()
			)
		);
	}


	public void add(DataDouble data) {
		super.add(data);
		_check();
	}

	public boolean remove(DataDouble data) {
		boolean result = super.remove(data);
		_check();
		return result;
	}
	
	DistanceFunction<? super DataDouble> getDistanceFunction() {
		return distanceFunction;
	}
	
	
	public void buildMtree(double [] trajectory, int traid) {
		DataDouble data = new DataDouble(trajectory, traid);
		add(data);
	}
	
	public void buildHistogram() {
		buildHistogram(this.root);
	}
	
	public void writeMtree(String folder) {
		writeMtreetoFile(this.root, 1, this.distanceFunction.getID(this.root.getData()), folder);
	}
}
