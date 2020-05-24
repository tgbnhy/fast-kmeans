package au.edu.rmit.trajectory.clustering.kmeans;

/*
 * this class is for the heap method
 */
public class kmeansHeap implements Comparable<kmeansHeap>{
	double bound;
	int idx;
	
	public double getBound() {
		return bound;
	}
	
	public int getidx() {
		return idx;
	}
	
	public kmeansHeap(int idx, double bound) {
		this.idx = idx;
		this.bound = bound;
	}
	
	@Override
    public int compareTo(kmeansHeap other) {
        return (int)(this.getBound() - other.getBound());
    }
	
}
