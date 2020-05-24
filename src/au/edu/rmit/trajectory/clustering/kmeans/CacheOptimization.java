package au.edu.rmit.trajectory.clustering.kmeans;

import java.util.concurrent.TimeUnit;

import com.github.benmanes.caffeine.cache.Cache;
import com.github.benmanes.caffeine.cache.Caffeine;

public class CacheOptimization {
	Cache<Integer, double[]> cache;
	public CacheOptimization() {
		// TODO Auto-generated constructor stub
	}
	
	public void CacheAllCentorid(int k, double centroid[][]) {
		cache = Caffeine.newBuilder().expireAfterWrite(10, TimeUnit.MINUTES).maximumSize(10_000)
				.build();
		for(int i=0; i<k; i++) {
			cache.put(i, centroid[i]);
		}
	}
	
	public double[] CacheAccess(int key) {
		double[] graph = cache.getIfPresent(key);
		return graph;
	}
	
	public void comparedCache(int testtime, int k, double centroid[][]) {
		long startTime1 = System.nanoTime();
		for(int i = 0; i<testtime; i++) {
			for(int j=0; j<k; j++) {
				double graph[] = centroid[j];
			}
		}	
		long endtime = System.nanoTime();
		double time = (endtime-startTime1)/1000000000.0;
		System.out.print("cache time compare: "+time+" ");
		startTime1 = System.nanoTime();
		for(int i = 0; i<testtime; i++) {
			for(int j=0; j<k; j++) {
				double graph[] = cache.getIfPresent(j);
			}
		}	
		endtime = System.nanoTime();
		time = (endtime-startTime1)/1000000000.0;
		System.out.println(time);
	}
	
	public void main() {
		
	}
	// cache, use caffeine to store the centroids, we also show the data and bound access is important component to the distance.
	// complexity analysis
}
