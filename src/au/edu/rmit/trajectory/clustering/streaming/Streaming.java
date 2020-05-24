package au.edu.rmit.trajectory.clustering.streaming;

import java.io.IOException;
import java.sql.Time;
import java.util.*;

import scala.collection.immutable.HashSet;


public class Streaming {
    private Simulator simulator;
    public Map<Integer, Set<Integer>> edgeInvertedIndex;
    public Map<Integer, int[]> trajectories;

    /**
     * @param startFrom we need to prepare some data before testing streaming. if startFrom = 3600, it means we load 3600s data
     *                  before the actual streaming testing
     * @param windowSize sliding window in time dimension
     */
    public Streaming(int startFrom, int windowSize, int speedup, String datafile) throws IOException {
        simulator= new Simulator(startFrom, windowSize, speedup, datafile);
    }

    /**
     * Start streaming simulation.
     * the simulator object keeps a dynamic List<Tuple>
     */
    public Streaming start(){
        simulator.start();
        return this;
    }

    /**
     * The stream which is continually updated in each second is simulator.stream
     * The trajectories and edgeInvertedIndex are not updated on a second basis.
     * When the clustering finish previous job and is ready for the next job,
     * you need call this method to prepare edgeInvertedIndex and trajectories from simulator.stream
     * Then you could offer the up-to-date trajectories and edgeInvertedIndex to clustering algorithm
     */
    public void updateIndexes() {
    	long time1 = System.nanoTime();
        edgeInvertedIndex = new HashMap<>();
        Map<Integer, Set<Integer>> trajs = new HashMap<>();
        for (Tuple t : simulator.stream){
            //update edge inverted index
            if (!edgeInvertedIndex.containsKey(t.edgeId))
                edgeInvertedIndex.put(t.edgeId, new TreeSet<Integer>());
            edgeInvertedIndex.get(t.edgeId).add(t.carId);
            //update trajectory data
            if (!trajs.containsKey(t.carId))
                trajs.put(t.carId, new TreeSet<>());
            trajs.get(t.carId).add(t.edgeId);
        }

        trajectories = new HashMap<>();
        for (Map.Entry<Integer, Set<Integer>> entry : trajs.entrySet())
            trajectories.put(entry.getKey(), entry.getValue().stream().mapToInt(i->i).toArray());
        long time2 = System.nanoTime();
    //    System.err.println(trajs.size());
        System.err.println((time2-time1)/1000000000.0);
    }

    /**
     * This must be called before new another Streaming object
     */
    public void close(){
        try {
            simulator.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
