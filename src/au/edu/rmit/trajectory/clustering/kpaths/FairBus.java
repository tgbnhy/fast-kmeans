package au.edu.rmit.trajectory.clustering.kpaths;
import org.jgrapht.alg.connectivity.*;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import org.openstreetmap.osmosis.core.container.v0_6.EntityContainer;
import org.openstreetmap.osmosis.core.container.v0_6.NodeContainer;
import org.openstreetmap.osmosis.core.container.v0_6.RelationContainer;
import org.openstreetmap.osmosis.core.container.v0_6.WayContainer;
import org.openstreetmap.osmosis.core.domain.v0_6.Tag;
import org.openstreetmap.osmosis.core.domain.v0_6.Way;
import org.openstreetmap.osmosis.core.task.v0_6.Sink;
import crosby.binary.osmosis.OsmosisReader;

public class FairBus implements Sink {
	
	/*
	 * compute the connectivity for a given path in the network.
	 */
	void computeConnectivity(String file) {
		
		int routeNumber = 0;
		// read every route, and store
		//
	//	for(int i=0; i< routeNumber; i++)
	//		Util.Intersection(arr1, arr2);
		
	}
	
	/*
	 * compute the lower bound of number of edges needs to be added to meet the connectivity constraint
	 */
	void computeLowerBoundNumberEdges() {
		
	}
	
	
	/*
	 * read the file and construct the transport network and road network from openstreetmap
	 */
	void constructGraph() {
		//construct the graph based on 
	}
	
	/*
	 * align bus stops to road network edge, and 
	 */
	void alignBusstopsToEdge(String stopFile, String routeFile, Map<Integer, String> edgeInfo) {
		Map<Integer, String> node = new HashMap<Integer, String>();// read from 
		Map<Integer, Set<Integer>> waysMap = new HashMap<Integer, Set<Integer>>();
		Map<Integer, ArrayList<String>> busroutesMap = new HashMap<Integer, ArrayList<String>>();
		
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(stopFile)));			
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the point by their id.
				String str = in.nextLine();
				if(str.contains("node")) {//this is the bus stops, 
					//must contain tags, stop at way "107430881"
				}
				if(str.contains("way")) {// this is train station.
					//read the first node,
				}
				
			}
			in.close();
		}		
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		try {// read all the ways,
			Scanner in = new Scanner(new BufferedReader(new FileReader(routeFile)));			
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the point by their id.
				String str = in.nextLine();
				if(str.contains("way")) {
					//read the first node,
				}
				
			}
			in.close();
		}		
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		try {// read all the relations, and call the node and ways to output the route.
			Scanner in = new Scanner(new BufferedReader(new FileReader(routeFile)));			
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the point by their id.
				String str = in.nextLine();
				if(str.contains("routes")) {
					//get all the nodes and way
				}
				
			}
			in.close();
		}		
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		//construct the transit network
		
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		//clean the data
		 InputStream inputStream = new FileInputStream("/Users/sw160/Desktop/transit_network/porto_transit_stops");
	        OsmosisReader reader = new OsmosisReader(inputStream);
	        reader.setSink(new FairBus());
	        reader.run();
	}
	
	
	void stopPairMatching() {
	
	}
	
	/*
	 * map stations to edge
	 */
	void mapStationsToEdge(String edgefile) {
		//assign the stations to the right edge
		//
	}
	
	void mapTransitRouteToPath(String edgefile) {
		//transfer each route to a set of edge in the road network
		//compute the distance with edge
		
	}

	@Override
	public void initialize(Map<String, Object> metaData) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void complete() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void close() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void process(EntityContainer entityContainer) {
		// TODO Auto-generated method stub
		if (entityContainer instanceof NodeContainer) {
            // Nothing to do here
        } else if (entityContainer instanceof WayContainer) {
            Way myWay = ((WayContainer) entityContainer).getEntity();
            for (Tag myTag : myWay.getTags()) {
                if ("highway".equalsIgnoreCase(myTag.getKey())) {
                    System.out.println(" Woha, it's a highway: " + myWay.getId());
                    break;
                }
            }
        } else if (entityContainer instanceof RelationContainer) {
            // Nothing to do here
        } else {
            System.out.println("Unknown Entity!");
        }
	}
}
