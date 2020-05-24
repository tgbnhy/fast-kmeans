package cleaning;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class nantong {

	public nantong() {
		// TODO Auto-generated constructor stub
	}
	public static void main(String[] args) {
		//	construct_graph(args[0], args[1], args[2], args[3],args[4]);
		//	generate_mapv_graph(args[0], args[1]);
		//	generate_mapv_graph_porto(args[0], args[1]);
		//	generate_mapv1_graph(args[0], args[1]);
		//	generate_mapv_graph_nantong(args[0], args[1]);
			generate_trajectory(args[0], args[1], args[2], args[3]);
		}
		/*
		 * generate the graph: edge, vertex1, vertex2.
		 */
		public static void construct_graph(String vertex, String edge, String output1, String dataset, String output2) {
			Map<String, String> vertex_id = new HashMap();
			Map<String, Integer> originid_new = new HashMap();
			int i=0;
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(vertex)));
				while (in.hasNextLine()) {
					String str = in.nextLine();
					String strr = str.trim();
					String[] abc = strr.split(";");
					System.out.println(abc[1]+"-"+abc[2]);
					vertex_id.put(abc[1]+"-"+abc[2], abc[0]+";"+i);
					originid_new.put((abc[0]), i++);
				}
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			Map<String, String> edge_vertexpair = new HashMap();
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(edge)));
				while (in.hasNextLine()) {
					String str = in.nextLine();
					String strr = str.trim();
					String[] abc = strr.split(";");
					String[] lats = abc[1].split(",");
					String[] longs = abc[2].split(",");
					String start = lats[0]+"-"+longs[0];
				//	System.out.println(start);
					String end = lats[lats.length-1]+"-"+longs[longs.length-1];
					String start_id="", end_id="";
					if(vertex_id.containsKey(start)) {
						 start_id = vertex_id.get(start);
						 start_id = start_id.split(";")[0];
					}
					if(vertex_id.containsKey(end)) {
						end_id = vertex_id.get(end);
						end_id = end_id.split(";")[0];
					}
					edge_vertexpair.put(abc[0], start_id+"\t"+end_id);
					write(output1, abc[0]+"\t"+start_id+"\t"+end_id+"\n");// generate the graph
				}
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			/*
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(dataset)));
				while (in.hasNextLine()) {
					String str = in.nextLine();
					String strr = str.trim();
					String[] abc = strr.split("\t");
					if(abc.length<2)
						continue;
					String[] edges = abc[1].split(",");
					String content = abc[0]+"\t";
					for(int edges_i=0; edges_i<edges.length; edges_i++) {
						int old_start = 0, old_end = 0;
						int label = 0;
						if(edge_vertexpair.containsKey(edges[edges_i])) {
							String vertex_pair = edge_vertexpair.get(edges[edges_i]);
							String a[] = vertex_pair.split("\t");
							System.out.println(a[0]);
							int start_id = Integer.valueOf(a[0]);
							int end_id = Integer.valueOf(a[1]);
							if(label == 1 && old_end == start_id) {
								content += Integer.toString(start_id)+",";
							}else if(label == 0){
								content += Integer.toString(start_id)+",";
							}
							label = 1;
							old_start = start_id;
							old_end = end_id;
						}
						content += "\n";
						write(output2, content);
					}
				}
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}*/
		}
		
		
		//C:\mapv-2.0.12\travis\data
		public static void basic_permutation() {
			
		}
		
		//write the information into files
		public static void write(String fileName, String content) {   
	        RandomAccessFile randomFile = null;  
	        try {       
	            randomFile = new RandomAccessFile(fileName, "rw");     
	            long fileLength = randomFile.length();
	            randomFile.seek(fileLength);     
	            randomFile.writeBytes(content);      
	        } catch (IOException e) {     
	            e.printStackTrace();     
	        } finally{  
	            if(randomFile != null){  
	                try {  
	                    randomFile.close();  
	                } catch (IOException e) {  
	                    e.printStackTrace();  
	                }  
	            }  
	        }  
	    }
		/*
		 * generate the graph for mapv based on the edge graph data
		 */
		public static void generate_mapv_graph_porto(String edge, String output) {
			write(output, "geometry\n");
			Map<String, String> edge_vertexpair = new HashMap();
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(edge)));
				while (in.hasNextLine()) {
					String content = "\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [";
					String str = in.nextLine();
					String strr = str.trim();
					String[] abc = strr.split(",");
					for(int i=0; i<abc.length/2;i++) {
						content+= "["+abc[i+abc.length/2]+","+abc[i]+"]";					
						if(i<abc.length/2-1)
							content+=",";
					}
					content+="]}\"";
					write(output, content+"\n");
				}
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		/*
		 * we generate the trajectories
		 */
		static void generate_trajectory(String streamfile, String output, String output1, String mapping) {
			//read the edge mapping
			Map<Integer, Integer> oldedgeNew = new HashMap<>();
			Map<String, Integer> cameridNewid = new HashMap<>();
			int count=0;
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(mapping)));
				while (in.hasNextLine()) {
					String str = in.nextLine();
					String strr = str.trim();
					String[] strr1 = str.split(",");
					int oldEdge = Integer.valueOf(strr1[0]);
					String[] aStrings = strr1[1].split("_");
					String id = aStrings[0]+"_"+aStrings[1];
					if(!cameridNewid.containsKey(id)) {
						cameridNewid.put(id, count);
						oldedgeNew.put(oldEdge, count);
						write(output, id+ "\t"+count+"\n");
						count++;
					}else {
						oldedgeNew.put(oldEdge, cameridNewid.get(id));
					}
				}
				in.close();
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			Map<Integer, ArrayList<Integer>> trajectory = new HashMap<>();
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(streamfile)));
				while (in.hasNextLine()) {
					String str = in.nextLine();
					String strr = str.trim();
					String[] strr1 = str.split(",");
					if(strr1[1].equals("13") || strr1[1].equals("3")) {
						continue;
					}
					int trajectory_id = Integer.valueOf(strr1[1]);
					int newEdge = oldedgeNew.get(Integer.valueOf(strr1[2]));
					ArrayList<Integer> trArrayList;
					if(trajectory.containsKey(trajectory_id)) {
						trArrayList = trajectory.get(trajectory_id);
					}else {
						trArrayList = new ArrayList<>();
					}
					trArrayList.add(newEdge);
					trajectory.put(trajectory_id, trArrayList);
				}
				in.close();
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			for(int traid: trajectory.keySet()) {
				ArrayList<Integer> traarray = trajectory.get(traid);
				String content ="";
				for(int id: traarray) {
					content += id+",";
				}
				write(output1, traid+ "\t"+content+"\n");
			}
		}
		
		/*
		 * generate the graph for mapv based on the edge graph data
		 */
		public static void generate_mapv_graph_nantong(String edge, String output) {
			write(output, "geometry\n");
			Map<String, String> edge_vertexpair = new HashMap();
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(edge)));
				while (in.hasNextLine()) {
					String content = "\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [";
					String str = in.nextLine();
					String strr = str.trim();
					String[] strr1 = str.split(";");
					String temp = strr1[1]+","+strr1[2];
					String[] abc = temp.split(",");
					if(abc.length<3)
						continue;
					for(int i=0; i<abc.length/2;i++) {
						content+= "["+abc[i+abc.length/2]+","+abc[i]+"]";					
						if(i<abc.length/2-1)
							content+=",";
					}
					content+="]}\"";
					write(output, content+"\n");
				}
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		/*
		 * generate the graph for mapv based on the edge graph data
		 */
		public static void generate_mapv_graph(String edge, String output) {
			write(output, "geometry\n");
			Map<String, String> edge_vertexpair = new HashMap();
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(edge)));
				while (in.hasNextLine()) {
					String content = "\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [";
					String str = in.nextLine();
					String strr = str.trim();
					String[] abc = strr.split(",");
					for(int i=0; i<abc.length/2;i++) {
					//	content+= "["+abc[i+abc.length/2]+","+abc[i]+"]";
						
						double a = (Double.valueOf(abc[i+abc.length/2]) + 0.012557898);//  this is the standard normalization for Beijing
						double b= Double.valueOf(abc[i]) + 0.0077440623;
						content+=  "["+Double.toString(a)+","+Double.toString(b)+"]";
						
						if(i<abc.length/2-1)
							content+=",";
					}
					content+="]}\"";
					write(output, content+"\n");
				}
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		/*
		 * generate the graph for mapv-"baidu-map-polyline-time" based on the edge graph data
		 */
		public static void generate_mapv1_graph(String edge, String output) {
		//	write(output, "geometry\n");
			Map<String, String> edge_vertexpair = new HashMap();
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(edge)));
				while (in.hasNextLine()) {
					String content = "";
					String str = in.nextLine();
					String strr = str.trim();
					String[] abc = strr.split(",");
					for(int i=0; i<abc.length/2;i++) {
						content+= abc[i+abc.length/2]+","+abc[i];
						if(i<abc.length/2-1)
							content+=",";
					}
				//	content+="]}\"";
					write(output, content+"\n");
				}
			}
			catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}

}
