package au.edu.rmit.trajectory.clustering.kmeans;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import au.edu.rmit.trajectory.clustering.kpaths.Util;
import au.edu.rmit.trajectory.clustering.streaming.DataReading;
import edu.wlu.cs.levy.cg.KeyDuplicateException;
import edu.wlu.cs.levy.cg.KeySizeException;

public class plotData {

	public plotData() {
		// TODO Auto-generated constructor stub
	}

	
	public static void savePlotDataFile(String folderName, int option, int k, String indexType, String output) throws IOException {
		File folder = new File(folderName);
		int datasetnum = 11;		
		String []content = new String[datasetnum];// dataset
		int dataCounter = 0;
		for (final File fileEntry : folder.listFiles()) {
	        if (!fileEntry.isDirectory()) {
	        	String fileName = fileEntry.getName();
	        	if(fileName.contains(Integer.toString(k)+"_"+indexType)) {
	        		if(fileName.contains("Bigcross")) {
	        			dataCounter = 0;
	        		}else if(fileName.contains("conf")){
	        			dataCounter = 1;
	        		}else if(fileName.contains("covt")) {
	        			dataCounter = 2;
	        		}else if(fileName.contains("euro")) {
	        			dataCounter = 3;
	        		}else if(fileName.contains("KeggDirect")) {
	        			dataCounter = 4;
	        		}else if(fileName.contains("KeggUndirect")) {
	        			dataCounter = 5;
	        		}else if(fileName.contains("NYC")) {
	        			dataCounter = 6;
	        		}else if(fileName.contains("Skin")) {
	        			dataCounter = 7;
	        		}else if(fileName.contains("power")) {
	        			dataCounter = 8;
	        		}else if(fileName.contains("Spatial")) {
	        			dataCounter = 9;
	        		}else if(fileName.contains("USCensus")) {
	        			dataCounter = 10;
	        		}
	        		content[dataCounter] = readContent(fileEntry, option);
	        		
	        	}	        	
	            System.out.println(fileEntry.getName());
	        }
	    }
		for(int i = 1;  i<=datasetnum; i++)
			Util.write(folderName+"/plot/"+output, Integer.toString(i)+ "\t"+ content[i-1]+"\n");
	}
	
    //read the lines
    public static String readContent(File file, int lineCountre) throws IOException{
        System.out.println("read file " + file.getCanonicalPath() );
        String content = "";
        try(BufferedReader br  = new BufferedReader(new FileReader(file))){
              String strLine;
              int counter = 0;
              while((strLine = br.readLine()) != null){
            	  if(counter++ == lineCountre) {
            		  content = strLine;
            		  break;
            	  }
              }
        }
        return content;
    }
    
    //read the specific column in lines
    public static String readContent(File file, int lineCountre, int start, int end) throws IOException{
        System.out.println("read file " + file.getCanonicalPath() );
        String content = "";
        try(BufferedReader br  = new BufferedReader(new FileReader(file))){
              String strLine;
              int counter = 0;
              while((strLine = br.readLine()) != null){
            	  if(counter++ == lineCountre) {
            		  String columns[] = strLine.split("\t");            		  
            		  for(int i = start; i<=end; i++)
            			  content += columns[i-1]+"\t";
            		  break;
            	  }
              }
        }
        return content;
    }
    
    // export the value that change with k.    
    static void extractk(String folderName, String datasetName, int option) throws IOException {
    	File folder = new File(folderName);
		int kvaluenum = 3;		
		String []content = new String[kvaluenum];// dataset
		for (final File fileEntry : folder.listFiles()) {
	        if (!fileEntry.isDirectory()) {
	        	String fileName = fileEntry.getName();
	        	if(fileName.contains(datasetName)) {
	        		if(fileName.contains("10_Ball")) {
	        			content[0] =  readContent(fileEntry, option);
	        		}else if(fileName.contains("100_Ball")) {
	        			content[1] =  readContent(fileEntry, option);
	        		}else if(fileName.contains("1000_Ball")) {
	        			content[2] =  readContent(fileEntry, option);
	        		}
	        	}
	        }
		}
		for(int i = 1;  i<=kvaluenum; i++)
			Util.write(folderName+"/plot/"+datasetName+"IncreasingK.txt", Integer.toString(i)+ "\t"+ content[i-1]+"\n");
    }
    
    // export the value that change with k.    
    static void extractDataScale(String folderName, String datasetName, int option, int k) throws IOException {
    	File folder = new File(folderName);
		int kvaluenum = 5;		
		String []content = new String[kvaluenum];// dataset
		for (final File fileEntry : folder.listFiles()) {
	        if (!fileEntry.isDirectory()) {
	        	String fileName = fileEntry.getName();	  
	        	String scaleset[] = fileName.split("_");
	        	if(fileName.contains(datasetName) && !fileName.contains("index") && scaleset.length>1) {
	        		String scale = scaleset[1];
	        		if(scale.equals("100")) {
	        			content[0] =  readContent(fileEntry, option);
	        		}else if(scale.equals("1000")) {
	        			content[1] =  readContent(fileEntry, option);
	        		}else if(scale.equals("10000")) {
	        			content[2] =  readContent(fileEntry, option);
	        		}else if(scale.equals("100000")) {
	        			content[3] =  readContent(fileEntry, option);
	        		}else if(scale.equals("1000000")) {
	        			content[4] =  readContent(fileEntry, option);
	        		}
	        	}
	        }
		}
		for(int i = 1;  i<=kvaluenum; i++)
			Util.write(folderName+"/plot/"+datasetName+"IncreasingData"+Integer.valueOf(k)+".txt", Integer.toString(i)+ "\t"+ content[i-1]+"\n");
    }
    
    /*
     * 
     */
    static void extrackIndex(String folderName, String datasetName, int option, int k) throws IOException {
    	File folder = new File(folderName);
		int kvaluenum = 3;		
		String []content = new String[kvaluenum];// dataset
		for (final File fileEntry : folder.listFiles()) {
	        if (!fileEntry.isDirectory()) {
	        	String fileName = fileEntry.getName();	  
	        	String scaleset[] = fileName.split("_");
	        	if(fileName.contains(datasetName) && scaleset.length>2) {
	        		String scale = scaleset[4];
	        		if(scaleset[3].equals(Integer.toString(k))) {
	        		if(scale.equals("BallMetric")) {
	        			content[0] =  readContent(fileEntry, option);
	        		}else if(scale.equals("Mtree")) {
	        			content[1] =  readContent(fileEntry, option);
	        		}else if(scale.equals("hkmtree")) {
	        			content[2] =  readContent(fileEntry, option);
	        		}
	        		}
	        	}
	        }
		}
		for(int i = 1;  i<=kvaluenum; i++)
			Util.write(folderName+"/plot/"+datasetName+"IndexType"+Integer.valueOf(k)+".txt", Integer.toString(i)+ "\t"+ content[i-1]+"\n");
    }
    
    /*
     * test
     */
    static void extracDimension(String folderName, String datasetName, int option, int k) throws IOException {
    	File folder = new File(folderName);
		int kvaluenum = 5;		
		String []content = new String[kvaluenum];// dataset
		for (final File fileEntry : folder.listFiles()) {
	        if (!fileEntry.isDirectory()) {
	        	String fileName = fileEntry.getName();	  
	        	String scaleset[] = fileName.split("_");
	        	if(fileName.contains(datasetName) && scaleset.length>2) {
	        		String scale = scaleset[2];
					if (scaleset[3].equals(Integer.toString(k))) {
						if (scale.equals("10")) {
							content[0] = readContent(fileEntry, option);
						} else if (scale.equals("20")) {
							content[1] = readContent(fileEntry, option);
						} else if (scale.equals("30")) {
							content[2] = readContent(fileEntry, option);
						} else if (scale.equals("40")) {
							content[3] = readContent(fileEntry, option);
						} else if (scale.equals("50")) {
							content[4] = readContent(fileEntry, option);
						}
					}
	        	}
	        }
		}
		for(int i = 1;  i<=kvaluenum; i++)
			Util.write(folderName+"/plot/"+datasetName+"Dimension"+Integer.valueOf(k)+".txt", Integer.toString(i)+ "\t"+ content[i-1]+"\n");
    }
    
    /*
     * extract the running time of different index with the increasing of capacity,
     */
    static void extrackIndexCapacity(String folderName, String datasetName, int option, int k) throws IOException {
    	File folder = new File(folderName);
		int kvaluenum = 6;		// capacity number
		String []content = new String[kvaluenum];// dataset
		String []content_ball = new String[kvaluenum];// dataset
		for (final File fileEntry : folder.listFiles()) {
	        if (!fileEntry.isDirectory()) {
	        	String fileName = fileEntry.getName();	  
	        	String scaleset[] = fileName.split("_");
	        	if(fileName.contains(datasetName) && scaleset.length>2) {
	        		String scale = scaleset[4];
	        		if(scaleset[3].equals(Integer.toString(k))) {
						if (fileName.contains("10.log")) {
							if (scale.equals("HKT"))
								content[0] = readContent(fileEntry, option, 14, 16);
							else if (scale.equals("BallMetric")) {
								content_ball[0] = readContent(fileEntry, option, 14, 16);
							}
						} else if (fileName.contains("20.log")) {
							if (scale.equals("HKT"))
								content[1] = readContent(fileEntry, option, 14, 16);
							else if (scale.equals("BallMetric")) {
								content_ball[1] = readContent(fileEntry, option, 14, 16);
							}
						}else if (fileName.contains("30.log")) {
							if (scale.equals("HKT"))
								content[2] = readContent(fileEntry, option, 14, 16);
							else if (scale.equals("BallMetric")) {
								content_ball[2] = readContent(fileEntry, option, 14, 16);
							}
						}else if (fileName.contains("40.log")) {
							if (scale.equals("HKT"))
								content[3] = readContent(fileEntry, option, 14, 16);
							else if (scale.equals("BallMetric")) {
								content_ball[3] = readContent(fileEntry, option, 14, 16);
							}
						}else if (fileName.contains("50.log")) {
							if (scale.equals("HKT"))
								content[4] = readContent(fileEntry, option, 14, 16);
							else if (scale.equals("BallMetric")) {
								content_ball[4] = readContent(fileEntry, option, 14, 16);
							}
						}else if (fileName.contains("60.log")) {
							if (scale.equals("HKT"))
								content[5] = readContent(fileEntry, option, 14, 16);
							else if (scale.equals("BallMetric")) {
								content_ball[5] = readContent(fileEntry, option, 14, 16);
							}
						}	        		
	        		}
	        	}
	        }
		}
		for(int i = 1;  i<=kvaluenum; i++)
			Util.write(folderName+"\\plot\\"+datasetName+"Capacity"+Integer.valueOf(k)+".txt", Integer.toString(i)+ "\t"+ content[i-1]+content_ball[i-1]+"\n");
    }
    
    
    
    
    // option 0: total running time, 1: ratio, 
 	// 3: assignment time, 4, ratio
 	// 6: refinement, 7, raito
 	// 9: computations, 10, ratio
 	// 12: boundCompares, 14: dataAccess, 16: boundUpdates, 18: memoryUsage
 	//read existing files and store the data in a folder, then add the new ran data to update the score.
    public static void main(String[] args) throws IOException, KeySizeException, KeyDuplicateException {
		for (int kvalue = 10; kvalue <= 1000; kvalue = kvalue * 10) {
			savePlotDataFile(args[0], 0, kvalue, "Ball", "runingtime" +Integer.valueOf(kvalue)+ ".txt");
			savePlotDataFile(args[0], 1, kvalue, "Ball", "runingtimeratio" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 3, kvalue, "Ball", "assigntime" +Integer.valueOf(kvalue)+ ".txt");
			savePlotDataFile(args[0], 4, kvalue, "Ball", "assigntimeratio" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 6, kvalue, "Ball", "refinetime" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 7, kvalue, "Ball", "refineratio" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 9, kvalue, "Ball", "computation" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 10, kvalue, "Ball", "computationratio" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 12, kvalue, "Ball", "boundCompares" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 14, kvalue, "Ball", "dataAccess" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 16, kvalue, "Ball", "boundUpdates" + Integer.valueOf(kvalue)+".txt");
			savePlotDataFile(args[0], 18, kvalue, "Ball", "memoryUsage" + Integer.valueOf(kvalue)+".txt");
		}
		
		extractk(args[0], "road", 1);
		
		extractDataScale(args[0], "NYC", 0, 100);
		
		extrackIndex(args[0], "road", 1, 100);
		
		extrackIndexCapacity(args[0], "NYC", 0, 10);// get the running time, all the index method.
		
		extracDimension(args[0], "covt", 0, 10);// get the running time, all the index method.
		extracDimension(args[0], "bigcross", 0, 10);// get the running time, all the index method.		
    }
}
