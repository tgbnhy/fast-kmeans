package cleaning;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.Scanner;

import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;
import scala.collection.generic.BitOperations.Int;

public class traclus {

	public traclus() {
		// TODO Auto-generated constructor stub
	}


	
	/* 
	 * This function clean and produce data with specific scale.
	 * 
	 */
	public static void cleanForTraClusInput(String folder, String output, int num) {		
		output = output+ Integer.toString(num);
		Util.write(output, "2\n"+num+"\n");
		int counter = 0;
		int length = 0; 
		String content = "";
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(folder)));
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by
				String str = in.nextLine();				
				if(str.equals("0 0 0")) {
					content = Integer.toString(counter) + " " + Integer.valueOf(length) + content;
					Util.write(output, content+"\n");
					counter ++;
					length = 0;
					content = "";
				} else {
					String strr = str.trim();
					String[] abc = strr.split(" ");
					content += " " + abc[0]+ " " + abc[1];
					length++;
				}
				if(counter>num)
					break;
			}
			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	
	
	public static void main(String[] args)  {
	//	convertTrajectory(args[0], args[1]);
	//	convertLineString(args[0], args[1], 10);
	//	for(int i=1000000; i<=1000000; i=i*10) {
		//	if(i==100000)
			//	continue;
		//	cleanPorto(args[0], args[1], i);
		//}
		for(int i=100; i<=100000; i=i*10)
			cleanForTraClusInput(args[0], args[1], i);
	}
}
