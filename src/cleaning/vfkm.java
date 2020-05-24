package cleaning;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.Scanner;

import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;
import scala.collection.generic.BitOperations.Int;

public class vfkm {

	public vfkm() {
		// TODO Auto-generated constructor stub
	}

	// 39.685082, 116.069502, 40.187159, 116.766942, Beijing,
	//  40.900374, -8.787604, 41.415925, -8.178577, Tdrive
	public static void convertTrajectory(String filename, String output) {
		Util.write(output, "-8.787604 -8.178577 40.900374 41.415925\n");
		Util.write(output, "116.069502 116.766942 39.685082 40.187159\n");
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(filename)));	
			double time = 0;
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split("\t");
				String[] points = abc[1].split("]");
				for(String point: points) {
					point = point.substring(2);
				//	point.replace('[', ' ');
				//	point.replace(']', ' ');
				//	point.replace(',', ' ');
					String latlong[] = point.split(",");
				//	System.out.println(point);
				//	point = points+" "+ Double.toString(time);
					time += 0.01;
					DecimalFormat dec = new DecimalFormat("#0.00");
					Util.write(output, latlong[0]+ " "+ latlong[1] + " "+ dec.format(time)+ "\n"); 
				}
				Util.write(output, "0 0 0\n");
			}
			in.close();
		}		
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	/* 
	 * This function converts the output of VFKM to LinString for visualization
	 * 
	 */
	public static void convertLineString(String folder, String output, int k) {
		for (int i = 0; i < k; i++) {
			String content = "\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [";
			try {
				Scanner in = new Scanner(new BufferedReader(new FileReader(folder+"\\vf_r_"+i+".txt")));
				while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
					String str = in.nextLine();
					String strr = str.trim();
					String[] abc = strr.split(" ");
					if (abc.length > 1) {
						content += "[" + abc[0] + "," + abc[1] + "],";
					}
				}
				in.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			content = content.substring(0, content.length() - 1);// remove the last character
			content += "]}\"";
			Util.write(output, content + "\n");
		}
	}
	
	/* 
	 * This function clean and produce data with specific scale.
	 * 
	 */
	public static void cleanPorto(String folder, String output, int num) {
	//	Util.write(output, "-8.787604 -8.178577 40.900374 41.415925 0 765782\n");
		output = output+ Integer.toString(num);
		int counter = 0;
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(folder)));
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by
				String str = in.nextLine();
				if(str.equals("0 0 0")) {
					Util.write(output, str+"\n");
					counter++;
				} else {
					String strr = str.trim();
					String[] abc = strr.split(" ");
					if(abc.length == 6)
						Util.write(output, str+"\n");
					if (Double.valueOf(abc[0]) >= -8.787604 && Double.valueOf(abc[0]) <= -8.178577
							&& Double.valueOf(abc[1]) >= 40.900374 && Double.valueOf(abc[1]) < 41.415925) {
						Util.write(output, str+"\n");
					}
				}
				if(counter>num)
					break;
			}
			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	/* 
	 * This function clean and produce data with specific scale.
	 * 
	 */
	public static void cleanTdrive(String folder, String output, int num) {
	//	Util.write(output, "-8.787604 -8.178577 40.900374 41.415925 0 765782\n");
		output = output+ Integer.toString(num);
		int counter = 0;
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(folder)));
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by
				String str = in.nextLine();
				if(str.equals("0 0 0")) {
					Util.write(output, str+"\n");
					counter++;
				} else {
					Util.write(output, str+"\n");
				}
				if(counter> num)
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
			cleanTdrive(args[0], args[1], i);
	}
}
