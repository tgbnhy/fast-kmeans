package edu.nyu.unik.expriments;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import kafka.common.GenerateBrokerIdException;

public class kcenter {

	public kcenter() {
		// TODO Auto-generated constructor stub
	}
	
	
	public static void readRawData(String filename, String pointfile, String output) throws FileNotFoundException {
		PrintStream fileOut = new PrintStream(output);
		System.setOut(fileOut);
		Map<String, double[]> points = new HashMap<>();
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(pointfile)));					
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split(";");
				double coordinate[] = new double[2];
				coordinate[0] = Double.valueOf(abc[1]);
				coordinate[1] = Double.valueOf(abc[2]);
				points.put(abc[0], coordinate);
			}
			in.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		int timestamp = 0;
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(filename)));					
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split("\t");
				if(abc.length==1)
					continue;
					
				System.out.print(timestamp+"\t");
				String[] vertexids = abc[1].split(",");
				System.out.print(vertexids.length+"\t");
				for(String vertexid: vertexids) {
					double coordinate[] = points.get(vertexid);
					System.out.print(coordinate[1]+","+coordinate[0]+"\t");
				}
				System.out.println();
				timestamp++;
			}
			in.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}	
//		timestamp tab n tab longitude_1,latitude_1 tab longitude_2,latitude_2 tab ... tab longitude_n,latitude_n
	}
	
	/*
	 * the query file for kcenter which shows the operations
	 */
	public void GenerateQueryfilekCenter() {
		
	}

	public static void main(String[] args) throws FileNotFoundException {
		readRawData(args[0], args[1], args[2]);
	}
}
