package cleaning;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Scanner;

import au.edu.rmit.trajectory.clustering.kpaths.Util;

public class skin {

	public skin() {
		// TODO Auto-generated constructor stub
	}
	
	/*
	 * 
	 */
	public static void main(String[] args)  {
		try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(args[0])));
			while (in.hasNextLine()) {
				String content = "";
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split("\t");
				for(int i=0; i<abc.length;i++) {
					content += abc[i]+",";					
				}
				content = content.substring(0, content.length()-1);
				Util.write(args[0]+"new", content+"\n");
			}
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
}
