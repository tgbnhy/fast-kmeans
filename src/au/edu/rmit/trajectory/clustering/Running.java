package au.edu.rmit.trajectory.clustering;

import java.io.IOException;
import java.sql.SQLException;

import au.edu.rmit.mtree.MTree;
import au.edu.rmit.mtree.MTree.KMeansHMTree;
import au.edu.rmit.trajectory.clustering.kpaths.KPathsOptimization;
import au.edu.rmit.trajectory.clustering.kpaths.KPaths;

public class Running {
	public static void main(String[] args) throws IOException, SQLException, InterruptedException {
		boolean indexPivot = true;
		boolean indexInverte = true;
	//	Process aProcess = new Process(args);
	//	aProcess.staticKpath(args, false);
	//1565595, porto
		//250997, tdrive
		double time[] = new double[5];
		int counter = 0;
	//	for(int i=0; i<100; i++)//test 100 times
		for(int scale = 556559 ; scale<=556559; scale*=10) {//create the index for different using different radius.
			args[2] = Integer.toString(scale);
			KPathsOptimization run2 = new KPathsOptimization(args);
			run2.staticKpath(false, 0); // test routes without CPEP
		//	run2.staticKpath(true, i); // test routes with CPEP, it works well for T-drive.
		// count the time to build index
			long starttime = System.nanoTime();
		//	run2.runIndexbuildQueue(10, 20); // the first parameter is radius, the second is the capacity, 
			long endtime = System.nanoTime();
			time[counter++] = (endtime-starttime)/1000000000.0;
		}
		
		for(int i= 0; i<5; i++)
			System.out.println("\n\n"+time[i]+"\n\n");// count the time
	}
}