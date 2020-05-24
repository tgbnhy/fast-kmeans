package au.edu.rmit.trajectory.clustering.kpaths;
/*
 * this class is used for recording the running statistics
 */
public class RunLog {
	private double assignmentTime;
	private double refinementTime;
	private int iterationTimes;
	private double simiComputationTime;
	private double ioTime;
	private double historgramTime;
	private double accumulatedTime;
	
	public RunLog() {
		assignmentTime = 0;
		refinementTime = 0;
		iterationTimes = 0;
		simiComputationTime = 0;
		ioTime = 0;
		historgramTime =0;
		accumulatedTime =0;
	}
	
	public void clear() {
		assignmentTime = 0;
		refinementTime = 0;
		iterationTimes = 0;
		simiComputationTime = 0;
		ioTime = 0;
		historgramTime =0;
		accumulatedTime =0;
	}
	
	public void addAssignmentTime(double time) {
		assignmentTime += time;
	}
	
	public void addRefinementTime(double time) {
		refinementTime += time;
	}
	
	public void addsimiComputationTime(double time) {
		simiComputationTime += time;
	}
	
	public void addIOTime(double time) {
		ioTime += time;
	}
	
	public void addHistorgramTime(double time) {
		historgramTime += time;
	}
	
	public void addAccumulatedTime(double time) {
		accumulatedTime += time;
	}
	
	public void setIterationtimes(int ite) {
		iterationTimes = ite;
	}
	
	public double getAssignmentTime() {
		return assignmentTime;
	}
	
	public double getRefinementTime() {
		return refinementTime;
	}
	
	public double getsimiComputationTime() {
		return simiComputationTime;
	}
	
	public double getAlltime() {
		return assignmentTime + refinementTime;
	}

	public void printLog() {
		double alltime = assignmentTime + refinementTime;
		System.out.println("\n==============================================================");
		System.out.println("it used "+iterationTimes+" iterations to find the centroids");
		System.out.print("the overall running time: ");
		System.out.printf("%.3f", alltime);
		System.out.println("s");
		System.out.print("\tthe assignment time: ");
		System.out.printf("%.3f", (assignmentTime-historgramTime));
		System.out.println("s");
		System.out.print("\t\tthe IO and Index access time: ");
		System.out.printf("%.3f", ioTime);
		System.out.println("s");
		System.out.print("\t\tthe similarity computation time: ");
		System.out.printf("%.3f", simiComputationTime);
		System.out.println("s");
		System.out.print("\tthe refinement time: ");
		System.out.printf("%.3f", (refinementTime+historgramTime));
		System.out.println("s");
		System.out.print("\t\tthe edge histogram building time: ");
		System.out.printf("%.3f", historgramTime);
		System.out.println("s");
		System.out.print("\t\tthe accumulated histogram time: ");
		System.out.printf("%.3f", accumulatedTime);
		System.out.println("s");
		System.out.println("==============================================================\n");
	}
}
