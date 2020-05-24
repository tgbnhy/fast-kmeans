package au.edu.rmit.trajectory.clustering.validity;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import au.edu.rmit.trajectory.clustering.kpaths.ClusterPath;
import au.edu.rmit.trajectory.clustering.kpaths.Util;

/**
 * Created by Josem on 17/02/2017.
 */
public class KMeansCVIs {
	Map<Integer, int[]> datamap;
	Map<Integer, double[]> idPoints;
	int distanceModel = 0;// the distance model
	double maxDis[];
    private Indice silhouette;
    private Indice dunn;
    private Indice bdSilhouette;
    private Indice bdDunn;
    private Indice davidBouldin;
    private Indice calinskiHarabasz;
    private Indice maximumDiameter;
    private Indice squaredDistance;
    private Indice averageDistance;
    private Indice averageBetweenClusterDistance;
    private Indice minimumDistance;
    private Indice XieB;
    private List<ClusterPath> clusters;

    public KMeansCVIs(int k, int dis, String file, Map<Integer, int[]> datamap, Map<Integer, double[]> idPoints) {
        distanceModel = dis;
        this.clusters = new ArrayList<ClusterPath>(k);
        this.datamap = datamap;
        this.idPoints = idPoints;
        try {
			Scanner in = new Scanner(new BufferedReader(new FileReader(file)));			
			while (in.hasNextLine()) {// load the trajectory dataset, and we can efficiently find the trajectory by their id.
				String str = in.nextLine();
				String strr = str.trim();
				String[] abc = strr.split(":");
				String[] traids = abc[1].split(",");
				Set<Integer> candidates = new HashSet<>();
				for(String ids: traids) {
					candidates.add(Integer.valueOf(ids));
				}
				this.clusters.add(new ClusterPath(candidates, Integer.valueOf(abc[0])));
			}
			in.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
    }
    
    public KMeansCVIs(int k, int dis, ArrayList<ClusterPath> CENTERS, Map<Integer, int[]> datamap, Map<Integer, double[]> idPoints, double[] maxDis) {
        distanceModel = dis;
        this.clusters = new ArrayList<ClusterPath>(k);
        this.datamap = datamap;
        this.idPoints = idPoints;
        this.clusters = CENTERS;
        this.maxDis = maxDis;
    //    System.out.println(CENTERS.get(0).getClusterTrajectories());
    }

    public String calculaIndices() {
        this.silhouette = calcularSilhouette();
        this.dunn = calcularDunn();
        this.bdSilhouette = calcularBDSilhouette();
        this.bdDunn = calcularBDDunn();
        this.davidBouldin = calcularDavidBouldin();
        this.calinskiHarabasz = calcularCalinskiHarabasz();
        this.maximumDiameter = calcularMaximumDiameter();
        this.squaredDistance = calcularSquaredDistance();
        this.averageDistance = calcularAverageDistance();
        this.averageBetweenClusterDistance = calcularAverageBetweenClusterDistance();
        this.minimumDistance = calcularMinimumDistance();
        this.XieB = calcularXieBeni();
        System.err.print(distanceModel + "_"+clusters.size()+"_"+datamap.size() + ":\t");
        String content = silhouette.getResultado()+" "+dunn.getResultado()+" "+bdSilhouette.getResultado()+" "
        + bdDunn.getResultado()+" "+davidBouldin.getResultado()+" "+
        		calinskiHarabasz.getResultado()+" "+maximumDiameter.getResultado()+" "+squaredDistance.getResultado()
        		+" "+ averageDistance.getResultado()+" "
        		+averageBetweenClusterDistance.getResultado()+" "+minimumDistance.getResultado()+" "+XieB.getResultado();
        System.err.printf("%.2f", silhouette.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", dunn.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", bdSilhouette.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", bdDunn.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", davidBouldin.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", calinskiHarabasz.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", maximumDiameter.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", squaredDistance.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", averageDistance.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", averageBetweenClusterDistance.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", minimumDistance.getResultado());
        System.err.print("\t");
        System.err.printf("%.2f", XieB.getResultado());
        System.err.println();
        return content;
    }
    
    /*
     * we normalize the distance 
     */
    private double distance(int punto, int punto2) { 
    	int []tra = datamap.get(punto);
    	int []clustra = datamap.get(punto2);
    	double min_dist = 0;
		if (distanceModel == 0)
			min_dist = Util.Intersection(tra, clustra)/maxDis[0];
		else if( distanceModel == 1) {
			min_dist = Util.EDR(tra, clustra)/maxDis[1];
		}else if (distanceModel == 2){
			min_dist = Util.DTW(tra, clustra, idPoints)/maxDis[2];
		}else if (distanceModel == 3){
			min_dist = Util.Frechet(tra, clustra, idPoints)/maxDis[3];
		}else if (distanceModel == 4){
			min_dist = Util.Hausdorff(tra, clustra, idPoints)/maxDis[4];
		}else if (distanceModel == 5){
			min_dist = Util.ERP(tra, clustra, idPoints)/maxDis[5];
		}
    	return min_dist;
	}
    
    /*
     * we cannot use one bencemark to measure all, then sum up all the distance
     */
    private double distance1(int punto, int punto2) {    
    	int []tra = datamap.get(punto);
    	int []clustra = datamap.get(punto2);
    	double min_dist;
		min_dist = Util.Intersection(tra, clustra)/maxDis[0]+
				Util.EDR(tra, clustra)/maxDis[1]+
				Util.DTW(tra, clustra, idPoints)/maxDis[2]+
				Util.Frechet(tra, clustra, idPoints)/maxDis[3]+
				Util.Hausdorff(tra, clustra, idPoints)/maxDis[4]+
				Util.ERP(tra, clustra, idPoints)/maxDis[5];
    	return min_dist;
	}
    
    private Indice calcularDunn() {
        double dunn = 0.0;
        double max = 0;
        double min = Double.MAX_VALUE;
        long startTime = System.currentTimeMillis();
        try {
            for (ClusterPath cluster : clusters) {
                for (int punto : cluster.getClusterTrajectories()) {
                    for (ClusterPath cluster2 : clusters) {
                        if (!cluster.equals(cluster2)) {
                            for (int punto2 : cluster.getClusterTrajectories()) {
                                if (punto != punto2) {
                                    double dist = distance(punto, punto2);
                                    if (min >= 0) {
                                        if (dist < min) {
                                            min = dist;
                                        }
                                    } else {
                                        min = dist;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (ClusterPath cluster : clusters) {
                for (int punto : cluster.getClusterTrajectories()) {
                    for (int punto2 : cluster.getClusterTrajectories()) {
                        if (punto !=punto2) {
                            double dist = distance(punto, punto2);
                            if (dist > max) {
                                max = dist;
                            }
                        }

                    }
                }
            }
            dunn = min / max;
        } catch (Exception e) {
            e.printStackTrace();
        }
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(dunn, elapsedTime);

    }



	private Indice calcularBDDunn() {
        double bdDunn = 0.0;
        double max = 0;
        double min = Double.MAX_VALUE;

        long startTime = System.currentTimeMillis();
        try {

            for (ClusterPath cluster : clusters) {
                for (ClusterPath cluster2 : clusters) {
                    if (!cluster.equals(cluster2) && cluster.getTrajectoryID() >= 0 && cluster2.getTrajectoryID() >= 0) {
                        double dist = distance(cluster.getTrajectoryID(), cluster2.getTrajectoryID());
                        if (min != 0) {
                            if (dist < min) {
                                min = dist;
                            }
                        } else {
                            min = dist;
                        }
                    }
                }
            }


            //get the maximum distance of the points to the centroid of the cluster they belong to
            for (ClusterPath cluster : clusters) {
                if (cluster.getTrajectoryID() >= 0) {
                    for (int punto : cluster.getClusterTrajectories()) {
                        double dist = distance(punto, cluster.getTrajectoryID());
                        if (dist > max) {
                            max = dist;
                        }
                    }
                }
            }
            //System.out.println("MINIMO: " + min);
            //System.out.println("MAXIMO: " + max);

            bdDunn = min / max;
        } catch (Exception e) {
            e.printStackTrace();
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(bdDunn, elapsedTime);

    }
	private Indice calcularXieBeni() {
        double bdDunn = 0.0;
        double max = 0;
        double min = Double.MAX_VALUE;

        long startTime = System.currentTimeMillis();
        double firstPart = 0;
        double secondPart = 0;
        try {
            for (ClusterPath cluster : clusters) {
            	
                for (int ids : cluster.getClusterTrajectories()) {
                      double dist = distance(cluster.getTrajectoryID(), ids);
                      firstPart += dist*dist;
                }
            }


            //get the maximum distance of the points to the centroid of the cluster they belong to
            for (ClusterPath cluster : clusters) {
                if (cluster.getTrajectoryID() >= 0) {
                	 for (ClusterPath cluster1 : clusters) {
                         if (cluster1.getTrajectoryID() >= 0) {
                             if(cluster.getTrajectoryID()!=cluster1.getTrajectoryID()) {
                            	 double dist = distance(cluster.getTrajectoryID(), cluster1.getTrajectoryID());
                            	 secondPart += dist*dist;
                             }
                         }
                     }
                }
            }
            secondPart *= datamap.size();
            //System.out.println("MINIMO: " + min);
            //System.out.println("MAXIMO: " + max);

            bdDunn = firstPart / secondPart;
        } catch (Exception e) {
            e.printStackTrace();
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(bdDunn, elapsedTime);

    }

    private Indice calcularSilhouette() {
        double silhouette;
        double a;
        double distA = 0;
        double b;
        double distB = 0;
        double cont;
        long startTime = System.currentTimeMillis();
        for (ClusterPath cluster : clusters) {
            for (int punto : cluster.getClusterTrajectories()) {
                for (ClusterPath cluster2 : clusters) {
                    if (!cluster.equals(cluster2)) {
                        for (int punto2 : cluster.getClusterTrajectories()) {
                            if (punto!=punto2) {
                                distB += distance(punto, punto2);
                            }
                        }
                    }
                }
            }
        }
        b = distB / clusters.size();

        cont = 0;
        for (ClusterPath cluster : clusters) {
            for (int punto : cluster.getClusterTrajectories()) {
                for (int punto2 : cluster.getClusterTrajectories()) {
                    if (punto!=punto2) {
                        distA += distance(punto, punto2);
                        cont++;
                    }
                }
            }
        }
        a = distA / clusters.size();
        //System.out.println("A: " + a);
        //System.out.println("B: " + b);

        silhouette = (b - a) / Math.max(a, b);
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(silhouette, elapsedTime);
    }

    private Indice calcularBDSilhouette() {
        double silhouette;
        double a;
        double distA = 0;
        double b;
        double distB = 0;
        double cont = 0;
        double score=0;
        long startTime = System.currentTimeMillis();

        for (ClusterPath cluster : clusters) {
        	for(int tra: cluster.getClusterTrajectories()) {
        		 double distnwe = 0;
        		 double minvalue = Double.MAX_VALUE;
        		 for (ClusterPath cluster1: clusters) {
        			 if(cluster!=cluster1) {
        				 for(int tra1: cluster1.getClusterTrajectories()) {
        					 distnwe += distance(tra, tra1);
        				 }
        			 }
        			 distnwe /= cluster1.getClusterTrajectories().size();
        		 }
        		 if(minvalue>distnwe)
        			 minvalue = distnwe;
        		 double a1 = 0;
        		 for(int tra1: cluster.getClusterTrajectories()) {
        			 distnwe += distance(tra, tra1);
        		 }
        		 a1/= cluster.getClusterTrajectories().size();
        		 score += (distnwe-a1)/Math.max(distnwe, a1);
        	}
            if (cluster.getTrajectoryID() >= 0) {
                for (ClusterPath cluster2 : clusters) {
                    if (cluster2.getTrajectoryID() >= 0) {
                        if (!cluster.equals(cluster2)) {
                            distB += distance(cluster.getTrajectoryID(), cluster2.getTrajectoryID());
                            cont++;
                        }
                    }
                }
            }
        }

        b = distB / cont;

        cont = 0;
        for (ClusterPath cluster : clusters) {
            if (cluster.getTrajectoryID() >= 0) {
                for (int punto : cluster.getClusterTrajectories()) {
                    distA += distance(punto, cluster.getTrajectoryID());
                    cont++;
                }
            }
        }
        a = distA / cont;
        //System.out.println("A: " + a);
        //System.out.println("B: " + b);

        silhouette = (b - a) / Math.max(a, b);
        
        score /= datamap.size();
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(score, elapsedTime);
    }

    private Indice calcularDavidBouldin() {
        int numberOfClusters = clusters.size();
        double david = 0.0;

        long startTime = System.currentTimeMillis();

        if (numberOfClusters == 1) {
            throw new RuntimeException(
                    "Impossible to evaluate Davies-Bouldin index over a single cluster");
        } else {
            // counting distances within
            double[] withinClusterDistance = new double[numberOfClusters];

            int i = 0;
            for (ClusterPath cluster : clusters) {
                for (int punto : cluster.getClusterTrajectories()) {
                    withinClusterDistance[i] += distance(punto, cluster.getTrajectoryID());
                }
                withinClusterDistance[i] /= cluster.getClusterTrajectories().size();
                i++;
            }


            double result = 0.0;
            double max = Double.NEGATIVE_INFINITY;

            try {
                for (i = 0; i < numberOfClusters; i++) {
                    //if the cluster is null
                    if (clusters.get(i).getTrajectoryID() >= 0) {

                        for (int j = 0; j < numberOfClusters; j++)
                            //if the cluster is null
                            if (i != j && clusters.get(j).getTrajectoryID() >= 0) {
                                double val = (withinClusterDistance[i] + withinClusterDistance[j])
                                        / distance(clusters.get(i).getTrajectoryID(), clusters.get(j).getTrajectoryID());
                                if (val > max)
                                    max = val;
                            }
                    }
                    result = result + max;
                }
            } catch (Exception e) {
                System.out.println("Excepcion al calcular DAVID BOULDIN");
                e.printStackTrace();
            }
            david = result / numberOfClusters;
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(david, elapsedTime);
    }

    private Indice calcularCalinskiHarabasz() {
        double calinski = 0.0;
        double squaredInterCluter = 0;
        double aux;
        double cont = 0;

        long startTime = System.currentTimeMillis();

        try {
            for (ClusterPath cluster : clusters) {
                if (cluster.getTrajectoryID() != 0) {
                    for (ClusterPath cluster2 : clusters) {
                        if (cluster2.getTrajectoryID() >= 0) {
                            if (!cluster.equals(cluster2)) {
                                aux = distance(cluster.getTrajectoryID(), cluster2.getTrajectoryID());
                                squaredInterCluter += aux * aux;
                                cont++;
                            }
                        }
                    }
                }
            }

            calinski = (this.calcularSquaredDistance().getResultado()) / (squaredInterCluter / cont);
        } catch (Exception e) {
            System.out.println("Excepcion al calcular CALINSKI");
            e.printStackTrace();
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(calinski, elapsedTime);
    }

    private Indice calcularMaximumDiameter() {
        double maximumDiameter = 0;
        double aux;

        long startTime = System.currentTimeMillis();

        for (ClusterPath cluster : clusters) {
            for (int punto : cluster.getClusterTrajectories()) {
                for (int punto2 : cluster.getClusterTrajectories()) {
                    if (punto!=punto2) {
                        aux = distance(punto, punto2);
                        if (aux > maximumDiameter) {
                            maximumDiameter = aux;
                        }
                    }
                }
            }
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(maximumDiameter, elapsedTime);
    }

    private Indice calcularSquaredDistance() {
        double squaredDistance = 0;
        double aux;
        double cont = 0;

        long startTime = System.currentTimeMillis();

        for (ClusterPath cluster : clusters) {
            for (int punto : cluster.getClusterTrajectories()) {
                for (int punto2 : cluster.getClusterTrajectories()) {
                    if (punto!=punto2) {
                        aux = distance(punto, punto2);
                        squaredDistance += aux * aux;
                        cont++;
                    }
                }
            }
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(squaredDistance / cont, elapsedTime);
    }

    private Indice calcularAverageDistance() {
        double averageDistance;
        double distA = 0;
        double cont = 0;

        long startTime = System.currentTimeMillis();

        for (ClusterPath cluster : clusters) {
            for (int punto : cluster.getClusterTrajectories()) {
                for (int punto2 : cluster.getClusterTrajectories()) {
                    if (punto!=punto2) {
                        distA += distance(punto, punto2);
                        cont++;
                    }
                }
            }
        }
        averageDistance = distA / cont;


        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(averageDistance, elapsedTime);
    }

    private Indice calcularAverageBetweenClusterDistance() {
        double averageDistanceBetween;
        double distA = 0;
        double cont = 0;

        long startTime = System.currentTimeMillis();

        for (ClusterPath cluster : clusters) {
            for (int punto : cluster.getClusterTrajectories()) {
                for (ClusterPath cluster2 : clusters) {
                    if (!cluster.equals(cluster2)) {
                        for (int punto2 : cluster.getClusterTrajectories()) {
                            if (punto!=punto2) {
                                distA += distance(punto, punto2);
                                cont++;
                            }
                        }
                    }
                }
            }
        }
        averageDistanceBetween = distA / cont;

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(averageDistanceBetween, elapsedTime);

    }

    private Indice calcularMinimumDistance() {
        double minimumDistance = -1;
        double aux;

        long startTime = System.currentTimeMillis();

        for (ClusterPath cluster : clusters) {
            for (int punto : cluster.getClusterTrajectories()) {
                for (ClusterPath cluster2 : clusters) {
                    if (!cluster.equals(cluster2)) {
                        for (int punto2 : cluster.getClusterTrajectories()) {
                            if (punto!=punto2) {
                                if (minimumDistance == -1) {
                                    minimumDistance = distance(punto, punto2);
                                } else {
                                    aux = distance(punto, punto2);
                                    if (aux < minimumDistance)
                                        minimumDistance = aux;
                                }
                            }
                        }
                    }
                }
            }
        }
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        return new Indice(minimumDistance, elapsedTime);
    }


    public Indice getMaximumDiameter() {
        return maximumDiameter;
    }

    public void setMaximumDiameter(Indice maximumDiameter) {
        this.maximumDiameter = maximumDiameter;
    }

    public Indice getSilhouette() {
        return silhouette;
    }

    public void setSilhouette(Indice silhouette) {
        this.silhouette = silhouette;
    }

    public List<ClusterPath> getClusters() {
        return clusters;
    }

    public void setClusters(List<ClusterPath> clusters) {
        this.clusters = clusters;
    }

    public Indice getDunn() {
        return dunn;
    }

    public void setDunn(Indice dunn) {
        this.dunn = dunn;
    }

    public Indice getDavidBouldin() {
        return davidBouldin;
    }

    public void setDavidBouldin(Indice davidBouldin) {
        this.davidBouldin = davidBouldin;
    }

    public Indice getCalinskiHarabasz() {
        return calinskiHarabasz;
    }

    public void setCalinskiHarabasz(Indice calinskiHarabasz) {
        this.calinskiHarabasz = calinskiHarabasz;
    }

    public Indice getSquaredDistance() {
        return squaredDistance;
    }

    public void setSquaredDistance(Indice squaredDistance) {
        this.squaredDistance = squaredDistance;
    }

    public Indice getAverageDistance() {
        return averageDistance;
    }

    public void setAverageDistance(Indice averageDistance) {
        this.averageDistance = averageDistance;
    }

    public Indice getAverageBetweenClusterDistance() {
        return averageBetweenClusterDistance;
    }

    public void setAverageBetweenClusterDistance(Indice averageBetweenClusterDistance) {
        this.averageBetweenClusterDistance = averageBetweenClusterDistance;
    }

    public Indice getMinimumDistance() {
        return minimumDistance;
    }

    public void setMinimumDistance(Indice minimumDistance) {
        this.minimumDistance = minimumDistance;
    }

    public Indice getBdSilhouette() {
        return bdSilhouette;
    }

    public void setBdSilhouette(Indice bdSilhouette) {
        this.bdSilhouette = bdSilhouette;
    }

    public Indice getBdDunn() {
        return bdDunn;
    }

    public void setBdDunn(Indice bdDunn) {
        this.bdDunn = bdDunn;
    }
}



