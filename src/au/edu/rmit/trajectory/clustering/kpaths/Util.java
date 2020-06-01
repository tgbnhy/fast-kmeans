package au.edu.rmit.trajectory.clustering.kpaths;

import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;

import org.netlib.util.doubleW;

import kafka.utils.Exit;
import scala.Int;
import com.google.*;
import com.google.common.math.Stats;

import java_cup.internal_error;

public class Util {

	public Util() {
		// TODO Auto-generated constructor stub
	}

	// write the information into files
	public static void write(String fileName, String content) {
		RandomAccessFile randomFile = null;
		try {
			randomFile = new RandomAccessFile(fileName, "rw");
			long fileLength = randomFile.length();
			randomFile.seek(fileLength);
			randomFile.writeBytes(content);
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (randomFile != null) {
				try {
					randomFile.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	// write the information into files
		public static void rewrite(String fileName, String content) {
			RandomAccessFile randomFile = null;
			try {
				randomFile = new RandomAccessFile(fileName, "rws");
				long fileLength = randomFile.length();
				randomFile.seek(fileLength);
				randomFile.writeBytes(content);
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (randomFile != null) {
					try {
						randomFile.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
		}
	
	/*
	 * we compute the distance on Euclidean space.
	 */
	public static double EuclideanDis(double []a, double []b, int dimension) {
		double distance = 0;
		for(int i = 0; i<dimension; i++ ) {
			distance += Math.pow(Math.abs(a[i]- b[i]), 2);// the l2-norm
		}

		distance = Math.sqrt(distance);

		Double x = new Double(distance);
		boolean isNan = x.isNaN();
		if(isNan) {
		/*	System.out.println(Arrays.toString(b));
			System.out.println(Arrays.toString(a));
			System.out.println(distance);*/
		}
		return distance;
	}
	/*
	 * the data needs to be sorted before the intersection, the edge-based distance (EBD)
	 */
	public static int Intersection(int arr1[], int arr2[], int m, int n) {
		int i = 0, j = 0;
		int dist = 0;
		while (i < m && j < n) {
			if (arr1[i] < arr2[j])
				i++;
			else if (arr2[j] < arr1[i])
				j++;
			else
			{
				dist++;
				i++;
				j++;
			}
		}
		return Math.max(m, n)-dist;
	}
	
	/*
	 * the data needs to be sorted before the intersection, the edge-based distance (EBD)
	 */
	public static int Intersection(int arr1[], int arr2[]) {
		int i = 0, j = 0;
		int dist = 0;
		while (i < arr1.length && j < arr2.length) {
			if (arr1[i] < arr2[j])
				i++;
			else if (arr2[j] < arr1[i])
				j++;
			else
			{
				dist++;
				i++;
				j++;
			}
		}
		return Math.max(arr1.length, arr2.length)-dist;
	}
	
	/*
	 * the data needs to be sorted before the intersection, the edge-based distance (EBD)
	 */
	public static int[] getIntersection(int arr1[], int arr2[], int k) {
		int centroidList[] = new int[k];
		for(int i=0; i<k; i++) {
			if(arr1[i]==1 && arr2[i]==1) {
				centroidList[i]=1;
			}else {
				centroidList[i]=0;
			}
		}
		return centroidList;
	}
	
	/*
	 * Edit distance on real sequence
	 */
	public static double EDR(int[] T1, int T2[]) {
		if (T1 == null || T1.length == 0) {
			if (T2 != null)
				return T2.length;
			else
				return 0;
		}
		if (T2 == null || T2.length == 0) {
			if (T1 != null)
				return T1.length;
			else
				return 0;
		}
		int[][] dpInts = new int[T1.length + 1][T2.length + 1];
		dpInts[0][0] = 0;
		for (int i = 1; i <= T1.length; ++i) {
			dpInts[i][0] = i;
		}
		for (int j = 1; j <= T2.length; ++j) {
			dpInts[0][j] = j;
		}
		for (int i = 1; i <= T1.length; ++i) {
			for (int j = 1; j <= T2.length; ++j) {
				int subCost = 1;
				if (T1[i - 1] == T2[j - 1])
					subCost = 0;
				dpInts[i][j] = min(dpInts[i - 1][j - 1] + subCost, dpInts[i - 1][j] + 1, dpInts[i][j - 1] + 1);
			}
		}
		return dpInts[T1.length][T2.length];
	}
	
	private static int min(int a, int b, int c) {
		if (a > b)
			a = b;
		if (a > c)
			a = c;
		return a;
	}
	
	/*
	 * the dynamic time warping
	 */
	public static double DTW(int[] T1, int[] T2, Map<Integer, double[]> idPoints) {
		if (T1.length == 0 && T2.length == 0)
			return 0;
		if (T1.length == 0 || T2.length == 0)
			return Integer.MAX_VALUE;
		double[][] dpInts = new double[T1.length + 1][T2.length + 1];
		for (int i = 1; i <= T1.length; ++i) {
			dpInts[i][0] = Integer.MAX_VALUE;
		}
		for (int j = 1; j <= T2.length; ++j) {
			dpInts[0][j] = Integer.MAX_VALUE;
		}
		for (int i = 1; i <= T1.length; ++i) {
			for (int j = 1; j <= T2.length; ++j) {
				dpInts[i][j] = EuclideanDis(T1[i - 1], T2[j - 1], idPoints)
						+ min(dpInts[i - 1][j - 1], dpInts[i - 1][j], dpInts[i][j - 1]);
			}
		}
		return dpInts[T1.length][T2.length];
	}
	
    private static double min(double a, double b, double c) {
        if (a > b) 
        	a = b;
        if (a > c) 
        	a = c;
        return a;
    }
    
    //The ERP distance function
    public static double ERP(int[] r, int[] s, Map<Integer, double[]> idPoints)
	{
		double[][] erpMetric = new double[r.length+ 1][s.length + 1];
		double rAB = 0;
		for (int i = 0; i < r.length; i++)
		{
			double []point1xy = idPoints.get(r[i]);
			rAB += Math.sqrt(Math.abs(point1xy[0]*point1xy[0]+ point1xy[1]*point1xy[1]));
		}
		double sAB = 0;
		for (int i = 0; i < s.length; i++){
			double []point1xy = idPoints.get(s[i]);
			sAB += Math.sqrt(Math.abs(point1xy[0]*point1xy[0]+ point1xy[1]*point1xy[1]));
		}
		for (int i = 0; i <= r.length; i++){
			erpMetric[i][0] = sAB;
		}
		for (int i = 0; i <= s.length; i++){
			erpMetric[0][i] = rAB;
		}
		erpMetric[0][0] = 0;
		for (int i = 1; i <= r.length; i++){
			for (int j = 1; j <= s.length; j++){
				double []point1xy = idPoints.get(r[i-1]);
				double []point2xy = idPoints.get(s[j-1]);
				erpMetric[i][j] = min(
				        erpMetric[i - 1][j - 1]
				                + EuclideanDis(r[i-1], s[j - 1], idPoints),
				        erpMetric[i - 1][j]
				                + Math.sqrt(point1xy[0]*point1xy[0] + point1xy[1]*point1xy[1]),
				        erpMetric[i][j - 1]
				                + Math.sqrt(point2xy[0]*point2xy[0] + point2xy[1]*point2xy[1]));
			}
		}
		return erpMetric[r.length][s.length];
	}
    
    /*
     * the hausdorff distance
     */
    public static double Hausdorff(int[] t1, int[] t2, Map<Integer, double[]> idPoints) {
        double[][] dist_matrix;
        dist_matrix = new double[t2.length][t1.length];
        double result = 0.0D;
        ArrayList<Double> minDistances1 = new ArrayList<Double>();
        ArrayList<Double> minDistances2 = new ArrayList<Double>();
        int i;
        for (i = 0; i < dist_matrix.length; ++i) {
            for (int j = 0; j < dist_matrix[0].length; ++j) {
                dist_matrix[i][j] = EuclideanDis(t1[j], t2[i], idPoints);
            }
        }

        int j;
        double min;
        for (i = 0; i < dist_matrix.length; ++i) {
            min = Double.MAX_VALUE;
            for (j = 0; j < dist_matrix[0].length; ++j) {
                if (dist_matrix[i][j] <= min) {
                    min = dist_matrix[i][j];
                }
            }
            minDistances1.add(min);
        }

        for (i = 0; i < dist_matrix[0].length; ++i) {
            min = Double.MAX_VALUE;
            for (j = 0; j < dist_matrix.length; ++j) {
                if (dist_matrix[j][i] <= min) {
                    min = dist_matrix[j][i];
                }
            }
            minDistances2.add(min);
        }
        Collections.sort(minDistances1);
        Collections.sort(minDistances2);
        double value1 = minDistances1.get(minDistances1.size() - 1);
        double value2 = minDistances2.get(minDistances2.size() - 1);
        result = Math.max(value1, value2);
        return result;
    }
    
    /*
     * the discrete frechet distance
     */
	public static double Frechet(int[] t1, int[] t2, Map<Integer, double[]> idPoints) {
		double[][] ca = new double[t2.length][t1.length];
		for (int i = 0; i < t2.length; ++i) {
			for (int j = 0; j < t1.length; ++j) {
				ca[i][j] = -1.0D;
			}
		}
		return c(t2.length - 1, t1.length - 1, ca, t1, t2, idPoints);
	}

	private static double c(int i, int j, double[][] ca, int[] t1, int[] t2, Map<Integer, double[]> idPoints) {
		if (ca[i][j] > -1.0D)
			return ca[i][j];
		if (i == 0 && j == 0) {
			ca[i][j] = EuclideanDis(t1[0], t2[0], idPoints);
		} else if (j == 0) {
			ca[i][j] = Math.max(c(i - 1, 0, ca, t1, t2, idPoints), EuclideanDis(t2[i], t1[0], idPoints));
		} else if (i == 0) {
			ca[i][j] = Math.max(c(0, j - 1, ca, t1, t2, idPoints), EuclideanDis(t2[0], t1[j], idPoints));
		} else {
			ca[i][j] = Math.max(Math.min(Math.min(c(i - 1, j, ca, t1, t2, idPoints), c(i - 1, j - 1, ca, t1, t2, idPoints)), 
							c(i, j - 1, ca, t1, t2, idPoints)),
					EuclideanDis(t2[i], t1[j], idPoints));
		}
		return ca[i][j];
	}
	
	/*
	 * compute the distance between two starting points of id.
	 */
	static double EuclideanDis(int point1, int point2, Map<Integer, double[]> idPoints) {
		double []point1xy = idPoints.get(point1);
		double []point2xy = idPoints.get(point2);
		return Math.sqrt(Math.pow((point1xy[0]-point2xy[0]), 2) + Math.pow((point1xy[1]-point2xy[1]), 2));
	}
	
	/*
	 * we use Guava stat to profile the dataset, and return values
	 */
	static double[] datasetProfiling(double[][] dataMatix) {
		
	//	Stats.meanOf(dataMatix);
		return null;
	}
	
	
	public static double calculateSD(ArrayList<Double> numArray, double mean)
    {
        double standardDeviation = 0.0;
        int length = numArray.size();
        for(double num: numArray) {
            standardDeviation += Math.pow(num - mean, 2);
        }

        return Math.sqrt(standardDeviation/length);
    }
	
	public static double calculateSDInt(ArrayList<Integer> numArray, double mean)
    {
        double standardDeviation = 0.0;
        int length = numArray.size();
        for(int num: numArray) {
            standardDeviation += Math.pow(num - mean, 2);
        }

        return Math.sqrt(standardDeviation/length);
    }
	
	public static double calculateMean(ArrayList<Double> numArray)
    {
        double sum = 0.0;
        int length = numArray.size();

        for(double num : numArray) {
            sum += num;
        }

        double mean = sum/length;
        return mean;
    }
}
