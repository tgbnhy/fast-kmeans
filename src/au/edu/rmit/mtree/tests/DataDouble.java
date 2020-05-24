package au.edu.rmit.mtree.tests;

import au.edu.rmit.mtree.DistanceFunctions.EuclideanCoordinate;

public class DataDouble implements EuclideanCoordinate, Comparable<DataDouble> {
	
	private final double[] values;
	private final int hashCode;
	
	public DataDouble(double... values) {
		this.values = values;
		
		int hashCode = 1;
		for(double value : values) {
			hashCode = 31*hashCode + (int)value;
		}
		this.hashCode = hashCode;
	}
	
	public DataDouble(double[] values, int id) {
		this.values = values;
		this.hashCode = id;
	}
	
	@Override
	public int dimensions() {
		return values.length;
	}

	@Override
	public double get(int index) {
		return values[index];
	}
	
	
	@Override
	public int hashCode() {
		return hashCode;
	}
	
	@Override
	public boolean equals(Object obj) {
		if(obj instanceof DataDouble) {
			DataDouble that = (DataDouble) obj;
			if(this.dimensions() != that.dimensions()) {
				return false;
			}
			for(int i = 0; i < this.dimensions(); i++) {
				if(this.values[i] != that.values[i]) {
					return false;
				}
			}
			return true;
		} else {
			return false;
		}
	}
	
	@Override
	public int compareTo(DataDouble that) {
		int dimensions = Math.min(this.dimensions(), that.dimensions());
		for(int i = 0; i < dimensions; i++) {
			double v1 = this.values[i];
			double v2 = that.values[i];
			if(v1 > v2) {
				return +1;
			}
			if(v1 < v2) {
				return -1;
			}
		}
		
		if(this.dimensions() > dimensions) {
			return +1;
		}
		
		if(that.dimensions() > dimensions) {
			return -1;
		}
		
		return 0;
	}

	@Override
	public int getID() {
		return hashCode;
	}

	@Override
	public int[] getData() {
		return null;
	}
	
	@Override
	public double[] getData1() {
		return values;
	}
}