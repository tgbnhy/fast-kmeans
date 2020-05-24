package skyline0623.balltree;

public class DimensionNotMatchException extends Exception{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public String toString(){
		return "The two operators' dimension are not matched.";
	}
}
