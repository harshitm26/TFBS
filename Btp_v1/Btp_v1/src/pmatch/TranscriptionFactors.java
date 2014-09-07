package pmatch;

import java.io.Serializable;

public class TranscriptionFactors implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String name;
    private double[][] weightMatrix;

    //constructor
    public TranscriptionFactors(String name, double[][] weightMatrix) {
        this.name = name;
        this.weightMatrix = weightMatrix;
    }

    //get name of TF
    public String getName() {
        return this.name;
    }

    //get Position Weight Matrix of TF
    public double[][] getWtMatrix() {
        return this.weightMatrix;
    }
    
    //get number of rows in frequency matrix of a TF
    public int length() {
        return this.weightMatrix.length;
    }
}
