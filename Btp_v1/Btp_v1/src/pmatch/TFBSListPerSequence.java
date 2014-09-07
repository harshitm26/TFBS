package pmatch;

import java.io.IOException;
import java.util.*;

public class TFBSListPerSequence {
	private String sequence;
	private ArrayList<TFBindingSites> TFBSList;
	private double dScoreCutOffMatrix;
	private double dScoreCutOffCore;
	private static int conservedSitesNumber;
	
	//constructor
	public TFBSListPerSequence(String sequence, double dScoreCutOffMatrix, double dScoreCutOffCore, int conservedSitesNumber){
		this. sequence = sequence;
		this.dScoreCutOffMatrix = dScoreCutOffMatrix;
		this.dScoreCutOffCore = dScoreCutOffCore;
		this.conservedSitesNumber = conservedSitesNumber;
		TFBSList = new ArrayList<TFBindingSites>();
	}
	
	//get list of TF Binding Sites
	public ArrayList<TFBindingSites> getTFBSList(){
		return this.TFBSList;
	}
	
	//computes D-Score
	public double dScoreCalculator(double maxWeight, double weightDiff){
		if(maxWeight == 0){
			return Double.POSITIVE_INFINITY;
		}
		return ((maxWeight - weightDiff)/maxWeight);
	}
	
	//returns TFBS from the list
	public TFBindingSites getTFBS(String name){
		if(name == null){
			return null;
		}
		for(TFBindingSites tfbs : TFBSList){
			if(tfbs.getName().equals(name)){
				return tfbs;
			}
		}
		return null;
	}
	
	//computes maximum weight for each row
	public static double[] maximums(double[][] weightMatrix) {
		if(weightMatrix == null){
			return null;
		}
        double max;
        double[] maximum = new double[weightMatrix.length];
        for(int i=0; i<weightMatrix.length; i++) {
            max = 0;
            for(int j=0; j<weightMatrix[i].length; j++) {
                if(weightMatrix[i][j] > max){
                	max = weightMatrix[i][j];
                }
            }
            maximum[i] = max;
        }
        return maximum;
    }
	
	//returns sum of maximum weights
	public static double getMaxWeight(double[] maximum) throws IOException{
		if(maximum == null){
			return Double.POSITIVE_INFINITY;
		}
		double sum = 0;
		for(int i=0; i<maximum.length;i++){
			sum += maximum[i];
		}
		return sum;
	}
	
	static class columnComparator implements Comparator {
	    private int columnToSortOn;
	   
	    //contructor to set the column to sort on.
	    columnComparator(int columnToSortOn) {
	      this.columnToSortOn = columnToSortOn;
	    }

	    // Implement the abstract method which tells
	    // how to order the two elements in the array.
	    public int compare(Object o1, Object o2) {
	    // cast the object args back to string arrays
	        double[] row1 = (double[])o1;
	        double[] row2 = (double[])o2;

	        // compare the desired column values
	        if(row1[columnToSortOn] < row2[columnToSortOn]){
	        	return -1;
	        }
	        else if(row1[columnToSortOn] == row2[columnToSortOn]){
	        	return 0;
	        }
	        else{
	        	return 1;
	        }
//	        return row1[columnToSortOn].compareTo(row2[columnToSortOn]);
	    }
	}
	
	
	public static double[][] findCoreSites(double[] maximum, int conservedSitesNumber){
		double[][] arrayIndices = new double[conservedSitesNumber][2];
		for(int i=0;i<maximum.length;i++){
			if(arrayIndices[0][0] < maximum[i]){
				arrayIndices[0][1] = i;
				arrayIndices[0][0] = maximum[i];
//				swap(0,minindex(maximum[i]));
				Arrays.sort(arrayIndices, new columnComparator(1));
			}
		}
		return arrayIndices;
	}
	
	//returns sum of weight differences over the entire length
	public static double getWeightDiff(double[][] weightMatrix, String s, double[] maximum) {
        double ret = 0;
        for(int i=0; i<s.length(); i++) {
            ret += Math.abs(maximum[i]-weightMatrix[i][getID(s.charAt(i))]);
        }
        return ret;
    }
	
	 //returns sum of weight difference over only 5 core sites( i.e. 5 locations where PWM score is maximum)
	public static double getWeightDiffCore(double[][] weightMatrix, String s, double[] maximum) {
		double ret = 0;
		double arrayIndices[][] = new double[conservedSitesNumber][2];
		arrayIndices = findCoreSites(maximum, conservedSitesNumber);
		for(int i=0; i<conservedSitesNumber; i++){
			int j = (int)arrayIndices[i][1];
			System.out.println(s.charAt(j) + " " + maximum[j]);
			ret += Math.abs(maximum[j] - weightMatrix[j][getID(s.charAt(j))]);
		}
		System.out.println("");
		return ret;
	}
	
	
	//get index by char
    public static int getID(char ch) {
        if(ch == 'a' || ch == 'A'){
        	return 0;
        }
        if(ch == 'c' || ch == 'C'){
        	return 1;
        }
        if(ch == 'g' || ch == 'G'){
        	return 2;
        }
        if(ch == 't' || ch == 'T'){
        	return 3;
        }
        System.out.println("Not among {a,c,g,t}: "+ch);
        return -1;
    }
	
	//predicting TFBS
	public void predictTFBS() throws ClassNotFoundException, IOException {
        int tfLength, foundSites;
        double maxWeight, dScore, dScoreCoreSites, weightDiff, weightDiffCoreSites;
        //get stored TFs with PWM
        ArrayList<TranscriptionFactors> tFactors = PositionWeightMatrix.retrieve();
        for (TranscriptionFactors tf : tFactors ) {
            tfLength = tf.length();
            double[][] weightMatrix = tf.getWtMatrix();
            //set the maximums of each row
            double[] maximum = maximums(weightMatrix);
            //sum of maximums
            maxWeight = getMaxWeight(maximum);
            //sites found per TF
            foundSites = 0;
            for(int j=0; j+tfLength <=sequence.length(); j++) {
            	weightDiff = getWeightDiff(tf.getWtMatrix(), sequence.substring(j,j+tfLength), maximum);
            	weightDiffCoreSites = getWeightDiffCore(tf.getWtMatrix(), sequence.substring(j,j+tfLength), maximum);
            	
                dScore = dScoreCalculator(maxWeight, weightDiff);
                dScoreCoreSites = dScoreCalculator(maxWeight, weightDiffCoreSites);
                
                if(dScore < 0 || dScore > 1){
                	System.out.println("INVALID D-Score: "+ dScore + " " + tf.getName());
                }
                if(dScoreCoreSites <0 || dScoreCoreSites>1){
                	System.out.println("INVALID D-Score Core Site:" + dScoreCoreSites +" "+ tf.getName());
                }
                
                if(dScore >= dScoreCutOffMatrix && dScoreCoreSites >= dScoreCutOffCore) {
                	//first TF in the sequence
                	if(foundSites == 0){
                    	TFBSList.add(new TFBindingSites(tf.getName()));
                	}
                	//retrieve TFBS from the list
                	TFBindingSites tfbs = getTFBS(tf.getName());
                	if(tfbs != null){
                		//add starting location
                		tfbs.addLocation(j);
                		//increment count
                		tfbs.incrementMatchCount();
                	}
                }
            }
        }
    }
}
