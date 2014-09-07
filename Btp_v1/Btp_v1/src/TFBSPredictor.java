import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import pmatch.*;


public class TFBSPredictor {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 */
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		// TODO Auto-generated method stub
		//To get the program running time ... to compare for different d-score values
		long startTime = System.currentTimeMillis();	
		
		String geneName = "TEMP";
		String fileLocation = "/home/harshit/Desktop/.course/7thsemester/btp/BTP/GENES/";
        String outputLocation = "/home/harshit/Desktop/.course/7thsemester/btp/BTP/output1/";
        String matrixLocation = "/home/harshit/Desktop/.course/7thsemester/btp/BTP/matrix.dat";
        
		BufferedReader in = new BufferedReader(new FileReader(fileLocation + geneName.toUpperCase()));
		String sequence = in.readLine();
		in.close();
		double dScoreCutOffMatrix = 1;
		double dScoreCutOffCore = 1;
		int conservedSitesNumber = 5;
		
		PositionWeightMatrix.setSavedLocPath(outputLocation + "pwm.dat");
		//to generate PWMs, one-time only , uncomment below two lines if running for the first time and comment out for the next rounds for faster execution
//		PositionWeightMatrix pwm = new PositionWeightMatrix(new File(matrixLocation));
//		pwm.PWMCreator();
		
		TFBSListPerSequence t = new TFBSListPerSequence(sequence,dScoreCutOffMatrix, dScoreCutOffCore, conservedSitesNumber);
		t.predictTFBS();
		System.out.println("Predictions completed");
		ArrayList<TFBindingSites> tfbsList = t.getTFBSList();
		
		FileWriter fstream = new FileWriter(outputLocation + geneName.toUpperCase() + "_TFBSv2.tsv");
        BufferedWriter out = new BufferedWriter(fstream);
        
		for(TFBindingSites tfbs: tfbsList){
			for(Integer location: tfbs.getLocationList()){
				out.write(tfbs.getName() + "\t" + tfbs.getMatchCount() + "\t" + location + "\n");
			}
		}
		long endTime = System.currentTimeMillis();
		long programDuration = (endTime - startTime)/1000;
		System.out.println("The program ran for: "+programDuration + "seconds");
		out.close();
	}
}
