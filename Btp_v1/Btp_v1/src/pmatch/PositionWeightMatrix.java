package pmatch;

import java.io.*;
import java.util.*;

public class PositionWeightMatrix{
	private File freqMatrix;                          //matrix.dat
    private ArrayList<TranscriptionFactors> tfList;   //List of TFs with PWM
    private static String savedLocPath;

    //constructor
    public PositionWeightMatrix(File f) {
    	freqMatrix = f;
        tfList = new ArrayList<TranscriptionFactors>();
    }
    
    //set save location path
    public static void setSavedLocPath(String savedLocPath){
    	PositionWeightMatrix.savedLocPath = savedLocPath;
    }
    
    //get save location path
    public static String getSavedLocPath(String savedLocPath){
    	return savedLocPath;
    }
    
   //save TFs in the input file
    public static void save(ArrayList<TranscriptionFactors> tfList) throws IOException {
        ObjectOutputStream savedStream = new ObjectOutputStream(new FileOutputStream(savedLocPath));
        savedStream.writeObject(tfList);
        savedStream.close();
    }
    
   //retrieve TFs for predicting binding sites
    public static ArrayList<TranscriptionFactors> retrieve() throws IOException, ClassNotFoundException {
        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(savedLocPath));
        ArrayList<TranscriptionFactors> tfList = (ArrayList<TranscriptionFactors>) ois.readObject();
        return tfList;
    }
    
  //parse frequency matrix row from string to array
    public static double[] parseFrequency(String freqString){
    	if(freqString == null){
    		return null;
    	}
    	double[] frequencyArray = new double[4];
    	StringTokenizer st = new StringTokenizer(freqString);
    	double sum=0;
    	for(int i=0;i<4;i++){
    		frequencyArray[i] = Double.valueOf(st.nextToken()).doubleValue();
    		sum += frequencyArray[i];
    	}
    	for(int i=0;i<4;i++){
    		frequencyArray[i] /= sum;
    	}    	
    	return frequencyArray;
        
    }

    //checks if string is numeric type
    public static boolean isNumber(String s){
    	if(s == null){
    		return false;
    	}
    	try{
    		Integer.parseInt(s);
    	}catch(NumberFormatException nfe){  
    		return false;
    	}
    	return true;
    }
    
    //computes PWM and adds to the list of TFs
    private void fillPWM(String tfName, String frequencyMap){
    	if(tfName == null || frequencyMap == null){
    		return;
    	}
    	StringTokenizer st = new StringTokenizer(frequencyMap, "\n");
        double[][] weightMatrix = new double[st.countTokens()][4];
        double[] frequencyArray;
        double maximalInfo, netNucleotideEntropy;
        netNucleotideEntropy =  Math.log(0.25)/Math.log(2);
        String tmpBuffer;
        for(int i=0; i<weightMatrix.length; i++) {
            tmpBuffer = st.nextToken();
            //get a row of frequency matrix
            frequencyArray = parseFrequency(tmpBuffer);
            maximalInfo = 0;
            //maximalInfo estimation
            if(frequencyArray != null){
                for(int j=0; j<4; j++) {
                    if(frequencyArray[j]!=0){
                    	maximalInfo += frequencyArray[j]*Math.log(frequencyArray[j])/Math.log(2);
                    }
                }
                maximalInfo -= netNucleotideEntropy;
            }
            //generating Position Weight Matrix
            for(int j=0; j<4; j++) {
                weightMatrix[i][j] = maximalInfo*frequencyArray[j];
            }
        }
        //add TF with PWM to the list
        this.tfList.add(new TranscriptionFactors(tfName, weightMatrix));
    }
    
    //creates Position Weight Matrix for freqMatrix("matrix.dat"), stores result as list of TranscriptionFactors
    public void PWMCreator() throws IOException{
    	BufferedReader br = new BufferedReader(new FileReader(freqMatrix));
    	String matrixReadLine, tfName, frequencyMap;
    	tfName = null;
    	frequencyMap = "";
    	while((matrixReadLine = br.readLine()) != null){
    		//beginning of a TF in "matrix.dat"
    		if(matrixReadLine.equals("XX")){
    			//first TF in the file
    			if(tfName == null || frequencyMap.equals("")){
    				continue;
    			}
    			//fills the weight matrix of TF by utilizing frequencyMap
    			fillPWM(tfName,frequencyMap);
    			//initialize
    			tfName = null;
    			frequencyMap = "";
    		}
    		else if (matrixReadLine.length() >= 3){
    			//check for name of TF
    			if (matrixReadLine.substring(0, 2).equals("NA")){
    				tfName = matrixReadLine.substring(3, matrixReadLine.length()).trim();
    			}
    			//check for frequency row
    			else if(isNumber(matrixReadLine.substring(0, 2))){
    				frequencyMap = frequencyMap + matrixReadLine.substring(3, matrixReadLine.length()) + "\n";
    			}	
    		}	
    	}
    	//save in object stream
    	save(tfList);
    }
}
