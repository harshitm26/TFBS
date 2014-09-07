package pmatch;

import java.util.*;

public class TFBindingSites {
	private String name;                       //TF name
	private ArrayList<Integer> locationList;   //Starting position of TF in a sequence
	private int matchCount;
	
	//constructor
	public TFBindingSites(String name){
		this.name = name;
		this.locationList = new ArrayList<Integer>();
		this.matchCount = 0;
	}
	
	public String getName(){
		return name;
	}
	
	public ArrayList<Integer> getLocationList(){
		return locationList;
	}
	
	public void addLocation(int location){
		this.locationList.add(location);
	}
	
	public int getMatchCount(){
		return matchCount;
	}
	
	public void incrementMatchCount(){
		matchCount++;
	}
}
