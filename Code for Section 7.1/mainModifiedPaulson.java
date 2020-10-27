//This is the main file to run expertiments related to the PPP in Section 7.1 
package comparePaulsonandModifiedPaulson;
import java.util.ArrayList;
public class mainModifiedPaulson {
	public static void main(String[] args) {
		int numOfK = 10;
		double[][] recordResults = new double[numOfK][3];
		
		
		for(int i = 1 ; i <= numOfK; i++ ) {
			int k = i * 100;
			
			int repeatTime = 1000;
			double sumCPUTime = 0d;
			int sumCorrectness=0,sumSampleSize=0;
			for(int j = 0 ; j < repeatTime; j++) {
				//Input Parameters: k n_0 delta lambda alpha
				ModifiedPaulson y = new ModifiedPaulson(k,50,0.1,0.05,0.05);
				long startTime = System.currentTimeMillis();
				ArrayList<Integer> RESULT = new ArrayList<Integer>(y.getResult());

				long endTime = System.currentTimeMillis();
				sumCPUTime = sumCPUTime+endTime -  startTime;
				sumCorrectness = sumCorrectness+RESULT.get(0);
				sumSampleSize = sumSampleSize+RESULT.get(1);
			}
			recordResults[i-1][0] = (sumCPUTime/repeatTime);
			recordResults[i-1][1] = (sumCorrectness*1.0/repeatTime);
			recordResults[i-1][2] = (sumSampleSize*1.0/repeatTime);
			System.out.println(i+" "+recordResults[i-1][0]+" "+recordResults[i-1][1]+" "+recordResults[i-1][2]);
		}
	}
}
