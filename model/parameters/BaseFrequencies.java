
package yeswecan.model.parameters;

import java.util.Hashtable;
import yeswecan.phylo.States;
import yeswecan.phylo.States;
import yeswecan.Constants;

/**
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 */
public class BaseFrequencies  {
    
//    private Hashtable<Character, Double> frequencies;
//    
//    public BaseFrequencies(double[] frequencyValues){ 
//        //assumes frequencyValues are in the same order as in States.BASES
//        for (int i = 0; i < States.NT_STATES; i++) {
//            frequencies.put( States.BASES[i], frequencyValues[i] );
//        }
//    }//constructor
    
    private double[] frequencies;
    
    public BaseFrequencies( double[] frequencies ){
        //check these sum to 1
        double sum = 0.0;
        for (int i = 0; i < frequencies.length; i++) {
            sum += frequencies[i];
        }
        if (sum < 1.0-Constants.EPSILON || sum > 1.0+Constants.EPSILON) {
            throw new RuntimeException("BaseFrequencies: Base frequencies do not sum to 1");
        }

        this.frequencies = frequencies;
    }
    
    
    public double[] get(){ return frequencies; }
    
    
}//class
