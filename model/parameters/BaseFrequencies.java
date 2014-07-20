
package yeswecan.model.parameters;

import java.util.Hashtable;
import yeswecan.phylo.States;
import yeswecan.phylo.States;

/**
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 */
public class BaseFrequencies extends Parameter {
    
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
        this.frequencies = frequencies;
    }
    
    
    public double[] get(){ return frequencies; }
    
    
}//class
