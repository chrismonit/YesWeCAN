
package yeswecan.model;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * class for housing various means of calculating exp(Qt)
 */
public class TransitionProbabilityMatrix {
    
    RateMatrix Q;
    public TransitionProbabilityMatrix( RateMatrix Q ){ 
        this.Q = Q;
        
        
    }//constructor
    
    
    
}//class
