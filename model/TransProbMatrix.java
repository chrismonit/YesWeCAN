
package yeswecan.model;
import org.apache.commons.math3.linear.EigenDecomposition;


/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * created 8/7/15. This class differs from TransitionProbabilityMatrix in that it
 * will make use of library code to perform eigendecompsotion, rather than trying
 * to do this myself
 */
public class TransProbMatrix {
    
    public TransProbMatrix(RateMatrix Q){
        // perform eigen decomposition
        
        //EigenDecomposition eigen = new EigenDecomposition(Q);
        
    }
    
    
    
    /*
    public static double[][] getP(double t){
    
        // return P(t)
        
    }
    */
    
    public static void main(String[] args){
        System.out.println("Testing TransProbMatrix");
        
        
        //TransProbMatrix P = new TransProbMatrix();
    }
    
}
