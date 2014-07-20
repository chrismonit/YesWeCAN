
package yeswecan.model.parameters;

import java.util.Hashtable;
import yeswecan.phylo.States;
/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 */
public class TsTvRatio extends Parameter {
    
    double kappa;
    
    private Hashtable<Character, String> baseTypes = new Hashtable();
    
    public TsTvRatio(double kappa){
        this.kappa = kappa;
        
        baseTypes.put('A', "purine");
        baseTypes.put('C', "pyrimidine");
        baseTypes.put('G', "purine");
        baseTypes.put('T', "pyrimidine");
        
    }//constrcutor
    
    
    public double getKappaIfTransition(int base_i, int base_j){

        if ( baseTypes.get( States.BASES[base_i] ) == baseTypes.get( States.BASES[base_j] ) ){
            return kappa;
        }
        else{
            return 1.0;
        }
    }
    
    public static void main(String[] args){
        System.out.println("Testing TransitionOrTransversion");
        TsTvRatio t = new TsTvRatio(2.0);
        int i = 0;
        int j = 2;
        double output = t.getKappaIfTransition(i, j);
        System.out.println("value = " + output);
    }//main for testing
    
}//class
