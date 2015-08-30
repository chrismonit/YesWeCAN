
package yeswecan.model.parameters;

import java.util.Hashtable;
import yeswecan.phylo.States;
import swmutsel.model.parameters.TsTvRatio;
import yeswecan.Constants;
/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 */
public class TsTvRatioAdvanced extends TsTvRatio {
    
    
    private Hashtable<Character, Character> baseTypes = new Hashtable();
    

    
    public TsTvRatioAdvanced(TsTvRatio kappaInstance){
        this(kappaInstance.get());
    }
    
    public TsTvRatioAdvanced(double kappa){
        super(kappa);
        
        // TODO faster to compare ints/bits than chars?
        // but still need an Integer instance for the hashtable rather than primitive types, might slow it down again?
        baseTypes.put('A', 'R'); //purine
        baseTypes.put('C', 'Y'); //pyrimidine
        baseTypes.put('G', 'R'); //purine
        baseTypes.put('T', 'Y'); //pyrimidine
    }//constrcutor
    
    
    public double getKappaIfTransition(int base_i, int base_j){

        if ( baseTypes.get( States.BASES[base_i] ) == baseTypes.get( States.BASES[base_j] ) ){
            return this.get();
        }
        else{
            return 1.0;
        }
    }
    
    public static void main(String[] args){
        System.out.println("Testing TransitionOrTransversion");
        TsTvRatioAdvanced t = new TsTvRatioAdvanced(2.0);
        int i = 0;
        int j = 2;
        double output = t.getKappaIfTransition(i, j);
        System.out.println("value = " + output);
    }//main for testing
    
}//class
