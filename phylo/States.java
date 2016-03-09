
package yeswecan.phylo;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 */
public class States {
    
    public static final int NT_STATES = 4;
    
    //public static final char[] BASES = { 'T', 'C', 'A', 'G' };
    
    public static final char[] BASES = { 'A', 'C', 'G', 'T' }; 
    //determines the order in which bases will be organised. 
    //Note that PAL uses ACGT and therefore deviating from this could cause issues in the pruning algorithm

    /*
        Number of states here is number of states in total
        ie for nucleotide there are 4 states total, 64 for codons etc
        N = 4 means there are N-1 non-originalStates, so an array or length N-1 is returned
    */
    public static int[] getMutationStates(int originalState, int numberOfStates){
        int[] mutationStates = new int[numberOfStates-1];
        
        int mutState = 0;
        for (int jIndex = 0; jIndex < mutationStates.length; jIndex++) {
            
            if (mutState == originalState){ // move onto next possible mutation state
                mutState++;
            }
            
            mutationStates[jIndex] = mutState;            
            mutState++;
      
        }
        
        return mutationStates;
    }

    // draw a state (0/1/2/3 for nt states) from a probability distribution
    public static int draw(double[] distribution, double randomUniform) {
        // NB elements of distribution MUST sum to one!
        // randomUniform must be between 0 and 1
        //ArrayPrinter.print(distribution, ",");
        double current = 0.0;
        for (int i = 0; i < distribution.length; i++) {
            if (randomUniform > current && randomUniform < current + distribution[i]) {
                return i;
            } else {
                current += distribution[i];
            }
        }
        return distribution.length - 1;
    }
    
}//class
