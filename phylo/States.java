
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
    
}//class
