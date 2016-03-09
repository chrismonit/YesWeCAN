
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

    public static int[] getMutationStates(int originalState, int numberOfStates){
        int[] mutationStates = new int[numberOfStates];
        for (int j = 0; j < mutationStates.length; j++) {
            if (j == originalState){
                continue;
            }else{
                mutationStates[j] = j;
            }
        }
        return mutationStates;
    }
    
}//class
