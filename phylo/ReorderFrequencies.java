/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.


    PAML orders nucleotide frequencies as TCAG
    PAL orders them ACGT
    This class has methods for mapping between the two.
 */

package yeswecan.phylo;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class ReorderFrequencies {
    
    public static double[] pamlToAlpha(double[] tcag){
        double[] acgt = new double[4];
        acgt[0] = tcag[2];
        acgt[1] = tcag[1];
        acgt[2] = tcag[3];
        acgt[3] = tcag[0];
        return acgt;
    }
    
    public static double[] alphaToPaml(double[] acgt){
        double[] tcag = new double[4];
        tcag[2] = acgt[0];
        tcag[1] = acgt[1];
        tcag[3] = acgt[2];
        tcag[0] = acgt[3];
        return tcag;
    }
    
    

    // TODO could use these methods in the above
    private static int pamlToAlpha(int tcag){
        switch (tcag){
            case 0: // T
                return 3;
            case 1: // C
                return 1;
            case 2: // A
                return 0;
            case 3: // G
                return 2;
            default:
                return -1;
        }
    }
    
    
    
    private static int alphaToPaml(int acgt){
        switch (acgt){
            case 0: // A
                return 2;
            case 1: // C
                return 1;
            case 2: // G
                return 3;
            case 3: // T
                return 0;
            default:
                return -1;
        }
    
    }
    
    public static int[] pamlToAlpha(int[] tcagArray){
        int[] acgtArray = new int[tcagArray.length];
        for (int i = 0; i < acgtArray.length; i++) {
            acgtArray[i] = ReorderFrequencies.pamlToAlpha(tcagArray[i]);
        }
        return acgtArray;
    }
    
    public static int[] alphaToPaml(int[] acgtArray){
        int[] tcagArray = new int[acgtArray.length];
        for (int i = 0; i <tcagArray.length; i++) {
            tcagArray[i] = ReorderFrequencies.alphaToPaml(acgtArray[i]);
        }
        return tcagArray;
    }
}
