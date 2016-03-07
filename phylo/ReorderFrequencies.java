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
    
    
    public static int[] pamlToAlpha(int[] tcag){
        int[] acgt = new int[4];
        acgt[0] = tcag[2];
        acgt[1] = tcag[1];
        acgt[2] = tcag[3];
        acgt[3] = tcag[0];
        return acgt;
    }
    

    public static int[] alphaToPaml(int[] acgt){
        int[] tcag = new int[4];
        tcag[2] = acgt[0];
        tcag[1] = acgt[1];
        tcag[3] = acgt[2];
        tcag[0] = acgt[3];
        return tcag;
    }
    
}
