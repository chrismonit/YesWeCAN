/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model;

import pal.datatype.Codons;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.utils.ArrayPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CodonScaler implements RatioScaler {
    
    private CodonFrequencies codonFrequencies;
    private double[] probNonsyn;
    private int CODON_POSITIONS = 3;
    private int NUM_CODONS = 64;
    
    double[][] nonSynProportions;
    
    
    public CodonScaler(CodonFrequencies codonFrequencies){
        this.codonFrequencies = codonFrequencies;
        
        probNonsyn = new double[3];
        // compute the three pn values
            
        nonSynProportions = new double[][]{
            { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 1.0, 0.6666666666666666, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 1.0, 0.6666666666666666, 1.0, 0.6666666666666666, 1.0, 0.6666666666666666, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.6666666666666666, 1.0, 0.6666666666666666, 1.0  }, 
            { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  },
            { 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.0, 0.0, 0.0, 0.0, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.3333333333333333, 0.3333333333333333, 1.0, 0.3333333333333333, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666 }
        };
        
        for (int iCodonPos = 0; iCodonPos < 3; iCodonPos++) {
            
            double sum = 0.0;
            for (int iCodon = 0; iCodon < NUM_CODONS; iCodon++) { 
      
                int[] codon = Codons.getNucleotideStatesFromCodonIndex(iCodon); // codon bases ordered ACGT = 0123
                double codonFreq = this.codonFrequencies.getFrequency(ReorderFrequencies.alphaToPaml(codon)); // codonFrequencies class expects bases to be ordered TCAG = 0123
                //System.out.println("iCodon "+iCodon+" codon "+ArrayPrinter.toString(codon,",")+" codonFreq "+codonFreq+" proportion "+nonSynProportions[iCodonPos][iCodon]);
            
                sum += codonFreq * nonSynProportions[iCodonPos][iCodon];
            }
            probNonsyn[iCodonPos] = sum;
            System.out.println("sum "+sum);
        }
        
    }
    
    
    //TODO: should put in RatioScaler interface as static thing to be used by this and ProportionScaler classes
    private int[][] codonPositions = new int[][]{
            // a frame, b frame, c frame
            { 0, 2, 1 }, // alpha site
            { 1, 0, 2 }, // beta site
            { 2, 1, 0 }  // gamma site
    };
    
    
    
    
    public double get(double ratio, int siteType, int frame){
        
        int codonPos = codonPositions[siteType][frame];
        double probNonSyn = probNonsyn[codonPos];
        double probSyn = 1.0 - probNonSyn;
        return (probSyn * 1.0) + (probNonSyn * ratio);
    }
    
    
    public static void main(String[] args){
        String path = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/freq_test/sapiens.csv";
        CodonScaler scaler = new CodonScaler(new CodonFrequencies(path));
        //CodonScaler scaler = new CodonScaler(new CodonFrequencies()); // 1/64
        
        double ratio = 2.0;
        int siteType = 0;
        int frame = 0;
                
        System.out.println("get "+scaler.get(ratio, siteType, frame));
    }
    
}
