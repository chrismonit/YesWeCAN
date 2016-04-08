/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.codonawareness;

import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CodonSum {
    
    private CodonFrequencies codonFrequencies;
    
    private static int[][] otherCodonPositions = new int[][]{
        { 1, 2 }, //0
        { 0, 2 }, //1
        { 0, 1 }  //2
    }; // given some int i between 0 and 2 inclusive, access the other two ints between 0 and 2 inclusive
    
    private double[][] computedSums;
    
    public double getCodonSum(int codonPosition, int nucState){
        return this.computedSums[codonPosition][nucState];
    }
    
    public CodonSum(CodonFrequencies codonFrequencies){
        this.codonFrequencies = codonFrequencies;
        this.computedSums = new double[3][4];
        
        for (int iCodonPosition = 0; iCodonPosition < 3; iCodonPosition++) {
            
            for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {
                this.computedSums[iCodonPosition][iNucState] = sumCodon(iCodonPosition, iNucState, this.codonFrequencies);
            }
        }
        
    }// constructor
    
    
    private static double sumCodon(int positionInCodon, int nucState, CodonFrequencies codonFrequencies){
        int[] codonArray = new int[3];
        
        double sum = 0.0;
            
        for (int nBase = 0; nBase < States.NT_STATES; nBase++) { // nBase and mBase are the bases in the other two positions in the codon
            for (int mBase = 0; mBase < States.NT_STATES; mBase++) {

                codonArray[ positionInCodon ] = nucState;
                codonArray[ otherCodonPositions[positionInCodon][0] ] = nBase;
                codonArray[ otherCodonPositions[positionInCodon][1] ] = mBase;

                double codonFreq = codonFrequencies.getFrequency(ReorderFrequencies.alphaToPaml(codonArray)); // assuming here that bases are to be ordered ACGT, therefore we need to arrange as TCAG for interacting with CoodnFrequencies instance
                System.out.println("codon\t"+ArrayPrinter.toString(codonArray, ",")+"\tcodonFreq\t"+codonFreq);
                sum += codonFreq;
            }// for mBase
        } // for nBase
        return sum;
    }
    
    public static void main(String[] args){
        CodonFrequencies codonFrequencies = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/hiv.csv");
        CodonSum codonSum = new CodonSum(codonFrequencies);
        int positionInCodon = 0;
        int nucState = 2;
        System.out.println("sumCodon "+codonSum.sumCodon(positionInCodon, nucState, codonFrequencies));
        System.out.println("getSumCodon "+codonSum.getCodonSum(positionInCodon, nucState));
    }
    
}
