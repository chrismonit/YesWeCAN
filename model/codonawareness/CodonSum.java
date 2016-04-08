/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.codonawareness;

import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory; // only needed for testing
import pal.datatype.Codons;
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
                //System.out.println("codon\t"+ArrayPrinter.toString(codonArray, ",")+"\tcodonFreq\t"+codonFreq);
                sum += codonFreq;
            }// for mBase
        } // for nBase
        return sum;
    }
    
    /*
        In computing the numerators for q_ij, we sum over codons but include an omega term if
        the start and target states are nonsynonymous only:
        numerator = sum_{x_2} sum_{x_3} \omega_{f}^{0,1} \pi_{(i, x_2, x_3)} \pi_{(j, x_2, x_3)}
        To save computing the sums for each new omega, we can separate this sum into syn and nonsyn parts:
        
        nonsyn = \omega_{f}^{0,1} \times sum_{x_2} sum_{x_3}  \pi_{(i, x_2, x_3)} \pi_{(j, x_2, x_3) \text{if pair nonsynonymous}}
        syn =  sum_{x_2} sum_{x_3}  \pi_{(i, x_2, x_3)} \pi_{(j, x_2, x_3) \text{if pair synonymous}}
        This method lets you compute the sums over codons, and choose whether to include nonsyn or syn pairs.
        Omega term is not considered here.
    */
    private static double sumCodonProducts(int positionInCodon, int iNucState, int jNucState, 
        CodonFrequencies codonFrequencies, CodonTable codonTable, boolean wantSynonymousPairs){
        
        double sum = 0.0;
        
        int[] iCodon = new int[3];
        int[] jCodon = new int[3];
        
        for (int nBase = 0; nBase < States.NT_STATES; nBase++) {
            for (int mBase = 0; mBase < States.NT_STATES; mBase++) {
                
                iCodon[ positionInCodon ] = iNucState;
                iCodon[ otherCodonPositions[positionInCodon][0] ] = nBase;
                iCodon[ otherCodonPositions[positionInCodon][1] ] = mBase;
                
                jCodon[ positionInCodon ] = jNucState;
                jCodon[ otherCodonPositions[positionInCodon][0] ] = nBase;
                jCodon[ otherCodonPositions[positionInCodon][1] ] = mBase;                
                                
                int codonI_int = Codons.getCodonIndexFromNucleotideStates(iCodon);
                int codonJ_int = Codons.getCodonIndexFromNucleotideStates(jCodon);
                
                boolean pairSynonymous = codonTable.isSynonymous(codonI_int, codonJ_int);
                
                boolean synAndSyn = (wantSynonymousPairs && pairSynonymous); // we want synonymous codons AND this pair is synonymous
                boolean nonsynAndNonsyn = (!wantSynonymousPairs && !pairSynonymous); // we want nonsynonymous codons AND this pair is nonsynonymous
                                
                if (synAndSyn || nonsynAndNonsyn) { 
                    double iCodonFreq = codonFrequencies.getFrequency(ReorderFrequencies.alphaToPaml(iCodon));
                    double jCodonFreq = codonFrequencies.getFrequency(ReorderFrequencies.alphaToPaml(jCodon));
                    sum += iCodonFreq * jCodonFreq;
                }
                
            }
        }
        return sum;
    }
    
    
    
    
    
    public static void main(String[] args){
        CodonFrequencies codonFrequencies = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/hiv.csv");
        CodonSum codonSum = new CodonSum(codonFrequencies);
        int positionInCodon = 2;
        int nucState = 2;
        //System.out.println("sumCodon "+codonSum.sumCodon(positionInCodon, nucState, codonFrequencies));
        //System.out.println("getSumCodon "+codonSum.getCodonSum(positionInCodon, nucState));
        
        int iNucState = 0;
        int jNucState = 1;
        boolean wantSynonymous = false;
        
        CodonTable codonTable = codonTable = CodonTableFactory.createUniversalTranslator();
        
        System.out.println("sumCodonProducts sum "+codonSum.sumCodonProducts(positionInCodon, iNucState, jNucState, codonFrequencies, codonTable, wantSynonymous));
        
    }
    
}
