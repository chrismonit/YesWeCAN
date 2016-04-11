/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.codonawareness;

import pal.datatype.CodonTable;
import pal.datatype.Codons;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;

import pal.datatype.CodonTableFactory; // only needed for testing
import yeswecan.utils.ArrayPrinter; // only needed for testing
/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CodonSum {
    
    
    private static int[][] otherCodonPositions = new int[][]{
        { 1, 2 }, //0
        { 0, 2 }, //1
        { 0, 1 }  //2
    }; // given some int i between 0 and 2 inclusive, access the other two ints between 0 and 2 inclusive
    
    private double[][] computedSumCodons;
    
    private double[][][][] computedSumCodonProducts;
    
    public double getCodonSum(int codonPosition, int nucState){
        return this.computedSumCodons[codonPosition][nucState];
    }
    
    public double getCodonProductSum(int codonPosition, int iNucState, int jNucState, boolean wantSynonymous){
        return this.computedSumCodonProducts[wantSynonymousInt(wantSynonymous)][codonPosition][iNucState][jNucState];
    }
    
    private boolean wantSynonymousBoolean(int wantSynonymousInt){ // 0 means false (you want nonsyn), 1 means true (you want syn)
        return (wantSynonymousInt == 0) ? false : true;
    }
    
    private int wantSynonymousInt(boolean wantSynonymousBoolean){
        return (wantSynonymousBoolean) ? 1 : 0;
    }
    
    public CodonSum(CodonFrequencies codonFrequencies, CodonTable codonTable){
        int CODON_POSITIONS = 3;
        int SYN_NONSYN_STATES = 2;
        
        // 1) summing over individual codons
        
        this.computedSumCodons = new double[CODON_POSITIONS][States.NT_STATES];
        
        for (int iCodonPosition = 0; iCodonPosition < CODON_POSITIONS; iCodonPosition++) {
            
            for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {
                this.computedSumCodons[iCodonPosition][iNucState] = sumCodon(iCodonPosition, iNucState, codonFrequencies);
            }
        }
        
        // 2) summing over codon product pairs, dependent on whether pair is synonymous or not
        // 4d array has more cells than needed, because we don't include j!=i pairs. But easier to have it larger than arrange which of the remaining j states go in which element
        this.computedSumCodonProducts = new double[SYN_NONSYN_STATES][CODON_POSITIONS][States.NT_STATES][States.NT_STATES];
        
        for (int iSynonymous = 0; iSynonymous < SYN_NONSYN_STATES; iSynonymous++) {

            for (int iCodonPosition = 0; iCodonPosition < CODON_POSITIONS; iCodonPosition++) {

                for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {

                    for (int jNucState = 0; jNucState < States.NT_STATES; jNucState++){
                        
                        if (iNucState != jNucState){
                            double sum = sumCodonProducts(iCodonPosition, iNucState, jNucState, codonFrequencies, codonTable, wantSynonymousBoolean(iSynonymous));
                            this.computedSumCodonProducts[iSynonymous][iCodonPosition][iNucState][jNucState] = sum;
                            System.out.println("\tiSynonymous\t"+iSynonymous+"\tiCodonPosition\t"+iCodonPosition+"\tiNucState\t"+iNucState+"\tjNucState\t"+jNucState+"\tsum\t"+sum);
                        }

                    }// j
                }// i
            }// codon position
        }// syn/nonsyn 
    }// constructor
    
    private static double sumCodon(int positionInCodon, int nucState, CodonFrequencies codonFrequencies){
        int[] codonArray = new int[3];
        
        double sum = 0.0;
        //System.out.println("codon pos\t"+positionInCodon);
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
                    //System.out.println("iCodon\t"+ArrayPrinter.toString(iCodon, ",")+"\tjCodon\t"+ArrayPrinter.toString(jCodon, ",")+"\tiCodonFreq\t"+iCodonFreq+"\tjCodonFreq\t"+jCodonFreq+"\tproduct\t"+(iCodonFreq*jCodonFreq));
                    sum += iCodonFreq * jCodonFreq;
                }
                
            }
        }
        return sum;
    }
    
    
    
    
    
    public static void main(String[] args){
        CodonFrequencies codonFrequencies = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/hiv.csv");
        CodonTable codonTable = codonTable = CodonTableFactory.createUniversalTranslator();

        CodonSum codonSum = new CodonSum(codonFrequencies, codonTable);
        int positionInCodon = 0;
        int nucState = 2;
        //System.out.println("sumCodon "+codonSum.sumCodon(positionInCodon, nucState, codonFrequencies));
        //System.out.println("getSumCodon "+codonSum.getCodonSum(positionInCodon, nucState));
        
        int iNucState = 2;
        int jNucState = 1;
        boolean wantSynonymous = true;
        
        System.out.println("getCodonProductSum "+codonSum.getCodonProductSum(positionInCodon, iNucState, jNucState, wantSynonymous));
        
    }// main
    
}
