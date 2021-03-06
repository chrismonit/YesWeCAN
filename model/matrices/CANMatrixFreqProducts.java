/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.matrices;

import java.util.Arrays;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.datatype.Codons;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 */
public class CANMatrixFreqProducts extends RateMatrix {
    
    // indices in pentamer where codons start
    private final static int[][] codonStarts = new int[][]{ 
        // a frame, b frame, c frame
        { 2, 0, 1 }, // alpha site
        { 1, 2, 0 }, // beta site
        { 0, 1, 2 }  // gamma site
    };
    
    public CANMatrixFreqProducts(TsTvRatioAdvanced kappa, int siteType,
        Omega[] omegaArray, 
        CodonFrequencies[] codonFrequenciesArray, CodonTable codonTable){
        
        super(kappa, false); // we want to build on an unscaled K80 matrix
                
        // TODO this ought to be done in function class, outside of value method
        //Omega[] omegaArray = new Omega[]{w_A, w_B, w_C}; 
        
        // NB the unnormalised versions of these pi values are required for computing q_ij
        double[] unnormalisedContextProbSums = getUnnormalisedContextProbsSums(siteType, codonFrequenciesArray);
        double sumContextProbs = 0.0; // Z. Sum of all products of codon frequencies
        for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {
            sumContextProbs += unnormalisedContextProbSums[iNucState];
        }
        
        double[] normalisedBaseProbs = new double[States.NT_STATES];
        for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {
            normalisedBaseProbs[iNucState] = unnormalisedContextProbSums[iNucState]/sumContextProbs;
        }
        
        BaseFrequencies baseFreq = new BaseFrequencies();

        baseFreq.set(normalisedBaseProbs); // will normalise base frequencies
        super.setPi(baseFreq);
        
        //double[] piNormalised = baseFreq.get();
               
        for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {
            for (int jNucState = 0; jNucState < States.NT_STATES; jNucState++) {
                if (iNucState != jNucState){
                    //System.out.println("iNucState\t"+iNucState+"\tjNucState\t"+jNucState);
                    double q_ij = getQij(iNucState, jNucState, siteType, this.getKappa(), 
                            omegaArray, unnormalisedContextProbSums[iNucState], sumContextProbs, codonFrequenciesArray, codonTable
                    );
                    
                    this.setEntry(iNucState, jNucState, q_ij);
                }
            }
        }

        super.populateDiagonals();
          
    }// constructor
    
    
    // raw because not normalised
    // sum_{x_{-2}} \sum_{x_{-1}} \sum_{x_{1}} \sum_{x_{2}} F_{xxi} F_{xix} F_{ixx}
    public static double[] getUnnormalisedContextProbsSums(
        int siteType, CodonFrequencies[] codonFrequenciesArray){
        
        double[] unnormalisedContextProbs = new double[States.NT_STATES];
        
        for (int iState = 0; iState < unnormalisedContextProbs.length; iState++) {
            
            double sum = 0.0;
            
            for (int xm2 = 0; xm2 < States.NT_STATES; xm2++) {
                for (int xm1 = 0; xm1 < States.NT_STATES; xm1++) {
                    for (int xp1 = 0; xp1 < States.NT_STATES; xp1++) {
                        for (int xp2 = 0; xp2 < States.NT_STATES; xp2++) {

                            int[] pentamerI = new int[]{ xm2, xm1, iState, xp1, xp2 };
                            double product = 1.0;
                            for (int iFrame = 0; iFrame < 3; iFrame++) {
                                int[] codonI = Arrays.copyOfRange(
                                        pentamerI, codonStarts[siteType][iFrame], 
                                        codonStarts[siteType][iFrame]+3
                                );
                                
                                product *= codonFrequenciesArray[iFrame].getFrequency(ReorderFrequencies.alphaToPaml(codonI));

                            }// iFrame
                            
                            sum += product;
                        }//xp2
                    }//xp1
                }//xm1
            }//xm2
            unnormalisedContextProbs[iState] = sum;
        }// iState
        return unnormalisedContextProbs;
    }// getBaseFrequencyValues
    
    
    // TODO can probably make more efficient by separating sums of codon freqs from omegas
    private static double getQij(int iNucState, int jNucState, int siteType, 
        TsTvRatioAdvanced kappa, Omega[] omegas, double iSumUnnormalisedContextProb, double sumContextProbs,
        CodonFrequencies[] codonFrequenciesArray, CodonTable codonTable){

        double numerator = 0.0;

        //m = minus; p = plus
        for (int xm2 = 0; xm2 < States.NT_STATES; xm2++) {
            for (int xm1 = 0; xm1 < States.NT_STATES; xm1++) {
                for (int xp1 = 0; xp1 < States.NT_STATES; xp1++) {
                    for (int xp2 = 0; xp2 < States.NT_STATES; xp2++) {

                        double product = kappa.getKappaIfTransition(iNucState, jNucState);
                        
                        int[] pentamerI = new int[]{ xm2, xm1, iNucState, xp1, xp2 };
                        int[] pentamerJ = new int[]{ xm2, xm1, jNucState, xp1, xp2 };
                                            
                        for (int iFrame = 0; iFrame < 3; iFrame++) {

                            int[] codonI = Arrays.copyOfRange(
                                    pentamerI, codonStarts[siteType][iFrame], 
                                    codonStarts[siteType][iFrame]+3
                            );
                            int[] codonJ = Arrays.copyOfRange(
                                    pentamerJ, codonStarts[siteType][iFrame], 
                                    codonStarts[siteType][iFrame]+3
                            );
                            
                            int codonI_int = Codons.getCodonIndexFromNucleotideStates(codonI);
                            int codonJ_int = Codons.getCodonIndexFromNucleotideStates(codonJ);

                            
                            if (!codonTable.isSynonymous(codonI_int, codonJ_int)){
                                product *= omegas[iFrame].get();
                            }
                            product *= codonFrequenciesArray[iFrame].getFrequency(ReorderFrequencies.alphaToPaml(codonI));
                            product *= codonFrequenciesArray[iFrame].getFrequency(ReorderFrequencies.alphaToPaml(codonJ));
                            
                        }// iFrame

                        numerator += product;

                    }//xp2
                }// xp1
            }// xm1
        }// xm2
        
        return numerator/(sumContextProbs * iSumUnnormalisedContextProb);
  
    }// getQij
    
    
    
    public static void main(String[] args){
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(2.0);

        int siteType = 0;
        Omega w_A = new Omega(1.0);
        Omega w_B = new Omega(1.0);
        Omega w_C = new Omega(1.0);
        Omega[] omegas = new Omega[]{ w_A, w_B, w_C };
        BranchScaling scaling = new BranchScaling(3.0);
        CodonFrequencies codonFrequenciesHIV = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/hiv.csv");
        CodonFrequencies codonFrequencies64 = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/64.csv");
        CodonFrequencies codonFrequencies61 = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/61.csv");
        
        CodonTable table = CodonTableFactory.createUniversalTranslator();;
        
        CodonFrequencies[] codonFrequenciesArray = { codonFrequencies61, codonFrequencies61, codonFrequencies61 };

        CANMatrixFreqProducts can = new CANMatrixFreqProducts(kappa, siteType,
        omegas,
        codonFrequenciesArray, table);

        MatrixPrinter.PrintMatrix(can.getData(), "can data");
        
        System.out.println(can.getBaseFrequencies().toString());
        
        System.out.println("\n~ Fin ~");
        System.exit(1);
        
//        double[] rawBaseFreqs = CANMatrixFreqProducts.getUnnormalisedContextProbsSums(siteType, codonFrequenciesArray);
//        System.out.println("rawBaseFreqs\t"+ArrayPrinter.toString(rawBaseFreqs, ","));
//        double sum = 0.0;
//        for (int k = 0; k < rawBaseFreqs.length; k++) {
//            sum += rawBaseFreqs[k];
//        }
//        System.out.println("sum (Z) "+sum);
                
        
    }// main
    
    
}
