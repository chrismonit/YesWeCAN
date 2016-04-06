/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model;

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
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class NewCodonAwareMatrix extends RateMatrix {
    
    private static int[][] codonPositions = new int[][]{
        // a frame, b frame, c frame
        { 0, 2, 1 }, // alpha site
        { 1, 0, 2 }, // beta site
        { 2, 1, 0 }  // gamma site
    };
    
    private static int[][] otherCodonPositions = new int[][]{
        { 1, 2 }, //0
        { 0, 2 }, //1
        { 0, 1 }  //2
    }; // given some int i between 0 and 2 inclusive, access the other two ints between 0 and 2 inclusive
    
    
    
    public NewCodonAwareMatrix(TsTvRatioAdvanced kappa, BaseFrequencies pi, int siteType,
        Omega w_A, Omega w_B, Omega w_C, BranchScaling scaling, CodonFrequencies codonFrequencies){
        
        super(kappa, pi);
        
        // matrix is now an HKY
        // need to extend to CAN properly
        
        double[][] matrixData = this.getData();
                
        MatrixPrinter.Print3RowMatrix(matrixData, "matrix data while still an HKY");
       
        Omega[] omegas = new Omega[]{w_A, w_B, w_C};
        
        for (int iState = 0; iState < States.NT_STATES; iState++) {
            for (int jState = 0; jState < States.NT_STATES; jState++) {
                matrixData[iState][jState] *= getRate(iState, jState, siteType, codonFrequencies, omegas);
            }
        }
        
        MatrixPrinter.Print3RowMatrix(matrixData, "matrix data after multiplying by new rates");

        super.setSubMatrix(matrixData, 0, 0); // replace data with new values
    }
    
    public static double getRate(int iNucState, int jNucState, int siteType, CodonFrequencies codonFrequencies, Omega[] omegas){
        double product = 1.0;
        
        for (int iFrame = 0; iFrame < 3; iFrame++) {
            product *= ( omegas[iFrame].get() * getTerm(iNucState, jNucState, siteType, iFrame, codonFrequencies) ) ;            
        }
        
        return product;
    }
    
    
    public static double getTerm(int iNucState, int jNucState, int siteType, int frame, CodonFrequencies codonFrequencies ){ 
        double numerator = 0.0;
        double denominator = 0.0;
        
        int[] iCodon = new int[3];
        int[] jCodon = new int[3];
        
        for (int nBase = 0; nBase < States.NT_STATES; nBase++) { // NB nBase/mBase have nothing to do with iNucState/jNucState
            for (int mBase = 0; mBase < States.NT_STATES; mBase++) {
                
                int codonPositionOfInterest = codonPositions[siteType][frame]; // the site of the hypothetical quintuplet which is actually changing
                
                iCodon[ codonPositionOfInterest ] = iNucState;
                iCodon[ otherCodonPositions[codonPositionOfInterest][0] ] = nBase;
                iCodon[ otherCodonPositions[codonPositionOfInterest][1] ] = mBase;
                
                jCodon[ codonPositionOfInterest ] = jNucState;
                jCodon[ otherCodonPositions[codonPositionOfInterest][0] ] = nBase;
                jCodon[ otherCodonPositions[codonPositionOfInterest][1] ] = mBase;                
                
                double iCodonFreq = codonFrequencies.getFrequency(ReorderFrequencies.alphaToPaml(iCodon));
                double jCodonFreq = codonFrequencies.getFrequency(ReorderFrequencies.alphaToPaml(jCodon));
                
//                System.out.println("iCodon "+ArrayPrinter.toString(iCodon, ",")+" freq "+iCodonFreq);
//                System.out.println("\tjCodon "+ArrayPrinter.toString(jCodon, ",")+" freq "+jCodonFreq);
//                System.out.println("");

                numerator += iCodonFreq * jCodonFreq;
                
                denominator += iCodonFreq;
            }
        }
//        System.out.println("numerator "+numerator + " denominator "+denominator+ " ratio "+numerator/denominator);
        return numerator / denominator; // watch out for divde by zero?
    }
    
    
    
    public static void main(String[] args){
        System.out.println("hello world");
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(1.0);
        BaseFrequencies pi = new BaseFrequencies(BaseFrequencies.getDefault());
        int siteType = 0;
        Omega w_A = new Omega(2.0);
        Omega w_B = new Omega(2.0);
        Omega w_C = new Omega(2.0);
        BranchScaling scaling = new BranchScaling(1.0);
        CodonFrequencies codonFrequencies = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/hiv.csv");
        
        NewCodonAwareMatrix can = new NewCodonAwareMatrix(kappa, pi, siteType, w_A, w_B, w_C, scaling, codonFrequencies);
        
        int iNucState = 0;
        int jNucState = 1;
        int siteType2 = 0;
        int frame = 0;
        
        //can.getTerm( iNucState,  jNucState,  siteType2,  frame, codonFrequencies);
        
        System.out.println("\nfin");
    }
    
}
