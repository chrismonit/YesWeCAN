/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model;

import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.datatype.Codons;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
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
        Omega w_A, Omega w_B, Omega w_C, BranchScaling scaling, CodonFrequencies codonFrequencies, CodonTable codonTable){
        
        super(kappa, pi);
        
        // matrix is now an HKY
        // need to extend to CAN properly
        
        double[][] matrixData = this.getData();

        MatrixPrinter.PrintMatrix(matrixData, "matrix data while still an HKY");
       
        Omega[] omegas = new Omega[]{w_A, w_B, w_C}; // this ought to be done in function class, outside of value method, to avoid overhead
        
        for (int iState = 0; iState < States.NT_STATES; iState++) {
            for (int jState = 0; jState < States.NT_STATES; jState++) {
                matrixData[iState][jState] *= getTermProducts(iState, jState, siteType, codonFrequencies, omegas, codonTable) * scaling.get();
            }
        }
        
        MatrixPrinter.PrintMatrix(matrixData, "matrix data after multiplying by new rates");

        super.setSubMatrix(matrixData, 0, 0); // replace data with new values
    }
    
    public static double getTermProducts(int iNucState, int jNucState, int siteType, 
            CodonFrequencies codonFrequencies, Omega[] omegas, CodonTable codonTable){
        double product = 1.0;
        
        for (int iFrame = 0; iFrame < 3; iFrame++) {
            product *= ( getTerm(iNucState, jNucState, siteType, iFrame, codonFrequencies, omegas[iFrame], codonTable) ) ; 
            //System.out.println("\n");
        }
        
        return product;
    }
    
    
    public static double getTerm(int iNucState, int jNucState, int siteType, int frame, 
            CodonFrequencies codonFrequencies, Omega omega, CodonTable codonTable){ 
        
        double numerator = 0.0; // equivalent to flux from codon I to codon J
        double denominator = 0.0; // equivalent to probability of codon I
        
        int[] iCodon = new int[3];
        int[] jCodon = new int[3];
        
        for (int nBase = 0; nBase < States.NT_STATES; nBase++) { // NB nBase/mBase have nothing to do with iNucState/jNucState. Iterating over the other bases in the codon, outside of the codon position of interest
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
                
                
                double omegaTerm = 1.0;
                int codonI_int = Codons.getCodonIndexFromNucleotideStates(iCodon);
                int codonJ_int = Codons.getCodonIndexFromNucleotideStates(jCodon);

                //System.out.println("codonI_int "+codonI_int+" codonJ_int "+codonJ_int);
                
                if (!codonTable.isSynonymous(codonI_int, codonJ_int)) { 
                    omegaTerm *= omega.get();
                }
                //System.out.println("omegaTerm "+omegaTerm);
                //System.out.println("iCodon "+ArrayPrinter.toString(iCodon, ",")+" freq "+iCodonFreq+"\t\tjCodon "+ArrayPrinter.toString(jCodon, ",")+" freq "+jCodonFreq+"\tomegaTerm "+omegaTerm);

                numerator += omegaTerm * iCodonFreq * jCodonFreq;
                denominator += iCodonFreq;
                //System.out.println("");

            }
            //System.out.println("---");
        }
        //System.out.println("numerator "+numerator + " denominator "+denominator+ " ratio "+numerator/denominator);
        //System.out.println("");

        return numerator / (denominator + Constants.EPSILON); // adding small term to avoid an unlikely divide by zero problem 
    }
    
    
    public static double computePi(int iNucState, int siteType,  
        CodonFrequencies codonFrequencies, CodonTable codonTable){
        
        double product = 1.0;
        int[] iCodon = new int[3];
        
        for (int iFrame = 0; iFrame < 3; iFrame++) {
            double sum = 0.0;
            
            for (int nBase = 0; nBase < States.NT_STATES; nBase++) { // NB nBase/mBase have nothing to do with iNucState/jNucState. Iterating over the other bases in the codon, outside of the codon position of interest
                for (int mBase = 0; mBase < States.NT_STATES; mBase++) {

                    int codonPositionOfInterest = codonPositions[siteType][iFrame]; // the site of the hypothetical quintuplet which is actually changing

                    iCodon[ codonPositionOfInterest ] = iNucState;
                    iCodon[ otherCodonPositions[codonPositionOfInterest][0] ] = nBase;
                    iCodon[ otherCodonPositions[codonPositionOfInterest][1] ] = mBase;

                    double iCodonFreq = codonFrequencies.getFrequency(ReorderFrequencies.alphaToPaml(iCodon));
                    sum += iCodonFreq;
                }// for mBase
            } // for nBase
            product *= sum;
        }// for iFrame
        
        return product;
    }
    
    
    public static void main(String[] args){
        System.out.println("hello world");
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(2.0);
        //double[] frequencies = new double[]{.1,.2,.3,.4};
        double[] frequencies = BaseFrequencies.getDefault();
        BaseFrequencies pi = new BaseFrequencies(frequencies);
        int siteType = 2;
        Omega w_A = new Omega(2.0);
        Omega w_B = new Omega(3.0);
        Omega w_C = new Omega(4.0);
        BranchScaling scaling = new BranchScaling(1.0);
        CodonFrequencies codonFrequencies = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/hiv.csv");
        
        CodonTable table = CodonTableFactory.createUniversalTranslator();;
        NewCodonAwareMatrix can = new NewCodonAwareMatrix(kappa, pi, siteType, w_A, w_B, w_C, scaling, codonFrequencies, table);
        
        int iNucState = 0;
        int jNucState = 2;
        //int siteType2 = 0;
        //int frame = 0;
                
        //double term = can.getTerm( iNucState,  jNucState,  siteType2,  frame, codonFrequencies, w_A);
        //System.out.println("can.getTerm "+term);
        
        Omega[] omegas = new Omega[]{ w_A, w_B, w_C };

        double termProducts = can.getTermProducts(iNucState, jNucState, siteType, codonFrequencies, omegas, table);
        System.out.println("termProducts "+termProducts);

        
        System.out.println("\ntesting original CAN\n");
        
        CodonAwareMatrix origCan = new CodonAwareMatrix(kappa, pi, new ProportionScaler(), siteType, w_A, w_B, w_C, scaling);
        MatrixPrinter.PrintMatrix(origCan.getData(), "original CAN matrix");
        
        System.out.println("\nfin");
    }
    
}
