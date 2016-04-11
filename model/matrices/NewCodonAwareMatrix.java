/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.matrices;

import yeswecan.model.matrices.RateMatrix;
import yeswecan.model.matrices.CodonAwareMatrix;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.datatype.Codons;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
import yeswecan.model.codonawareness.ProportionScaler;
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
    
    
    
    
    public NewCodonAwareMatrix(TsTvRatioAdvanced kappa, int siteType,
        Omega w_A, Omega w_B, Omega w_C, BranchScaling scaling, CodonFrequencies codonFrequencies, CodonTable codonTable){
        
        super(kappa, false); // we want to build on an unscaled K80 matrix
        
        double[][] matrixData = this.getData();
       
        Omega[] omegas = new Omega[]{w_A, w_B, w_C}; // this ought to be done in function class, outside of value method, to avoid overhead
        
        for (int iState = 0; iState < States.NT_STATES; iState++) {
            for (int jState = 0; jState < States.NT_STATES; jState++) {
                if (iState != jState){
                    matrixData[iState][jState] *= getTermProducts(iState, jState, siteType, codonFrequencies, omegas, codonTable) * scaling.get();
                }
            }
        }
        super.setSubMatrix(matrixData, 0, 0); // replace data with new values

        //MatrixPrinter.PrintMatrix(matrixData, "matrix data after multiplying by new rates");
                
        super.populateDiagonals();
        
        double[] piValues = new double[States.NT_STATES];
        for (int i = 0; i < piValues.length; i++) {
            piValues[i] = computePi(i, siteType, codonFrequencies);
        }
        super.setPi(new BaseFrequencies(piValues));
        
        super.scale(); // needs pi values to be set first
    }// constructor
    
 
    
     
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
        NewCodonAwareMatrix can = new NewCodonAwareMatrix(kappa, siteType, w_A, w_B, w_C, scaling, codonFrequencies, table);
        
        int iNucState = 0;
        int jNucState = 2;
        //int siteType2 = 0;
        //int frame = 0;
                
        //double term = can.getTerm( iNucState,  jNucState,  siteType2,  frame, codonFrequencies, w_A);
        //System.out.println("can.getTerm "+term);
        
        Omega[] omegas = new Omega[]{ w_A, w_B, w_C };

        //double termProducts = can.getTermProducts(iNucState, jNucState, siteType, codonFrequencies, omegas, table);
        //.out.println("termProducts "+termProducts);

        
        System.out.println("\ntesting original CAN\n");
        
        CodonAwareMatrix origCan = new CodonAwareMatrix(kappa, pi, new ProportionScaler(), siteType, w_A, w_B, w_C, scaling);
        MatrixPrinter.PrintMatrix(origCan.getData(), "original CAN matrix");
        
        System.out.println("\nfin");
    }
    
}
