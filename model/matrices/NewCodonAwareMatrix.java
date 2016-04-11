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
import yeswecan.model.codonawareness.CodonSum;
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
    
    private CodonSum codonSum;
    private BranchScaling scaling;
    
    private static int[][] codonPositions = new int[][]{
        // a frame, b frame, c frame
        { 0, 2, 1 }, // alpha site
        { 1, 0, 2 }, // beta site
        { 2, 1, 0 }  // gamma site
    };
    
    
    public NewCodonAwareMatrix(TsTvRatioAdvanced kappa, int siteType,
        Omega w_A, Omega w_B, Omega w_C, BranchScaling scaling, CodonFrequencies codonFrequencies, CodonSum codonSum){
        
        super(kappa, false); // we want to build on an unscaled K80 matrix
               
        this.codonSum = codonSum;
        this.scaling = scaling;
        
        Omega[] omegas = new Omega[]{w_A, w_B, w_C}; // this ought to be done in function class, outside of value method, to avoid overhead
                
        //MatrixPrinter.PrintMatrix(this.getData(), "matrix data before doing anything");

        
        for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {
            for (int jNucState = 0; jNucState < States.NT_STATES; jNucState++) {
                if (iNucState != jNucState){
                    
                    double q_ij = getQij(iNucState, jNucState, siteType, this.getKappa(), scaling, omegas);
                    this.setEntry(iNucState, jNucState, q_ij);
                }
            }
        }
        
        
        //MatrixPrinter.PrintMatrix(this.getData(), "matrix data after setting q_ij");
                
        super.populateDiagonals();
        
        //MatrixPrinter.PrintMatrix(this.getData(), "matrix data after setting diagonals");


        double[] pi = getNormalisedPiValues(codonSum);
        super.setPi(new BaseFrequencies(pi));
        
        super.scale(); // needs pi values to be set first, as scaling depends on eq freqs
        //MatrixPrinter.PrintMatrix(this.getData(), "matrix data after scaling");

    }// constructor
    
 
    
    private double getQij(int iNucState, int jNucState, int siteType, 
            TsTvRatioAdvanced kappa, BranchScaling scaling, Omega[] omegas){
        
        double product = 1.0;
        
        product *= kappa.getKappaIfTransition(iNucState, jNucState);
        product *= scaling.get();
        
        for (int iFrame = 0; iFrame < 3; iFrame++) {
        
            int positionOfInterest = this.codonPositions[siteType][iFrame];
            
            double nonsynSum = this.codonSum.getCodonProductSum(positionOfInterest, iNucState, jNucState, false);
            double synSum = this.codonSum.getCodonProductSum(positionOfInterest, iNucState, jNucState, true);
            
            double numerator = (omegas[iFrame].get() * nonsynSum) + synSum;
            double denominator = this.codonSum.getCodonSum(positionOfInterest, iNucState);
            
            product *= (numerator/denominator);
            
        }// for iFrame
        
        return product;
    }
    
    
    private static double[] getNormalisedPiValues(CodonSum codonSum){
        double[] piValues = new double[States.NT_STATES];
        double sum = 0.0;
        for (int iNucState = 0; iNucState < piValues.length; iNucState++) {
            double pi = computeRawPi(iNucState, codonSum);
            piValues[iNucState] = pi;
            sum += pi;
        }
        
        for (int iNucState = 0; iNucState < piValues.length; iNucState++) {
            piValues[iNucState] /= sum;
        }
        
        return piValues;
    }
    
    private static double computeRawPi(int nucState, CodonSum codonSum){ // pi values before normalising
        
        double product = 1.0;
        for (int iCodonPosition = 0; iCodonPosition < 3; iCodonPosition++) {
            product *= codonSum.getCodonSum(iCodonPosition, nucState);
        }
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
        
        CodonSum codonSum = new CodonSum(codonFrequencies, table);
        NewCodonAwareMatrix can = new NewCodonAwareMatrix(kappa, siteType, w_A, w_B, w_C, scaling, codonFrequencies, codonSum);
        
        //System.out.println("pi:");
        //ArrayPrinter.print(can.getBaseFrequencies().get(), "\t");

//        System.out.println("\ntesting original CAN\n");
//        CodonAwareMatrix origCan = new CodonAwareMatrix(kappa, pi, new ProportionScaler(), siteType, w_A, w_B, w_C, scaling);
//        MatrixPrinter.PrintMatrix(origCan.getData(), "original CAN matrix");
        
        System.out.println("\nfin");
    }
    
}
