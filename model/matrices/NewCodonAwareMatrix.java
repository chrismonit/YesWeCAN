/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.matrices;

import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import yeswecan.model.codonawareness.CodonSum;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class NewCodonAwareMatrix extends RateMatrix {
    
    private final CodonSum codonSum;
    
    private final static int[][] codonPositions = new int[][]{
        // a frame, b frame, c frame
        { 0, 2, 1 }, // alpha site
        { 1, 0, 2 }, // beta site
        { 2, 1, 0 }  // gamma site
    };
    
    
    public NewCodonAwareMatrix(TsTvRatioAdvanced kappa, int siteType,
        Omega w_A, Omega w_B, Omega w_C, BranchScaling scaling, CodonSum codonSum){
        
        super(kappa, false); // we want to build on an unscaled K80 matrix
               
        this.codonSum = codonSum;
        
        Omega[] omegas = new Omega[]{w_A, w_B, w_C}; // TODO this ought to be done in function class, outside of value method, to avoid overhead
                
        for (int iNucState = 0; iNucState < States.NT_STATES; iNucState++) {
            for (int jNucState = 0; jNucState < States.NT_STATES; jNucState++) {
                if (iNucState != jNucState){
                    
                    double q_ij = getQij(iNucState, jNucState, siteType, this.getKappa(), scaling, omegas);
                    this.setEntry(iNucState, jNucState, q_ij);
                }
            }
        }
        
        //MatrixPrinter.PrintMatrix(this.getData(), "unscaled diagonals only");
        super.populateDiagonals();
        

        double[] pi = getNormalisedPiValues(codonSum);
        BaseFrequencies baseFreq = new BaseFrequencies();
        baseFreq.set(pi);
        super.setPi(baseFreq);
        
        super.scale(); // needs pi values to be set first, as scaling depends on eq freqs

    }// constructor
    
 
    
    private double getQij(int iNucState, int jNucState, int siteType, 
            TsTvRatioAdvanced kappa, BranchScaling scaling, Omega[] omegas){
        
        double product = 1.0;
        
        product *= kappa.getKappaIfTransition(iNucState, jNucState);
        product *= scaling.get();
        //System.out.println("kappa_if\t"+kappa.getKappaIfTransition(iNucState, jNucState)+"\tscaling\t"+scaling.get());
        //System.out.println("iNucState\t"+iNucState+"\tjNucState\t"+jNucState);
        for (int iFrame = 0; iFrame < 3; iFrame++) {
        
            int positionOfInterest = this.codonPositions[siteType][iFrame];
            
            double nonsynSum = this.codonSum.getCodonProductSum(positionOfInterest, iNucState, jNucState, false);
            double synSum = this.codonSum.getCodonProductSum(positionOfInterest, iNucState, jNucState, true);
            
            double numerator = (omegas[iFrame].get() * nonsynSum) + synSum;
            double denominator = this.codonSum.getCodonSum(positionOfInterest, iNucState);
            
            //System.out.println("iFrame\t"+iFrame+"\tposition\t"+positionOfInterest+"\tnonsynSum\t"+nonsynSum+"\tsynSum\t"+synSum+"\tomega\t"+omegas[iFrame].get()+"\tnumerator\t"+numerator+"\tdenominator\t"+denominator+"\tratio\t"+numerator/denominator);
            
            product *= (numerator/denominator);
            
        }// for iFrame
        //System.out.println("");
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
        System.out.println("new can class");
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(2.0);
        //double[] frequencies = new double[]{.1,.2,.3,.4};
        double[] frequencies = BaseFrequencies.getDefault();
        BaseFrequencies pi = new BaseFrequencies(frequencies);
        int siteType = 1;
        Omega w_A = new Omega(2.0);
        Omega w_B = new Omega(3.0);
        Omega w_C = new Omega(4.0);
        BranchScaling scaling = new BranchScaling(3.0);
        CodonFrequencies codonFrequencies = new CodonFrequencies("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/hiv.csv");
        
        CodonTable table = CodonTableFactory.createUniversalTranslator();;
        
        CodonSum codonSum = new CodonSum(codonFrequencies, table);
        NewCodonAwareMatrix can = new NewCodonAwareMatrix(kappa, siteType, w_A, w_B, w_C, scaling, codonSum);
        
        MatrixPrinter.PrintMatrix(can.getData(), "Q");
        
        
        System.out.println("\nfin");
    }
    
}
