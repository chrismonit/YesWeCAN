
package yeswecan.model;

import java.lang.ArithmeticException;

import java.lang.Math; //only for testing, creating a Q matrix
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import yeswecan.model.parameters.BaseFrequencies;
import yeswecan.model.parameters.TsTvRatio;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * class for housing various means of calculating exp(Qt)
 * Method of eigendecomposition of Q comes from Yang (2006) p68, as do the variable names
 */
public class TransitionProbabilityMatrix {
    
    RealMatrix U;
    RealMatrix lambda;
    RealMatrix invU;

    
    public TransitionProbabilityMatrix( RateMatrix Q ){ 
               
        Array2DRowRealMatrix bigPiSqrt = makeBigPiSqrt( Q.getBaseFrequencies() );
        Array2DRowRealMatrix invBigPiSqrt = inverseBigPiSqrt( bigPiSqrt );
        
        Array2DRowRealMatrix B = makeB( bigPiSqrt, Q, invBigPiSqrt );
        
        EigenDecomposition decompositionB = new EigenDecomposition(B);
        
        RealMatrix R = decompositionB.getV();
        this.lambda = decompositionB.getD(); //eigenvalues of B, which are the same as the eigenvalues of Q
        RealMatrix invR = decompositionB.getVT(); // inverse(R) == transpose(R)
        
        // could throw an error if decompositionB.hasComplexEigenValues
        this.U = invBigPiSqrt.multiply(R);
        this.invU = invR.multiply(bigPiSqrt);
 
        MatrixPrinter.PrintMatrix( bigPiSqrt.getData(), "bigPiSqrt");
        
        MatrixPrinter.PrintMatrix( invBigPiSqrt.getData(), "invBigPiSqrt");

        MatrixPrinter.PrintMatrix( B.getData(), "B matrix, similar to Q (i.e. should have same eigenvalues");
    }//constructor
    
    //TODO
    public double[][] getTransitionProbabilityMatrix(double branchlength){
        return new double[0][0];
    
    
    };
    
    
    
    private Array2DRowRealMatrix makeBigPiSqrt(BaseFrequencies frequencies){

        double[][] piDiag = new double[frequencies.get().length][frequencies.get().length];
        for (int i = 0; i < frequencies.get().length; i++) {
            piDiag[i][i] = Math.sqrt( frequencies.get()[i] ); //may need to assert that the freqeuncies are all >= 0
        }
        return new Array2DRowRealMatrix(piDiag);
    }
    
    private Array2DRowRealMatrix inverseBigPiSqrt( Array2DRowRealMatrix bigPiSqrt ){
        /* PI^(-1/2) = diag{ 1/sqrt(pi_A), 1/sqrt(pi_C), 1/sqrt(pi_G), 1/sqrt(pi_T) }
        */
        Array2DRowRealMatrix invBigPiSqrt = new Array2DRowRealMatrix( bigPiSqrt.getData() );
        for (int i = 0; i < bigPiSqrt.getRowDimension(); i++) {
            double inverseOfDiagonalElement;
            try{
                inverseOfDiagonalElement = 1.0 / invBigPiSqrt.getEntry(i, i);
            }
            catch (ArithmeticException e){ 
                System.out.println("Divide by zero error in TransitionProbabilityMatrix. Likely due to having a pi_x == 0."); 
                inverseOfDiagonalElement = 0.0; 
                System.exit(1); 
            }
            invBigPiSqrt.setEntry(i, i, inverseOfDiagonalElement);
        }
        return invBigPiSqrt;
    }
    
    private Array2DRowRealMatrix makeB( Array2DRowRealMatrix bigPiSqrt,
                                        Array2DRowRealMatrix Q, 
                                        Array2DRowRealMatrix invBigPiSqrt){
    
       return bigPiSqrt.multiply(Q).multiply(invBigPiSqrt);
    }
    
    
    public static void main(String[] args){
        System.out.println("Testing TransitionProbabilityMatrix");
        TsTvRatio k = new TsTvRatio(1.0);
        BaseFrequencies freqs = new BaseFrequencies(new double[]{ 0.25, 0.25, 0.25, 0.25 });
        
        RateMatrix Q = new RateMatrix( k, freqs );
        
        TransitionProbabilityMatrix P = new TransitionProbabilityMatrix(Q);
        
        
        //start doing proper maths test to make sure this is working. refer to linear algebra theory!!
        
        
    }
    
    
    
    
    
}//class
