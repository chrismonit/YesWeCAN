
package yeswecan.model;

import com.sun.org.apache.bcel.internal.Constants;
import java.lang.ArithmeticException;

import java.lang.Math; //only for testing, creating a Q matrix
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import yeswecan.model.parameters.BaseFrequencies;
import yeswecan.model.parameters.TsTvRatio;
import yeswecan.phylo.States;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * class for housing various means of calculating exp(Qt)
 * Method of eigendecomposition of Q comes from Yang (2006) p68, as do the variable names
 */
public class YangEigenDecomposition implements ProbMatrixGenerator {
    
    private RealMatrix U;
    private RealMatrix lambda;
    private RealMatrix invU;

    
    public YangEigenDecomposition( RateMatrix Q ){ 
               
        RealMatrix bigPiSqrt = makeBigPiSqrt( Q.getBaseFrequencies() );
        RealMatrix invBigPiSqrt = inverseBigPiSqrt( bigPiSqrt );
        
        RealMatrix B = makeB2( bigPiSqrt, Q, invBigPiSqrt );
        
        EigenDecomposition decompB = new EigenDecomposition(B);
        
        RealMatrix R = decompB.getV();
        this.lambda = decompB.getD(); //eigenvalues of B, which are the same as the eigenvalues of Q
        RealMatrix invR = decompB.getVT(); // inverse(R) == transpose(R)
        
        // could throw an error if decompB.hasComplexEigenValues
        if (decompB.hasComplexEigenvalues()) {
            throw new RuntimeException("YangEigenDecomposition: Matrix B has complex Eigenvalues");
        }
        
        
        this.U = invBigPiSqrt.multiply(R);
        this.invU = invR.multiply(bigPiSqrt);
 
        MatrixPrinter.PrintMatrix( bigPiSqrt.getData(), "bigPiSqrt");       
        MatrixPrinter.PrintMatrix( invBigPiSqrt.getData(), "invBigPiSqrt");
        MatrixPrinter.PrintMatrix( B.getData(), "B matrix, similar to Q (i.e. should have same eigenvalues");
    }    

    private static RealMatrix makeB(RateMatrix Q){
        double[][] B = new double[States.NT_STATES][States.NT_STATES];
        double[][] QData = Q.getData();
        
        double[] pi = Q.getBaseFrequencies().get();
        
        for (int i = 0; i < States.NT_STATES; i++) {
            for (int j = 0; j < States.NT_STATES; j++) {
                B[i][j] = QData[i][j] * Math.sqrt(pi[i]) * 1.0/Math.sqrt(pi[j]);
            }// j
        }// i
        
        return new Array2DRowRealMatrix(B);
    }// makeB
    
    
    private static RealMatrix makeBigPiSqrt(BaseFrequencies frequencies){

        double[][] piDiag = new double[frequencies.get().length][frequencies.get().length];
        for (int i = 0; i < frequencies.get().length; i++) {
            piDiag[i][i] = Math.sqrt( frequencies.get()[i] ); //may need to assert that the freqeuncies are all >= 0
        }
        return new Array2DRowRealMatrix(piDiag);
    }
    
    private static RealMatrix inverseBigPiSqrt( RealMatrix bigPiSqrt ){
        /* PI^(-1/2) = diag{ 1/sqrt(pi_A), 1/sqrt(pi_C), 1/sqrt(pi_G), 1/sqrt(pi_T) }
        */
        RealMatrix invBigPiSqrt = new Array2DRowRealMatrix( bigPiSqrt.getData() );
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
    
    
    private RealMatrix makeB2( RealMatrix bigPiSqrt,
                                        RealMatrix Q, 
                                        RealMatrix invBigPiSqrt){
    
       return bigPiSqrt.multiply(Q).multiply(invBigPiSqrt);
    }
    
    
    public RealMatrix getP(double t){
        return null;
    }
    
    
    
    
}//class
