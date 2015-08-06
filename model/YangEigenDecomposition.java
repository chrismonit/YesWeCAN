
package yeswecan.model;

import java.lang.ArithmeticException;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import yeswecan.model.parameters.BaseFrequencies;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.States;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * Method of eigendecomposition of Q comes from Yang (2006) p68, as do the variable names
 *    
 * Diagonalise B = R * Λ * R⁻¹ - Yang CME - pg. 69, eq. 2.27
 *
 * where B = Π^(½) * Q * Π^(-½)
 *
 * so
 *
 * Q = (Π^(-½) * R) * Λ * (R⁻¹ * Π^(½))       eq. 2.28
 *
 * and
 *
 * U   = Π^(-½) * R
 * U⁻¹ = R⁻¹ * Π^(½)
 *
 * Note, R⁻¹ == transpose(R)
 *
 * and so
 *
 * P(t) =  U * exp(Λ * t) * U⁻¹
 *
 */

public class YangEigenDecomposition implements ProbMatrixGenerator {
    
    private RateMatrix Q; // necessary to have a field for Q so we can have a getQ method
    
    private RealMatrix U;
    private RealMatrix lambda;
    private RealMatrix invU;
    
    public YangEigenDecomposition( RateMatrix Q ){ 
        
        this.Q = Q;
        
        RealMatrix B = makeB( Q );
        
        EigenDecomposition decompB = new EigenDecomposition(B);
        
        /**
         * According to Yang 2006 (citing Kelly 1979), 
         * it can be proved that a reversible Q matrix must have real eigenvalues.
         * But just in case:
         */
        if (decompB.hasComplexEigenvalues()) {
            throw new RuntimeException("YangEigenDecomposition: Matrix B has complex Eigenvalues");
        }
        
        // B = R * LAMBDA * R^-1
        
        RealMatrix R = decompB.getV();
        this.lambda = decompB.getD(); //eigenvalues of B, which are the same as the eigenvalues of Q
        RealMatrix invR = decompB.getVT(); // inverse(R) == transpose(R)
 
        double[] pi = Q.getBaseFrequencies().get();
        
        // U = PI^(-1/2) * R
        // U^-1 = R^-1 * PI^(1/2)
        
        double[][] UData = new double[States.NT_STATES][States.NT_STATES];
        double[][] invUData = new double[States.NT_STATES][States.NT_STATES];

        for (int i = 0; i < States.NT_STATES; i++) {
            for (int j = 0; j < States.NT_STATES; j++) {
                UData[i][j] = R.getEntry(i, j) * 1.0/Math.sqrt(pi[i]);
                invUData[i][j] = invR.getEntry(i, j) * Math.sqrt(pi[j]);
            }
        }
        this.U = new Array2DRowRealMatrix(UData);
        this.invU = new Array2DRowRealMatrix(invUData);
    }    

    private static RealMatrix makeB(RateMatrix Q){
        double[][] B = new double[States.NT_STATES][States.NT_STATES];
        double[][] QData = Q.getData();
        
        double[] pi = Q.getBaseFrequencies().get();
        
        for (int i = 0; i < States.NT_STATES; i++) {
            for (int j = 0; j < States.NT_STATES; j++) {
                B[i][j] = QData[i][j] * Math.sqrt(pi[i]) / Math.sqrt(pi[j]); 
            }// j
        }// i
        
        return new Array2DRowRealMatrix(B);
    }// makeB
    
    
//    private static RealMatrix makeBigPiSqrt(BaseFrequencies frequencies){
//
//        double[][] piDiag = new double[frequencies.get().length][frequencies.get().length];
//        for (int i = 0; i < frequencies.get().length; i++) {
//            piDiag[i][i] = Math.sqrt( frequencies.get()[i] ); //may need to assert that the freqeuncies are all >= 0
//        }
//        return new Array2DRowRealMatrix(piDiag);
//    }
//    
//    private static RealMatrix inverseBigPiSqrt( RealMatrix bigPiSqrt ){
//        /* PI^(-1/2) = diag{ 1/sqrt(pi_A), 1/sqrt(pi_C), 1/sqrt(pi_G), 1/sqrt(pi_T) }
//        */
//        RealMatrix invBigPiSqrt = new Array2DRowRealMatrix( bigPiSqrt.getData() );
//        for (int i = 0; i < bigPiSqrt.getRowDimension(); i++) {
//            double inverseOfDiagonalElement;
//            try{
//                inverseOfDiagonalElement = 1.0 / invBigPiSqrt.getEntry(i, i);
//            }
//            catch (ArithmeticException e){ 
//                System.out.println("Divide by zero error in TransitionProbabilityMatrix. Likely due to having a pi_x == 0."); 
//                inverseOfDiagonalElement = 0.0; 
//                System.exit(1); 
//            }
//            invBigPiSqrt.setEntry(i, i, inverseOfDiagonalElement);
//        }
//        return invBigPiSqrt;
//    }
//    
//    
//    private RealMatrix makeB2( RealMatrix bigPiSqrt,
//                                        RealMatrix Q, 
//                                        RealMatrix invBigPiSqrt){
//    
//       return bigPiSqrt.multiply(Q).multiply(invBigPiSqrt);
//    }
    
    
    public RealMatrix getP(double t){
        
        // apply function to lambda: exp(LAMBDA * t)
        
        double[][] lambdaTransformedData = new double[States.NT_STATES][States.NT_STATES];
        for (int i = 0; i < States.NT_STATES; i++) {
                lambdaTransformedData[i][i] = Math.exp( this.lambda.getEntry(i,i) * t );           
        }
        
        RealMatrix lambdaTransformed = new Array2DRowRealMatrix(lambdaTransformedData);
        
        // do multiplication to get P
        RealMatrix Pt = this.U.multiply( lambdaTransformed.multiply(this.invU) );
                
        return Pt;
    }
    
    public RateMatrix getQ(){
        return this.Q;
    }
    
    
    
}//class
