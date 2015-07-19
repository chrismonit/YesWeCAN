
package yeswecan.model;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix; // may not need this by the end, useful in development
import yeswecan.utils.MatrixPrinter;

import java.lang.Math;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * created 8/7/15. This class differs from TransitionProbabilityMatrix in that it
 * will make use of library code to perform eigendecompsotion, rather than trying
 * to do this myself
 * 
 * from Yang (Comp Mol Evo book, 2006) p69
 * Q = V*D*(V^-1)
 * P(t) = exp(Qt) = U * exp(Dt) * (V^-1) = U * diag{ exp(d_1*t), exp(d_2*t) ... } * (V^-1)
 * where (V^-1) = V_transpose (Yang uses ^-1 while apache commons docs refers to transpose)
 */
public class ProbMatrixGenerator {
    
    
    private EigenDecomposition decomp;
            
    public ProbMatrixGenerator(RateMatrix Q){
        // perform eigen decomposition
        
        this.decomp = new EigenDecomposition(Q);

        //get rid of this later
//        RealMatrix diag = decomp.getD();
//        double[][] data = diag.getData();
//        MatrixPrinter.PrintMatrix(data, "diag");
        
    }
 
    public RealMatrix getP(double t){

        RealMatrix V = decomp.getV();
        
        double[][] diag = decomp.getD().getData();  
        for (int i = 0; i < diag.length; i++) {
            diag[i][i] = Math.exp( diag[i][i] * t );
        }
        RealMatrix transformD = new Array2DRowRealMatrix(diag);
        
        RealMatrix VT = decomp.getVT();
        
        RealMatrix P_t = V.multiply(transformD.multiply(VT));

        return P_t;
    }
    
}
