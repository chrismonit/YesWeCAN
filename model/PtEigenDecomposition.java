
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
public class PtEigenDecomposition implements ProbMatrixGenerator {
    
    private RateMatrix Q;
    private EigenDecomposition decomp;
            
    public PtEigenDecomposition(RateMatrix Q){
        // perform eigen decomposition
        this.Q = Q;
        this.decomp = new EigenDecomposition(Q);

        //get rid of this later
//        RealMatrix diag = decomp.getD();
//        double[][] data = diag.getData();
//        MatrixPrintser.PrintMatrix(data, "diag");
        
    }
    
    
 
    public RealMatrix getP(double t){

        RealMatrix V = decomp.getV();
        
        //MatrixPrinter.PrintMatrix(V.getData(), "V");
        
        double[][] diag = decomp.getD().getData(); 
        //MatrixPrinter.PrintMatrix(diag, "diag before transform");

        
        for (int i = 0; i < diag.length; i++) {
            diag[i][i] = Math.exp( diag[i][i] * t );
        }
        
        //MatrixPrinter.PrintMatrix(diag, "diag after transform");
        
        RealMatrix transformD = new Array2DRowRealMatrix(diag);
        
        RealMatrix VT = decomp.getVT();
        
        //MatrixPrinter.PrintMatrix(VT.getData(), "VT");

        
        RealMatrix P_t = V.multiply(transformD.multiply(VT));
        
        //System.out.println("eigenvalue complex:   " + decomp.hasComplexEigenvalues());

        return P_t;
    }
    
    public RateMatrix getQ(){
        return this.Q;
    }    
    

    
}// class
