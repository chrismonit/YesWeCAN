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

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CodonAwareMatrix extends RateMatrix {
    
    public CodonAwareMatrix(TsTvRatioAdvanced kappa, BaseFrequencies pi, RatioScaler scaler, int siteType,
            Omega w_A, Omega w_B, Omega w_C, BranchScaling scaling){
        super(kappa, pi);
        
        // matrix is now an HKY
        // need to extend to CAN properly
        
        double[][] matrixData = this.getData();
        
        for (int i = 0; i < pi.get().length; i++) {
            for (int j = 0; j < pi.get().length; j++) {
                
                matrixData[i][j] *= 
                        scaler.get(w_A.get(), siteType, 0) *
                        scaler.get(w_B.get(), siteType, 1) *
                        scaler.get(w_C.get(), siteType, 2) *
                        (1.0/scaling.get());
       
            }
        }
        super.setSubMatrix(matrixData, 0, 0); // replace data with new values
        //this.scale(); // having changed the matrix content we need to scale it again
    }
    
}
