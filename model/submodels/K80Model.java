/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.submodels;

import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class K80Model extends SubstitutionModel {
    
    protected TsTvRatioAdvanced kappa;
    
    public K80Model(TsTvRatioAdvanced kappa){
        this.kappa = kappa;

        super.clearParameters();
        super.addParameters(kappa);
    }
    
    public TsTvRatioAdvanced getKappa() {
        return kappa;
    }
    
}
