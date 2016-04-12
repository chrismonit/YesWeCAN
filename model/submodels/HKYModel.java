/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.submodels;

import swmutsel.model.parameters.BaseFrequencies;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 * 
 * Based on swmutsel's SwMut.java and its base class SubstitutionModel.java
 * 
 * It is the responsibility of this class to house the parameter instances
 * which comprise the model in question. The classes from swmutsel's Parameters
 * package require that such a class exists.
 * 
 * 
 * TODO this class can inherit from K80 model
 */
public class HKYModel extends SubstitutionModel {
    
    protected TsTvRatioAdvanced kappa;
    protected BaseFrequencies pi;
    
    
    public HKYModel(TsTvRatioAdvanced kappa, BaseFrequencies pi){
        this.kappa = kappa;
        this.pi = pi;
        super.clearParameters();
        super.addParameters(kappa, pi);
    }
    
    public TsTvRatioAdvanced getKappa() {
        return kappa;
    }

    public BaseFrequencies getPi() {
        return pi;
    }

}
