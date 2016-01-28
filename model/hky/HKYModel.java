/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.hky;

import com.google.common.collect.Lists; // this is within the swmutsel jar
import java.util.Collections;
import java.util.List;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.TsTvRatio;
import yeswecan.model.parameters.TsTvRatioAdvanced;

import yeswecan.Constants;
import yeswecan.model.SubstitutionModel;

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
 */
public class HKYModel extends SubstitutionModel {
    
    
    private TsTvRatioAdvanced kappa;
    private BaseFrequencies pi;
    
    
    public HKYModel(){
        this(new TsTvRatioAdvanced(Constants.DEFAULT_KAPPA), 
                new BaseFrequencies(Constants.DEFAULT_PI));
    }
    
    public HKYModel(TsTvRatioAdvanced kappa, BaseFrequencies pi){
        this.kappa = kappa;
        this.pi = pi;
        super.clearParameters();
        super.addParameters(kappa, pi);
    }
    
    //TODO this is incorrent - using other methods you set kappa and pi, but using this constructor, you don't
    public HKYModel(List<Parameter> p){
        super.clearParameters();
        super.setParameters(p);
    }
    
//    public HKYModel(HKYModel model){
//        this.kappa = model.kappa;
//        this.pi = model.pi;
//        super.clearParameters();
//        super.addParameters(kappa, pi);
//    }
    
    
    public TsTvRatioAdvanced getKappa() {
        return kappa;
    }

    public BaseFrequencies getPi() {
        return pi;
    }

}
