/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model;

import com.google.common.collect.Lists;
import java.util.Collections;
import java.util.List;
import swmutsel.model.parameters.Parameter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public abstract class SubstitutionModel {  // smaller version of equivalent class, swmutsel's SubstitutionModel.java:
    
    private List<Parameter> parameters = Lists.newArrayList();
    
    private double lnL = Double.NEGATIVE_INFINITY; // default lnL

    public double getLnL(){
        return this.lnL;
    }
    
    public void setLnL(double lnL){
        this.lnL = lnL;
    }
    
    // from swmutsel version:
    public void addParameters(Parameter... parameters) {
        Collections.addAll(this.parameters, parameters);
    }

    public void clearParameters() {
        this.parameters.clear();
    }


    public List<Parameter> getParameters() {
        return this.parameters;
    }
    
    
}
