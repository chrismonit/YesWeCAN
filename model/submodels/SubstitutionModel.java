/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.submodels;

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
    
    private double lnL = Double.NaN; // keeping this blank so it must be overwritten

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

    public void setParameters(List<Parameter> p){
        this.parameters = p;
    }

    public List<Parameter> getParameters() {
        return this.parameters;
    }
    
    
}
