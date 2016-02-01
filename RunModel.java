/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import java.util.ArrayList;
import java.util.List;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.Probabilities;
import swmutsel.model.parameters.TsTvRatio;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public abstract class RunModel {
    
    abstract String[] getHeader();
    
    abstract double[] fit();
    
    abstract double[] calculate();
    
    public static ArrayList<Double> getParameterValues(List<Parameter> params){
        Class<?>[] types = new Class<?>[4];
        types[0] = String.class;
        
        ArrayList<Double> values = new ArrayList<Double>();
        for (Parameter p : params){
            if (p instanceof TsTvRatio){
                values.add( ((TsTvRatio)p).get() );
            }
            else if (p instanceof BranchScaling){
                values.add( ((BranchScaling)p).get() );
            }
            else if (p instanceof Omega){
                values.add( ((Omega)p).get() );
            }
            else if (p instanceof BaseFrequencies){
                for (double d : ((BaseFrequencies)p).get())
                    values.add(d);
            }
            else if (p instanceof Probabilities ){
                for (double d : ((Probabilities)p).get())
                    values.add(d);
            }
            else{
                throw new RuntimeException("Warning: unrecognised parameter instance in RunModel");
            }
            
        }// for
        
        return values;
    }

}
