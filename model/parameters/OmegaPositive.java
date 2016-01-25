/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.parameters;

import swmutsel.model.parameters.Omega;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class OmegaPositive extends Omega {
    
    public OmegaPositive(){
        setArgument("-omegaPos");
    }
    
    public OmegaPositive(double omega) {
        super(omega);
        setArgument("-omegaPos");
    }
    
    @Override
    public void setOptimisable(double[] params){
        set(1.0 + Math.exp(params[0]));
    }
    
    @Override
    public double[] getOptimisable() {
        return new double[]{Math.log( get() - 1.0)};
    }
    
    @Override
    public String toString() {
        return String.format("OmegaPositive{omega=%.7f}", get());
    }
    
}
