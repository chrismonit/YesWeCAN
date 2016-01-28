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
public class OmegaNegative extends Omega {
    
    public OmegaNegative(){
        setArgument("-omegaNeg");
    }
    
    public OmegaNegative(double omega) {
        super(omega);  
        setArgument("-omegaNeg");
    }
    
    // f(x) = \frac{1}{1+e^{-x}} (where x is the optimasble value and f(x) is the value used in the model (parameter space))
    @Override
    public void setOptimisable(double[] params){
        set( 1.0 / (1.0+Math.exp(-params[0])) );
    }
    
    // x = -\ln \frac{1}{f(x)-1} (where x is the optimisable value)
    @Override
    public double[] getOptimisable() {
        return new double[]{ -Math.log( (1.0/get()) - 1.0 )};
    }
    
    @Override
    public String toString() {
        return String.format("OmegaNegative{omega=%.7f}", get());
    }
    
}
