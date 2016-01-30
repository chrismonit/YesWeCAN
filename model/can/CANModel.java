/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.can;

import java.util.List;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import yeswecan.model.hky.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANModel extends HKYModel {
    
    private List<Omega> omegas;
    private BranchScaling scaling;

    
    public CANModel(HKYModel hky, BranchScaling scaling, List<Omega> omegas){
        this(hky.getKappa(), hky.getPi(), scaling, omegas);
    }
    
    public CANModel(TsTvRatioAdvanced kappa, BaseFrequencies pi, BranchScaling scaling, List<Omega> omegas){
        // NB the 0th omega has to be an unoptimisible 1.0 value
        
        super(kappa, pi);
        this.omegas = omegas;
        this.scaling = scaling;
        super.addParameters(this.scaling);
        
        for (Omega w : this.omegas){
            super.addParameters(w);
        }

    }
    
    
    public BranchScaling getScaling() {
        return this.scaling;
    }
    
    public Omega getOmega(int gene) {
        return this.omegas.get(gene);
    }
    
    public List<Omega> getOmegas(){
        return this.omegas;
    }
    
}
