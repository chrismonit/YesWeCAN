/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.can;

//import java.util.ArrayList;
import java.util.List;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
//import swmutsel.model.parameters.Parameter;
//import yeswecan.model.SubstitutionModel;
import yeswecan.model.hky.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANModel extends HKYModel {
    
//    private TsTvRatioAdvanced kappa;
//    private BaseFrequencies pi;
    private List<Omega> omegas;
    private BranchScaling scaling;
    
//    public CANModel(List<Parameter> parameters){
//        super.clearParameters();
//        super.setParameters(parameters);
//    }
    
    public CANModel(TsTvRatioAdvanced kappa, BaseFrequencies pi, BranchScaling scaling, List<Omega> omegas){
        // NB the 0th omega has to be an unoptimisible 1.0 value
        
//        this.kappa = kappa;
//        this.pi = pi;
        super(kappa, pi);
        this.omegas = omegas;
        this.scaling = scaling;
        super.clearParameters();
        super.addParameters(super.getKappa(), super.getKappa(), this.scaling);
        
        for (Omega w : this.omegas){
            super.addParameters(w);
        }

    }
    
    
//    public TsTvRatioAdvanced getKappa() {
//        return super.getKappa();
//    }
//
//    public BaseFrequencies getPi() {
//        return super.getPi();
//    }
    
    
    public BranchScaling getScaling() {
        return scaling;
    }
    
    public List<Omega> getOmegas() {
        return omegas;
    }
    
    
    
    
    
//    public static void main(String[] args){
//        List<Omega> omegas = new ArrayList<Omega>();
//        omegas.add(new Omega(1.0));
//        omegas.add(new Omega(2.0));
//        omegas.add(new Omega(3.0));
//
//        
//        CANModel can = new CANModel(
//            new TsTvRatioAdvanced(2.0), new BaseFrequencies(new double[]{.1,.2,.3,.4}), new BranchScaling(1.0), omegas 
//        );
//        
//        for (Parameter p : can.getParameters()){
//            System.out.println(p.toString());
//        }
//    
//    }
    
}
