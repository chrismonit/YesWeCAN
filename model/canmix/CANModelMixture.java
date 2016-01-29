/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.canmix;

//import java.util.ArrayList;
import java.util.List;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Probabilities;
import yeswecan.model.can.CANModel;
//import yeswecan.model.SubstitutionModel;
import yeswecan.model.hky.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANModelMixture extends CANModel {
    
//    private TsTvRatioAdvanced kappa;
//    private BaseFrequencies pi;
    //private List<Omega> omegas;
    //private List<List<Parameter>> omegaDistributions;
    //private ArrayList<ArrayList<Omega>> omegaDistributions;
    //private ArrayList<Omega> omegas
    private List<Probabilities> probabilities;
    int numSiteClasses; // used for accessing omegas
    //private BranchScaling scaling;
    
//    public CANModelMixture(List<Parameter> parameters){
//        super.clearParameters();
//        super.setParameters(parameters);
//    }
    
    public CANModelMixture(TsTvRatioAdvanced kappa, BaseFrequencies pi, BranchScaling scaling, 
            List<Omega> omegas, List<Probabilities> probabilities, int numSiteClasses){
        // NB the 0th omega has to be an unoptimisible 1.0 value
        
        super(kappa, pi, scaling, omegas);       
        //this.scaling = scaling;
        
        
        //this.omegaDistributions = omegaDistributions;
        this.probabilities = probabilities;
        this.numSiteClasses = numSiteClasses;
        
        super.clearParameters();
        
        super.addParameters(super.getKappa(), super.getPi(), super.getScaling());
        
        for (Probabilities p : this.probabilities){
            super.addParameters(p);
        }

    }
    
    
//    public TsTvRatioAdvanced getKappa() {
//        return kappa;
//    }
//
//    public BaseFrequencies getPi() {
//        return pi;
//    }
    
//    
//    public BranchScaling getScaling() {
//        return scaling;
//    }
    
    // need to implement linear index finder
    public Omega getOmega(int gene, int siteClass){
        //return this.omegaDistributions.get(gene).get(siteClass);
        int index = gene * this.numSiteClasses + siteClass;
        return super.getOmegas().get(index);
    }
    
    public double getProbability(int gene, int siteClass){
        return this.probabilities.get(gene).get()[siteClass];
    }
    
//    public List<Omega> getOmegas() {
//        return omegas;
//    }
    
    
    
    
    
//    public static void main(String[] args){
//        List<Omega> omegas = new ArrayList<Omega>();
//        omegas.add(new Omega(1.0));
//        omegas.add(new Omega(2.0));
//        omegas.add(new Omega(3.0));
//
//        
//        CANModelMixture can = new CANModelMixture(
//            new TsTvRatioAdvanced(2.0), new BaseFrequencies(new double[]{.1,.2,.3,.4}), new BranchScaling(1.0), omegas 
//        );
//        
//        for (Parameter p : can.getParameters()){
//            System.out.println(p.toString());
//        }
//    }
    
}
