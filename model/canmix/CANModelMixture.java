/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.canmix;

import java.util.ArrayList;
//import java.util.List;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
//import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.Probabilities;
//import yeswecan.model.SubstitutionModel;
import yeswecan.model.hky.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANModelMixture extends HKYModel {
    
//    private TsTvRatioAdvanced kappa;
//    private BaseFrequencies pi;
    //private List<Omega> omegas;
    //private List<List<Parameter>> omegaDistributions;
    private ArrayList<ArrayList<Omega>> omegaDistributions;
    private ArrayList<Probabilities> probabilities;
    private BranchScaling scaling;
    
//    public CANModelMixture(List<Parameter> parameters){
//        super.clearParameters();
//        super.setParameters(parameters);
//    }
    
    public CANModelMixture(TsTvRatioAdvanced kappa, BaseFrequencies pi, BranchScaling scaling, 
            ArrayList<ArrayList<Omega>> omegaDistributions, ArrayList<Probabilities> probabilities){
        // NB the 0th omega has to be an unoptimisible 1.0 value
        
        super(kappa, pi);       
        this.scaling = scaling;
        
        
        this.omegaDistributions = omegaDistributions;
        this.probabilities = probabilities;
        
        super.clearParameters();
        
        super.addParameters(super.getKappa(), super.getPi(), this.scaling);
        
        for (Probabilities p : this.probabilities){
            super.addParameters(p);
        }
        
        for (ArrayList<Omega> omegaList : this.omegaDistributions){
            for (Omega w : omegaList){
                super.addParameters(w);
            }
        }

    }
    
    
//    public TsTvRatioAdvanced getKappa() {
//        return kappa;
//    }
//
//    public BaseFrequencies getPi() {
//        return pi;
//    }
    
    
    public BranchScaling getScaling() {
        return scaling;
    }
    
    public Omega getOmega(int gene, int siteClass){
        return this.omegaDistributions.get(gene).get(siteClass);
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
