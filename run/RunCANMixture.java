/*
 * To change super license header, choose License Headers in Project Properties.
 * To change super template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.run;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import pal.tree.Tree;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.Probabilities;
import yeswecan.Constants;
import yeswecan.cli.CommandArgs;
import yeswecan.model.functions.CANFunctionMixture;
import yeswecan.model.submodels.CANModelMixture;
import yeswecan.model.submodels.HKYModel;
import yeswecan.model.parameters.OmegaNegative;
import yeswecan.model.parameters.OmegaPositive;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;


/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunCANMixture extends RunModel {
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    private GeneticStructure genStruct;
    
    private int numSiteClasses;
    
    public RunCANMixture(AdvancedAlignment alignment, Tree tree, CommandArgs input, int model){
        this.alignment = alignment;
        this.tree = tree;
        this. comArgs = input;
        
        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                             this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
        this.numSiteClasses = numberSiteClasses(model);
        
    }// constructor
    
    public static int numberSiteClasses(int mixtureModel){
        int numSiteClasses = -1;

        if (mixtureModel == Constants.M2_IDENTIFIER)
            numSiteClasses = Constants.NUM_M2_SITE_CLASSES;
        else
            numSiteClasses = Constants.NUM_M1_SITE_CLASSES;
        return numSiteClasses;
    }
    
//    public void fixed(){
//        CANModelMixture canMixInitial = makeMixture(this.comArgs, this.comArgs.getModel()); // NB will need to add complexModelIdentifier argument
//        for (Parameter p : canMixInitial.getParameters()){
//            System.out.print(p.toString());
//            System.out.print(" ");
//            System.out.println(p.isOptimisable());
//        }
//    }
    
    
    
    @Override
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "model", "lnL", "kappa", "A", "C", "G", "T", "sc");
        for (int iGene = 0; iGene < this.comArgs.getGeneNumber(); iGene++) {
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                columns.add(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.OMEGA_STRING + Integer.toString(jClass));
            }
            
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                columns.add(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.PROB_STRING + Integer.toString(jClass));
            }
        }
        return columns.toArray(new String[columns.size()]);
    }
    
    
    
    private double[] getValueArray(CANModelMixture canMix){
        List<Double> resultList = new ArrayList<Double>();
        resultList.add((double)this.comArgs.getModel());
        resultList.add(canMix.getLnL()); 


        resultList.add(canMix.getKappa().get());
        for (int i = 0; i < canMix.getPi().get().length; i++) {
            resultList.add(canMix.getPi().get()[i]);
        }
        resultList.add(canMix.getScaling().get());
        
        for (int iGene = 1; iGene < this.comArgs.getGeneNumber()+1; iGene++) {
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                resultList.add( canMix.getOmega(iGene, jClass).get() );
            }
            
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                resultList.add( canMix.getProbability(iGene, jClass) );
            }
        }
        
        double[] resultArray = new double[resultList.size()];
        for (int i = 0; i < resultList.size(); i++) {
            resultArray[i] = resultList.get(i).doubleValue();
        }
        return resultArray;

    }
    
    
 
    
    @Override
    public  double[] getInitialValues(){ // NB first element does not contain lnL
        return getValueArray( makeMixture(this.comArgs, this.comArgs.getModel(), Constants.M2_IDENTIFIER) );
    }
    
    
    public static CANModelMixture makeMixture(CommandArgs comArgs, int mixtureModel, int complexModelIdentifier){
        
        HKYModel hky = RunHKY.makeHKY(comArgs);
        
        BranchScaling scaling = new BranchScaling(comArgs.scaling());
        if (comArgs.fix().contains(Constants.FIX_SCALING)) {
            scaling.setOptimisable(false);
        }
        
        List<Probabilities> probs = new ArrayList<Probabilities>();
        List<Omega> omegas = new ArrayList<Omega>();
                
        // neutral (for noncoding frames)
        
        Probabilities neutralProbs;
       
        OmegaNegative neutralW_0 = new OmegaNegative(1.0); // should never be used
        neutralW_0.setOptimisable(false);
        omegas.add(neutralW_0);
        
        Omega neutralW_1 = new Omega(1.0);
        neutralW_1.setOptimisable(false);
        omegas.add(neutralW_1);
        
        if (mixtureModel == complexModelIdentifier){
            OmegaPositive neutralW_2 = new OmegaPositive(1.0); 
            neutralW_2.setOptimisable(false);
            omegas.add(neutralW_2);
            
            neutralProbs = new Probabilities(new double[]{ 0.0, 1.0, 0.0 }); // all density on w_1. w_1 == 1.0
        }else{
            neutralProbs = new Probabilities(new double[]{ 0.0, 1.0 });
        }
        
        neutralProbs.setOptimisable(false);
        probs.add(neutralProbs);
        
        // for coding frames
                
        for (int iGene = 0; iGene < comArgs.getGeneNumber(); iGene++) {
       
            OmegaNegative geneW_0 = new OmegaNegative(comArgs.omega0()[iGene]);
            // fix if needs fixing
            
            if (comArgs.fix().contains(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.OMEGA_STRING + Constants.SITE_CLASS_0)) //+1 for zero based
               geneW_0.setOptimisable(false);
        
            omegas.add(geneW_0);

            Omega geneW_1 = new Omega(1.0);
            geneW_1.setOptimisable(false); // w_1 always fixed to 1
            omegas.add(geneW_1);
            
            Probabilities geneProbs;

            if (mixtureModel == Constants.M2_IDENTIFIER){
                OmegaPositive geneW_2 = new OmegaPositive(comArgs.omega2()[iGene]); 
                if (comArgs.fix().contains(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.OMEGA_STRING + Constants.SITE_CLASS_2)) //+1 for zero based
                   geneW_2.setOptimisable(false);
                omegas.add(geneW_2);

                geneProbs = new Probabilities(new double[]{ comArgs.prob0()[iGene], comArgs.prob1()[iGene], comArgs.prob2()[iGene]});
            }else{
                geneProbs = new Probabilities(new double[]{ comArgs.prob0()[iGene], comArgs.prob1()[iGene] });
            }
            
            if (comArgs.fix().contains(Integer.toString(iGene+1)+Constants.WITIHIN_FIELD_SEPARATOR+Constants.PROB_STRING)) // zero based
                geneProbs.setOptimisable(false);
                
            probs.add(geneProbs);
           
        } // for iGene
        
        return new CANModelMixture(hky, scaling, omegas, probs, numberSiteClasses(mixtureModel)); // need to call numberSiteClasses rather than using numSiteClasses field so this method can be static
    }// make mixture
    
    
    @Override
    public double[] fit(){
        
        CANModelMixture canMix = makeMixture(this.comArgs, this.comArgs.getModel(), Constants.M2_IDENTIFIER);
        CANFunctionMixture optFunction = 
                new CANFunctionMixture(this.alignment, this.tree, genStruct, canMix, 
                        this.numSiteClasses
                );
        Optimise opt = new Optimise();
        long start = System.currentTimeMillis();
        CANModelMixture result = (CANModelMixture)opt.optNMS(optFunction, canMix);
        long time = System.currentTimeMillis() - start;
        
        double[] mles = getValueArray(result);        
        return mles;
    }
    
    @Override
    public double[] calculate(){
                
        CANModelMixture canMix = makeMixture(this.comArgs, this.comArgs.getModel(), Constants.M2_IDENTIFIER);
        double[] optimisableParams = Mapper.getOptimisable(canMix.getParameters()); // map parameters to optimisation space, so FunctionHKY.value canMix use them
        CANFunctionMixture calculator = 
                new CANFunctionMixture(this.alignment, this.tree, this.genStruct, 
                        canMix, this.numSiteClasses
                );
        
        
        double[] resultArray = getInitialValues();
        long start = System.currentTimeMillis();
        resultArray[1] = calculator.value(optimisableParams);
        long time = System.currentTimeMillis() - start;
        System.out.println("overall time (s): "+ time/1000.0);

        return resultArray;
    }
         
    
    
    
}
