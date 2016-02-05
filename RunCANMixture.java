/*
 * To change super license header, choose License Headers in Project Properties.
 * To change super template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import pal.tree.Tree;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Probabilities;
import yeswecan.cli.CommandArgs;
import yeswecan.model.canmix.CANFunctionMixture;
import yeswecan.model.canmix.CANModelMixture;
import yeswecan.model.hky.HKYModel;
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
    
    public RunCANMixture(AdvancedAlignment alignment, Tree tree, CommandArgs input){
        this.alignment = alignment;
        this.tree = tree;
        this. comArgs = input;
        
        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                             this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
        
    }// constructor
    
    private static int numberSiteClasses(int mixtureModel){
        int numSiteClasses = -1;

        if (mixtureModel == Constants.M2_IDENTIFIER)
            numSiteClasses = Constants.NUM_M2_SITE_CLASSES;
        else
            numSiteClasses = Constants.NUM_M1_SITE_CLASSES;
        return numSiteClasses;
    }
    
    @Override
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "lnL", "kappa", "A", "C", "G", "T");
        for (int iGene = 0; iGene < this.comArgs.getGeneNumber(); iGene++) {
            for (int jClass = 0; jClass < 10; jClass++) {
                columns.add(Integer.toString(iGene) + "_" + Constants.OMEGA_STRING + Integer.toString(jClass));
            }
            
            for (int jClass = 0; jClass < 10; jClass++) {
                columns.add(Integer.toString(iGene) + "_" + Constants.PROB_STRING + Integer.toString(jClass));
            }
        }
        return columns.toArray(new String[columns.size()]);
    }
    
    
    @Override
    public double[] getInitialValues(){ // NB first element does not contain lnL
        ArrayList<Double> values = RunModel.getParameterValues(makeMixture(this.comArgs, this.comArgs.getModel()).getParameters());
        double[] resultArray = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            resultArray[i] = values.get(i);
        }
        return resultArray;
    }
    
    
    public CANModelMixture makeMixture(CommandArgs comArgs, int mixtureModel){
        
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
        
        if (mixtureModel == Constants.M2_IDENTIFIER){
            OmegaPositive neutralW_2 = new OmegaPositive(1.0); 
            neutralW_2.setOptimisable(false);
            omegas.add(neutralW_2);
            
            neutralProbs = new Probabilities(new double[]{ 0.0, 1.0, 0.0 }); // all density on w_1
        }else{
            neutralProbs = new Probabilities(new double[]{ 0.0, 1.0 });
        }
        
        neutralProbs.setOptimisable(false);
        probs.add(neutralProbs);
        
        // for coding frames
                
        for (int iGene = 0; iGene < comArgs.getGeneNumber(); iGene++) {
       
            OmegaNegative geneW_0 = new OmegaNegative(comArgs.omega0()[iGene]);
            // fix if needs fixing
            if (comArgs.fix().contains("0"+Constants.OMEGA_STRING+Integer.toString(iGene)))
               geneW_0.setOptimisable(false);
        
            omegas.add(geneW_0);

            Omega geneW_1 = new Omega(1.0);
            geneW_1.setOptimisable(false); // w_1 always fixed to 1
            omegas.add(geneW_1);
            
            Probabilities geneProbs;

            if (mixtureModel == Constants.M2_IDENTIFIER){
                OmegaPositive geneW_2 = new OmegaPositive(comArgs.omega2()[iGene]); 
                if (comArgs.fix().contains("2"+Constants.OMEGA_STRING+Integer.toString(iGene)))
                   geneW_2.setOptimisable(false);
                omegas.add(geneW_2);

                geneProbs = new Probabilities(new double[]{ comArgs.prob0()[iGene], comArgs.prob1()[iGene], comArgs.prob2()[iGene]});
            }else{
                geneProbs = new Probabilities(new double[]{ comArgs.prob0()[iGene], comArgs.prob1()[iGene] });
            }
            
            if (comArgs.fix().contains(Constants.PROB_STRING+Integer.toString(iGene)))
                probs.add(geneProbs);
           
        } // for iGene
        
        return new CANModelMixture(hky, scaling, omegas, probs, numberSiteClasses(mixtureModel));
    }// make mixture
    
    
        @Override
    public double[] fit(){
  
        CANModelMixture canMix = makeMixture(this.comArgs, this.comArgs.getModel());
        CANFunctionMixture optFunction = 
                new CANFunctionMixture(this.alignment, this.tree, genStruct, canMix, 
                        numberSiteClasses(this.comArgs.getModel())
                );
        Optimise opt = new Optimise();
        CANModelMixture result = (CANModelMixture)opt.optNMS(optFunction, canMix);
        
        ArrayList<Double> values = RunModel.getParameterValues(result.getParameters());
        values.add(0, result.getLnL()); // prepend
        double[] resultArray = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            resultArray[i] = values.get(i);
        }
        return resultArray;
 
    }
    
    @Override
    public double[] calculate(){
                
        CANModelMixture canMix = makeMixture(this.comArgs, this.comArgs.getModel());
        double[] optimisableParams = Mapper.getOptimisable(canMix.getParameters()); // map parameters to optimisation space, so FunctionHKY.value canMix use them
        CANFunctionMixture calculator = 
                new CANFunctionMixture(this.alignment, this.tree, this.genStruct, 
                        canMix, numberSiteClasses(this.comArgs.getModel()) 
                );
        
        
        ArrayList<Double> values = RunModel.getParameterValues(canMix.getParameters());
        double lnL = calculator.value(optimisableParams);
        values.add(0, lnL);
        double[] resultArray = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            resultArray[i] = values.get(i);
        }
        return resultArray;
    }
         
    
    
    
}
