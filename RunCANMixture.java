/*
 * To change super license header, choose License Headers in Project Properties.
 * To change super template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import java.util.ArrayList;
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
import yeswecan.phylo.AdvancedAlignment;


/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunCANMixture extends RunCAN {
//    private CommandArgs comArgs;
//    private AdvancedAlignment alignment;
//    private Tree tree;
    
//    private GeneticStructure genStruct;
    
    public RunCANMixture(AdvancedAlignment alignment, Tree tree, CommandArgs input){
//        super.alignment = alignment;
//        super.tree = tree;
//        super. comArgs = input;
        
//        super.genStruct = new GeneticStructure(super.comArgs.aFrame(),
//                                                            super.comArgs.bFrame(),
//                                                            super.comArgs.cFrame(),
//                                                            super.comArgs.lengths());
        
        super(alignment, tree, input); 
    }// constructor
    
    private static int numberSiteClasses(int mixtureModel){
        int numSiteClasses = -1;

        if (mixtureModel == Constants.M2_IDENTIFIER)
            numSiteClasses = Constants.NUM_M2_SITE_CLASSES;
        else
            numSiteClasses = Constants.NUM_M1_SITE_CLASSES;
        return numSiteClasses;
    }
    
    
    public CANModelMixture makeMixture(int mixtureModel){
        
        HKYModel hky = RunHKY.makeHKY(super.comArgs);
        
        BranchScaling scaling = new BranchScaling(super.comArgs.scaling());
        if (super.comArgs.fix().contains(Constants.FIX_SCALING)) {
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
                
        for (int iGene = 0; iGene < super.comArgs.getGeneNumber(); iGene++) {
       
            OmegaNegative geneW_0 = new OmegaNegative(super.comArgs.omega0()[iGene]);
            // fix if needs fixing
            if (super.comArgs.fix().contains("0"+Constants.FIX_OMEGA_STRING+Integer.toString(iGene)))
               geneW_0.setOptimisable(false);
        
            omegas.add(geneW_0);

            Omega geneW_1 = new Omega(1.0);
            geneW_1.setOptimisable(false); // w_1 always fixed to 1
            omegas.add(geneW_1);
            
            Probabilities geneProbs;

            if (mixtureModel == Constants.M2_IDENTIFIER){
                OmegaPositive geneW_2 = new OmegaPositive(super.comArgs.omega2()[iGene]); 
                if (super.comArgs.fix().contains("2"+Constants.FIX_OMEGA_STRING+Integer.toString(iGene)))
                   geneW_2.setOptimisable(false);
                omegas.add(geneW_2);

                geneProbs = new Probabilities(new double[]{ super.comArgs.prob0()[iGene], super.comArgs.prob1()[iGene], super.comArgs.prob2()[iGene]});
            }else{
                geneProbs = new Probabilities(new double[]{ super.comArgs.prob0()[iGene], super.comArgs.prob1()[iGene] });
            }
            
            if (super.comArgs.fix().contains(Constants.FIX_PROB_STRING+Integer.toString(iGene)))
                probs.add(geneProbs);
           
        } // for iGene
        
        return new CANModelMixture(hky, scaling, omegas, probs, numberSiteClasses(mixtureModel));
    }// make mixture
    
    
    public void calculateFixed(int mixtureModel){
                
        CANModelMixture canMix = makeMixture(mixtureModel);
        double[] optimisableParams = Mapper.getOptimisable(canMix.getParameters()); // map parameters to optimisation space, so FunctionHKY.value can use them
        CANFunctionMixture calculator = new CANFunctionMixture(super.alignment, super.tree, super.genStruct, canMix, numberSiteClasses(mixtureModel) );
        double lnL = calculator.value(optimisableParams);
        System.out.println("lnL: " + lnL + " "); // better to have it print the input parameters too, so you can see input and output together
    }
            
}
