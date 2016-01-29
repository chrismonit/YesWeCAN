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
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.Probabilities;
import yeswecan.cli.CommandArgs;
import yeswecan.model.can.CANModel;
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
    
    
    public CANModelMixture makeMixture(){
        
        HKYModel hky = RunHKY.makeHKY(comArgs);
        
        BranchScaling scaling = new BranchScaling(super.comArgs.scaling());
        if (super.comArgs.fix().contains(Constants.FIX_SCALING)) {
            scaling.setOptimisable(false);
        }
        
        List<Probabilities> probs = new ArrayList<Probabilities>();
        List<Omega> omegas = new ArrayList<Omega>();
        
        
        
        // neutral (for noncoding frames)
        
        
        Probabilities neutralProbs = new Probabilities(new double[]{0.0, 1.0, 0.0}); // all density on w_1 class, because w=1
        OmegaNegative neutralNeg = new OmegaNegative(0.0); // should never be used
        Omega neutralNeutral = new Omega(1.0);
        OmegaPositive neutralPos = new OmegaPositive(0.0); // should never be used
        probs.add(neutralProbs);
        omegas.add(neutralNeg);
        omegas.add(neutralNeutral);
        omegas.add(neutralPos);
        
        
        for (int i = 0; i < super.comArgs.getGeneNumber(); i++) {
            Probabilities geneProbs = new Probabilities(new double[]{ super.comArgs.prob0()[i], super.comArgs.prob1()[i], super.comArgs.prob2()[i]}); 
            OmegaNegative w_0 = new OmegaNegative(super.comArgs.omega0()[i]);
            
            Omega w_1 = new Omega(1.0);
            w_1.setOptimisable(false); // fix w_1 = 1 for every gene
            
            OmegaPositive w_2 = new OmegaPositive(super.comArgs.omega2()[i]);
            
            probs.add(geneProbs);
            omegas.add(w_0);
            omegas.add(w_1);
            omegas.add(w_2);
        }
        
        
        int siteClasses = 0; // CHANGE THIS!
        return new CANModelMixture(hky, scaling, omegas, probs, siteClasses);
    }
    
        
        
        
        
        
}
