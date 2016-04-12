/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.run;

import java.util.ArrayList;
import java.util.List;
import pal.tree.Tree;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
import yeswecan.cli.CommandArgs;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.model.submodels.CANModelSum;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunCANSum extends RunModel {
    
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    protected GeneticStructure genStruct;
    
    public RunCANSum(AdvancedAlignment alignment, Tree tree, CommandArgs input){
        
        this.comArgs = input;
        this.alignment = alignment;
        this.tree = tree;
        
        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                            this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
    }
    
    
    protected static CANModelSum makeCANSum(CommandArgs comArgs){
    
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(comArgs.kappa());
        if (comArgs.fix().contains(Constants.FIX_KAPPA)) {
            kappa.setOptimisable(false);
        }
        
        BranchScaling scaling = new BranchScaling(comArgs.scaling());
        if (comArgs.fix().contains(Constants.FIX_SCALING)){
            scaling.setOptimisable(false);
        }
        
                List<Omega> omegas = new ArrayList<Omega>();
        
        Omega neutralOmega = new Omega(1.0); // for frames where there is no gene
        neutralOmega.setOptimisable(false); // never want this to change in optimisation
        omegas.add(neutralOmega);
        
//         positions of omegas correspond to the genes they represent (i.e. gene 1 omega is first)
        for (int i = 0; i < comArgs.omegas().length; i++) {
            double omegaValue = comArgs.omegas()[i];
            Omega w = new Omega( omegaValue );
            if ( comArgs.fix().contains(Integer.toString(i+1)) ) { // +1 because 0th omega is neutral
                w.setOptimisable(false);
            }
            omegas.add(w);
        }
        
        //return new CANModelSum(kappa, )
    }
    
}
