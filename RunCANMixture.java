/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import yeswecan.cli.CommandArgs;
import yeswecan.model.CANModelMixture;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunCANMixture {
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
    
    
    public CANModelMixture makeCANMixture(){
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(this.comArgs.kappa());
        if (this.comArgs.fix().contains(Constants.FIX_KAPPA)) {
            kappa.setOptimisable(false);
        }
      
        double[] frequencies = new double[States.NT_STATES]; // will be in correct order, whatever that may be
        
        if (Boolean.parseBoolean(this.comArgs.tcag())){
            frequencies = ReorderFrequencies.pamlToAlpha(this.comArgs.pi());
        }
        else{
            frequencies = this.comArgs.pi();
        }

        BaseFrequencies pi = new BaseFrequencies(frequencies);
        
        if (this.comArgs.fix().contains(Constants.FIX_FREQUENCIES)) {
            pi.setOptimisable(false);
        }
        
        BranchScaling scaling = new BranchScaling(this.comArgs.scaling());
        if (this.comArgs.fix().contains(Constants.FIX_SCALING)) {
            scaling.setOptimisable(false);
        }
    
    }
    
        
        
        
        
        
}
