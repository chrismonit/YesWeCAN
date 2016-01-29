/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import java.util.ArrayList;
import java.util.List;
import pal.tree.Tree;
import swmutsel.model.parameters.Probabilities;
import yeswecan.cli.CommandArgs;
import yeswecan.model.can.CANModel;
import yeswecan.model.canmix.CANModelMixture;
import yeswecan.model.hky.HKYModel;
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
//        this.alignment = alignment;
//        this.tree = tree;
//        this. comArgs = input;
        
//        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
//                                                            this.comArgs.bFrame(),
//                                                            this.comArgs.cFrame(),
//                                                            this.comArgs.lengths());
        
        super(alignment, tree, input);
    }// constructor
    
    
    public CANModelMixture makeCANMixture(){
        
        HKYModel hky = RunHKY.makeHKY(comArgs);
        
        List<Probabilities> probabilities = new ArrayList<Probabilities>();
    
        
        
    }
    
        
        
        
        
        
}
