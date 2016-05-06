/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.sim;

import java.util.Random;
import pal.tree.Tree;
import yeswecan.cli.CommandArgs;
import yeswecan.model.submodels.CANModelMixture;
import yeswecan.phylo.GeneticStructure;
import yeswecan.run.RunCANMixture;

/**
 *
 * @author cmonit1
 */
public class SimFreqsMix extends SimFreqs {
        
        
    
        public SimFreqsMix(Tree tree, Random rand, GeneticStructure genStruct, 
            CommandArgs comArgs){
        
            super(tree, rand, genStruct, comArgs); 
            /* creates instance of FrequencySimulator rather than FrequencySimulatorMix
            Creates list of omegas based on -w argument, rather than -w1 etc  
            
            */          
            
            CANModelMixture canMix = RunCANMixture.makeMixture(comArgs, comArgs.getModel());
            
            
            
            // assing omegas
            // assign probabiltiies
            
            // construct freqsimmix instance
            
            // override header method
            // override get params method
            
        }
        
        
        
        
                
}
