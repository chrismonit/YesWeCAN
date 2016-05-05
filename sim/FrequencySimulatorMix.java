/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.sim;

import java.util.List;
import java.util.Random;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public class FrequencySimulatorMix extends FrequencySimulator {
    
    public FrequencySimulatorMix(Tree tree, Random rand, GeneticStructure genStruct,
        TsTvRatioAdvanced kappa, List<Omega> omegas, List<CodonFrequencies> codonFrequencies){

        
        super(tree, rand, genStruct, kappa, omegas, codonFrequencies);
        
        
    }
    
}
