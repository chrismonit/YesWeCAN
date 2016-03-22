/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.cli.CommandArgs;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimFreqs {
    
    private int nSubs;
    private int nRepeats;
    private double equiBranchLength;
    
    private Tree tree;
    private Random rand;
    private GeneticStructure genStruct;
    private TsTvRatioAdvanced kappa;
    private List<Omega> omegas;
    private List<CodonFrequencies> codonFrequencies;
    
    
    public SimFreqs(Tree tree, Random rand, GeneticStructure genStruct, 
            CommandArgs comArgs){
        
        this.tree = tree;
        this.genStruct = genStruct;
        this.rand = rand;
        
        this.kappa = new TsTvRatioAdvanced(comArgs.kappa());
        
        this.omegas = new ArrayList<Omega>();
        this.omegas.add(new Omega(1.0)); // neutral, for noncoding regions/frames
        for (double omegaValue : comArgs.omegas()){
            this.omegas.add( new Omega(omegaValue) ); 
        }
        
        
        
        // need comArgs to accept one of these
        this.codonFrequencies = new ArrayList<CodonFrequencies>();
        
        
        // initialise the stuff above
        
        
        FrequencySimulator simulator = new FrequencySimulator(
            this.tree, this.rand, this.genStruct, this.kappa, this.omegas, this.codonFrequencies
        );
        
    }
    
    
    
    
    
}
