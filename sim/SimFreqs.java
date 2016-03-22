/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import pal.alignment.Alignment;
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
    

    
    private Tree tree;
    private Random rand;
    private GeneticStructure genStruct;
    private TsTvRatioAdvanced kappa;
    private List<Omega> omegas;
    private List<CodonFrequencies> codonFrequencies;
    
    private int nSubs;
    private int nRepeats;
    private double equiBranchLength;
    
    private FrequencySimulator simulator;
    
    public SimFreqs(Tree tree, Random rand, GeneticStructure genStruct, 
            CommandArgs comArgs){
        
        this.tree = tree;
        this.rand = rand;
        this.genStruct = genStruct;
        
        this.kappa = new TsTvRatioAdvanced(comArgs.kappa());
        
        this.omegas = new ArrayList<Omega>();
        this.omegas.add(new Omega(1.0)); // neutral, for noncoding regions/frames
        for (double omegaValue : comArgs.omegas()){
            this.omegas.add( new Omega(omegaValue) ); 
        }        
        
        // only allowing one set of frequencies for all genes for now at least
        this.codonFrequencies = new ArrayList<CodonFrequencies>();
        this.codonFrequencies.add( new CodonFrequencies(comArgs.getCodonFrequencyPath()) );
        
        this.nSubs = comArgs.getNuNumSubs();
        this.nRepeats = comArgs.getNuNumRepeats();
        this.equiBranchLength = comArgs.getEquiBranchLength();
        
        this.simulator = new FrequencySimulator(
            this.tree, this.rand, this.genStruct, this.kappa, this.omegas, this.codonFrequencies
        );
 
    }// constructor
    
    
    public Alignment simulate(){
        double Z = 0.0; // sum nu samples
        for (int i = 0; i < nRepeats; i++) {
            int[] random = simulator.getRandomSequence(genStruct.getTotalLength());
            int[] changeStops = simulator.changeStop(random);
            Z += simulator.simulateNu(changeStops, nSubs);
        }
        double nu = Z/(double)nRepeats;
    
        int[] randomSeq = simulator.getRandomSequence(genStruct.getTotalLength());
        int[] changeStopSeq = simulator.changeStop(randomSeq);
        
        int[] rootSeq = simulator.evolveBranch(changeStopSeq, 10.0, nu);        
        
        // sanity check
        if (!simulator.sequenceAcceptable(rootSeq)){
            throw new RuntimeException("Root sequence for simulation contains stop codons");
        }
        
        Alignment result = simulator.simulate(rootSeq, nu);
        return result;
        
    }
    
    
    
}
