/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;

import yeswecan.phylo.CodonFrequencies;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import pal.alignment.Alignment;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
import yeswecan.cli.CommandArgs;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimFreqs {
    

    
    protected Tree tree;
    protected Random rand;
    protected GeneticStructure genStruct;
    protected TsTvRatioAdvanced kappa;
    protected List<Omega> omegas;
    protected List<CodonFrequencies> codonFrequencies;
    
    protected int nSubs;
    protected int nRepeats;
    protected double equiBranchLength;
    protected boolean allowStops;
    
    protected int geneNumber; // just for get header and param values methods
    
    protected FrequencySimulator simulator;
    
    protected double meanNu; // for purposes of output only
    
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

        
        this.codonFrequencies = new ArrayList<CodonFrequencies>();
        codonFrequencies.add(new CodonFrequencies()); // default constructor has all freq = 1/64 for no gene case
        
        CodonFrequencies geneFrequencies = new CodonFrequencies(comArgs.getCodonFrequencyPath()); // only allowing one set of frequencies for all genes for now at least
        for (int iGene = 0; iGene < this.genStruct.getNumberOfGenes(); iGene++) {
            this.codonFrequencies.add( geneFrequencies ); // use reference to same instance for all genes
        }
                
        this.nSubs = comArgs.getNuNumSubs();
        this.nRepeats = comArgs.getNuNumRepeats();
        this.equiBranchLength = comArgs.getEquiBranchLength();
        this.allowStops = comArgs.getAllowStops();
        
        this.geneNumber = this.genStruct.getNumberOfGenes();
        
        this.simulator = new FrequencySimulator(
            this.tree, this.rand, this.genStruct, this.kappa, this.omegas, this.codonFrequencies
        );
 
    }// constructor
    
    
    public Alignment simulate(){
        
        // 1) Estimate Nu scaling parameter by simulation. Use mean of 'nRepeats' simulations
        
        double Z = 0.0; // sum
        for (int i = 0; i < nRepeats; i++) {
            int[] startSequence = simulator.getRandomSequence(genStruct.getTotalLength());
            
            if (!allowStops) {
                startSequence = simulator.changeStop(startSequence);
            }
            
            Z += simulator.simulateNu(startSequence, nSubs);
        }
        double nu = Z/(double)nRepeats;
        this.meanNu = nu;
        // 2) create root sequence
        
        int[] startSequence = simulator.getRandomSequence(genStruct.getTotalLength());
        
        if (!allowStops) {
            startSequence = simulator.changeStop(startSequence);
        }
        
        int[] rootSeq = simulator.evolveBranch(startSequence, equiBranchLength, nu);        
        
        // sanity check
        if (!simulator.sequenceAcceptable(rootSeq) && !allowStops){
            throw new RuntimeException("Root sequence for simulation contains stop codons");
        }
        
        // 3) simulate
        
        Alignment result = simulator.simulate(rootSeq, nu);
        return result;
        
    }
    
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "model", "kappa", "nuNumRepeats", "nuNumSubs", "equiBranchLength", "allowStops", "0_w");
        for (int i = 0; i < this.geneNumber; i++) {
            columns.add(Integer.toString(i+1) + "_" +Constants.OMEGA_STRING); // +1 for zero based correction
        }
        return columns.toArray(new String[columns.size()]);
        
    }
    
    public double[] getSimParameters(){
        double[] values = new double[ 7+this.geneNumber ]; // 4 are for the fields hard coded in getHeader
        values[0] = Constants.CODON_FREQ_IDENTIFIER;
        values[1] = this.kappa.get();
        values[2] = this.nRepeats;
        values[3] = this.nSubs;
        values[4] = this.equiBranchLength;
        values[5] = (this.allowStops) ? 1.0 : 0.0;
        
        for (int i = 0; i < this.geneNumber+1; i++) { // +1 here because I want to include no gene case
            values[i+6] = this.omegas.get(i).get();
        }
        return values;
    }
    
    public double getMeanNu(){
        return this.meanNu;
    }
    
}
