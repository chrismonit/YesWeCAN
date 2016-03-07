package yeswecan.sim;

import java.util.List;
import java.util.Random;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.GeneticStructure;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class FrequencySimulator {
    
    private Tree tree;
    private Random rand;
    private GeneticStructure genStruct;
    private List<Omega> omegas;
    private TsTvRatioAdvanced kappa;
    private List<CodonFrequencies> codonFrequencies;
    private CodonTable codonTable;
    
    public FrequencySimulator(Tree tree, String alignmentDestinationPath,
            Random rand, GeneticStructure genStruct, TsTvRatioAdvanced kappa,
            List<Omega> omegas, List<CodonFrequencies> codonFrequencies){
    
        this.tree = tree;
        this.rand = rand;
        this.genStruct = genStruct;
        this.omegas = omegas;
        this.kappa = kappa;
        this.codonFrequencies = codonFrequencies;
        this.codonTable = CodonTableFactory.createUniversalTranslator();
    }
    
    private double compute_r_ijl(int[] sequence, int j, int l){
        int[] genes = this.genStruct.getGenes(l);
        double product = 1.0; // multiplicative identity
        return Double.NEGATIVE_INFINITY;
        
        
        
    }
    
    
    public static void main(String[] args){
        CodonTable table = CodonTableFactory.createUniversalTranslator();
        for (int i = 0; i < 64; i++) {
            System.out.println(table.getAminoAcidCharFromCodonIndex(i));
        }
    }
    
    
    
}
