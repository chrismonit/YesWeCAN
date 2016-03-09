package yeswecan.sim;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.datatype.Codons;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.utils.ArrayPrinter;

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
    
    // TODO make this private after testing
    public double computeRate(int[] quintStates, int j, int site){
        // r_ijl = k * w_A * w_B * w_C * π_A * π_B * π_C
        
        int[] genes = this.genStruct.getGenes(site);
        double product = 1.0; // multiplicative identity
        
        //ArrayPrinter.print(genes, ",");
        
        product *= this.kappa.getKappaIfTransition(quintStates[2], j);
        for (int iFrame = 0; iFrame < 3; iFrame++) {
            int[] codonI_array = getCodon(quintStates, quintStates[2], iFrame, site%3);
            int codonI = Codons.getCodonIndexFromNucleotideStates(codonI_array);

            
            int[] codonJ_array = getCodon(quintStates, j, iFrame, site%3);

            int codonJ = Codons.getCodonIndexFromNucleotideStates(codonJ_array);

            if (!this.codonTable.isSynonymous(codonI, codonJ)) { 
                product *= this.omegas.get(genes[iFrame]).get();
            }
            
            
            CodonFrequencies geneCodonFreq = this.codonFrequencies.get(genes[iFrame]);
            int[] mappedToPaml = ReorderFrequencies.alphaToPaml(codonJ_array);
            double pi_J = geneCodonFreq.getFrequency(mappedToPaml); 
            
            product *= pi_J;
            
            
        }// iFrame
        return product;
    }
    
    
     /*Matrix A is accessed using the site type and the frame you want the codon for.
     * Each row off matrix B contains the co-ordinates of the nt positions for the codon you want, relative to position
     * Using the two in combination, we get the characters in the sequence which
     * correspond to the codon we want for each frame
     */

    private static final int[][] A = {  //i = site type (alpha, beta, gamma), j = frame (a, b, c)
        { 2, 0, 1 },
        { 1, 2, 0 },
        { 0, 1, 2 }
    };

    private static final int[][] B = {   //represent nt positions in CODON in the quintuplet, relative to 'position' variable
        { -2, -1, 0 }, //1st codon in quint
        { -1, 0, 1 }, //2nd codon in quint
        { 0, 1, 2 }  //3rd codon in quint
    };
    
    
    /*
        provide array of length 5 with nuceltodie states, with site of interest in the centre (index=2)
        centralState is the nuceltide state at quint[2]. 
            If the quint is not mutated: centralState == quint[2]
            If quint is mutated: centralState == j
        siteType is siteIndex%3 (and obviously not quint[2]%3)
    */
    
    private static int[] getCodon(int[] quint, int centralState, int frame, int siteType){
        int[] quintCopy = new int[quint.length]; // we create a copy so we can change the central state if needed
        System.arraycopy(quint, 0, quintCopy, 0, quint.length);
        quintCopy[2] = centralState; // if centralState == quint[2], then quintCopy == quint
        
        int[] codon = new int[3];
        codon[0] = quintCopy[ 2 + B[ A[siteType][frame] ][0] ];
        codon[1] = quintCopy[ 2 + B[ A[siteType][frame] ][1] ];
        codon[2] = quintCopy[ 2 + B[ A[siteType][frame] ][2] ];
        return codon;
    }
    
    
    public static Tree loadTree(String treePath){  
        try{
            return new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println("Failed to read in tree for Simulator");
            e.printStackTrace();
            System.exit(1);
        }
        return null;
    }
    
    
}
