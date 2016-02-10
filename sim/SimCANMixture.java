/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;

import java.util.Random;
import pal.alignment.Alignment;
import pal.alignment.AlignmentBuilder;
import pal.alignment.ConcatenatedAlignment;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.model.canmix.CANModelMixture;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimCANMixture extends SimModel {
    
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModelMixture canMix;        
    private Random rand;
    private AlignmentBuilder siteStates;
    private boolean printSubCounts;

    
    public SimCANMixture(Tree tree, Random rand, CANModelMixture can, GeneticStructure genStruct, boolean printSubCounts){
      
        this.tree = tree;
        this.genStruct = genStruct;
        this.rand = rand;
        this.canMix = can;

    }
    
    public Alignment simulate(){
        Alignment[] sites = new Alignment[this.genStruct.getTotalLength()];
        
        if (this.printSubCounts){
            System.out.println(super.subCountHeader);
        }
        
        Node root = tree.getRoot();
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            // determine evolutionary process (from site type and partition)
        
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
        
            Omega aOmega = this.canMix.getOmega(genes[0]); // genes[0] is the gene present in frame A
            Omega bOmega = this.canMix.getOmega(genes[1]);  
            Omega cOmega = this.canMix.getOmega(genes[2]);  
            
            
            
        }// for iSite
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        return allSites;
    }
    
}
