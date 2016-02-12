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
import pal.datatype.Nucleotides;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
import yeswecan.model.CodonAwareMatrix;
import yeswecan.model.ProbMatrixFactory;
import yeswecan.model.ProbMatrixGenerator;
import yeswecan.model.RatioScaler;
import yeswecan.model.RatioScalerFactory;
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

    
    public SimCANMixture(Tree tree, Random rand, CANModelMixture can, 
        GeneticStructure genStruct, boolean printSubCounts){
      
        this.tree = tree;
        this.genStruct = genStruct;
        this.rand = rand;
        this.canMix = can;
        this.printSubCounts = printSubCounts;

    }

    
    
    
    public Alignment simulate(){
        Alignment[] sites = new Alignment[this.genStruct.getTotalLength()];
        
        
        
        Node root = tree.getRoot();
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            // determine evolutionary process (from site type and partition)
        
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
        
            int aFrameSiteClass = SimModel.draw(this.canMix.getGeneProbabilities(genes[0]).get(), rand.nextDouble());
            int bFrameSiteClass = SimModel.draw(this.canMix.getGeneProbabilities(genes[1]).get(), rand.nextDouble());
            int cFrameSiteClass = SimModel.draw(this.canMix.getGeneProbabilities(genes[2]).get(), rand.nextDouble());
                        
            Omega aOmega = this.canMix.getOmega(genes[0], aFrameSiteClass); // genes[0] is the gene present in frame A
            Omega bOmega = this.canMix.getOmega(genes[1], bFrameSiteClass);  
            Omega cOmega = this.canMix.getOmega(genes[2], cFrameSiteClass);  
            
            System.out.println("MODEL"+Constants.DEL+"site"+Constants.DEL+"A"+Constants.DEL+"B"+Constants.DEL+"C"); // just a header for the class and w rows
            System.out.println("class"+Constants.DEL+iSite+Constants.DEL+aFrameSiteClass+Constants.DEL+bFrameSiteClass+Constants.DEL+cFrameSiteClass);
            System.out.println("w"+Constants.DEL+iSite+Constants.DEL+aOmega.get()+Constants.DEL+bOmega.get()+Constants.DEL+cOmega.get());
            
            RatioScaler ratioScaler = RatioScalerFactory.getRatioScaler();
            CodonAwareMatrix canQ = new CodonAwareMatrix(this.canMix.getKappa(), this.canMix.getPi(), ratioScaler, siteType, aOmega, bOmega, cOmega, this.canMix.getScaling());            
            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(canQ);
            
            // simulate according to process
            
            int rootState = SimModel.draw(this.canMix.getPi().get(), rand.nextDouble());

            this.siteStates = new AlignmentBuilder(this.tree.getExternalNodeCount()); // an 'alignment' for a single site, which will be populated with states by downTree.
            
            SubCount count = new SubCount();
            
            SimModel.downTree(tree, Pgen, root, rootState, count, this.siteStates, this.rand);
            
            if (this.printSubCounts){
                System.out.println(super.subCountHeader);
            }
            
            if (this.printSubCounts){
                System.out.println("count"+Constants.DEL +iSite + Constants.DEL + genStruct.getPartitionIndex(iSite) + Constants.DEL + iSite%3 + Constants.DEL + count.count); 
            }
            System.out.println("");

            // add this newly simulated site to the total set
            sites[iSite] = siteStates.generateAlignment(new Nucleotides());
                    
        }// for iSite
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        return allSites;
    }
    
}
