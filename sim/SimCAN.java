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
import yeswecan.model.can.CANModel;
import yeswecan.phylo.GeneticStructure;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimCAN extends SimModel {
    
    
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModel can;        
    private Random rand;
    private AlignmentBuilder siteStates;
    private boolean printSubCounts;

    
    public SimCAN(Tree tree, Random rand, CANModel can, GeneticStructure genStruct, boolean printSubCounts){
      
        this.tree = tree;
        this.genStruct = genStruct;
        this.rand = rand;
        this.can = can;
        this.printSubCounts = printSubCounts;

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
            Omega aOmega = this.can.getOmega(genes[0]); // genes[0] is the gene present in frame A
            Omega bOmega = this.can.getOmega(genes[1]);  
            Omega cOmega = this.can.getOmega(genes[2]);  
                        
            // make model
            
            RatioScaler ratioScaler = RatioScalerFactory.getRatioScaler();
            CodonAwareMatrix canQ = new CodonAwareMatrix(this.can.getKappa(), this.can.getPi(), ratioScaler, siteType, aOmega, bOmega, cOmega, this.can.getScaling());            
            //MatrixPrinter.PrintMatrix(canQ.getData(), "Q sim can", "");

            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(canQ);
            //MatrixPrinter.PrintMatrix(Pgen.getP(0.2).getData(), "P(0.2)");
            // simulate according to process
                        
            int rootState = SimModel.draw(this.can.getPi().get(), rand.nextDouble());

            this.siteStates = new AlignmentBuilder(this.tree.getExternalNodeCount()); // an 'alignment' for a single site, which will be populated with states by downTree.
            
            SubCount count = new SubCount();
            
            SimModel.downTree(tree, Pgen, root, rootState, count, this.siteStates, this.rand);
            
            if (this.printSubCounts){
                System.out.println(iSite + Constants.DEL + genStruct.getPartitionIndex(iSite) + Constants.DEL + iSite%3 + Constants.DEL + count.count);
            }
            
            // add this newly simulated site to the total set
            sites[iSite] = siteStates.generateAlignment(new Nucleotides());
        }
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        return allSites;
    }
    
    
    
    
}
