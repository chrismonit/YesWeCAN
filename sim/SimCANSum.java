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
import yeswecan.model.matrices.CodonAwareMatrix;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.codonawareness.RatioScaler;
import yeswecan.model.codonawareness.RatioScalerFactory;
import yeswecan.model.matrices.CANMatrixSum;
import yeswecan.model.submodels.CANModelSum;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimCANSum extends SimModel {
    
    
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModelSum canSum;        
    private Random rand;
    private AlignmentBuilder siteStates;
    private boolean printSubCounts;

    
    public SimCANSum(Tree tree, Random rand, CANModelSum can, GeneticStructure genStruct, boolean printSubCounts){
      
        this.tree = tree;
        this.genStruct = genStruct;
        this.rand = rand;
        this.canSum = can;
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
            Omega aOmega = this.canSum.getOmega(genes[0]); // genes[0] is the gene present in frame A
            Omega bOmega = this.canSum.getOmega(genes[1]);  
            Omega cOmega = this.canSum.getOmega(genes[2]);  
                        
            // make model
            
            RatioScaler ratioScaler = RatioScalerFactory.getRatioScaler();
            CANMatrixSum canQ = new CANMatrixSum(this.canSum.getKappa(), siteType, aOmega, bOmega, cOmega, this.canSum.getScaling());            
            //MatrixPrinter.PrintMatrix(canQ.getData(), "Q sim canSum", "");

            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(canQ);
            //MatrixPrinter.PrintMatrix(Pgen.getP(0.2).getData(), "P(0.2)");
            // simulate according to process
                        
            int rootState = States.draw(canQ.getBaseFrequencies().get(), rand.nextDouble());

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
