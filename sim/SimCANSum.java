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
import yeswecan.Constants;
import yeswecan.model.codonawareness.CodonSum;
import yeswecan.model.functions.CANFunctionSum;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.CANMatrixSum;
import yeswecan.model.submodels.CANModelSum;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;

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
    
    private CodonSum codonSum;
    private CANMatrixSum[][] Q_matrices;

    
    public SimCANSum(Tree tree, Random rand, CANModelSum canModelSum, GeneticStructure genStruct, boolean printSubCounts, CodonSum codonSum){
      
        this.tree = tree;
        this.genStruct = genStruct;
        this.rand = rand;
        this.canSum = canModelSum;
        this.printSubCounts = printSubCounts;
        
        this.codonSum = codonSum;
        
        this.Q_matrices = CANFunctionSum.createUnscaledMatrices(genStruct, canModelSum, codonSum);
        double nu = CANFunctionSum.computeNu(this.genStruct, this.Q_matrices);
        CANFunctionSum.scaleMatrices(this.genStruct, this.Q_matrices, nu);
        CANFunctionSum.scaleMatrices(this.genStruct, this.Q_matrices, this.canSum.getScaling().get());
        
    }
    
    public Alignment simulate(){
        Alignment[] sites = new Alignment[this.genStruct.getTotalLength()];
        
        
        if (this.printSubCounts){
            System.out.println(super.subCountHeader);
        }
        
        Node root = tree.getRoot();
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            // determine evolutionary process (from site type and partition)
            
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite % 3;
            
            CANMatrixSum Q = Q_matrices[partition][siteType];
            
            //MatrixPrinter.PrintMatrix(Q.getData(), "scaled Q site "+iSite);
            
            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(Q);
            // simulate according to process
                        
            int rootState = States.draw(Q.getBaseFrequencies().get(), rand.nextDouble());

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
