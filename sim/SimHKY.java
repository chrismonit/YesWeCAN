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
import swmutsel.model.parameters.BaseFrequencies;
import yeswecan.Constants;
import yeswecan.model.ProbMatrixFactory;
import yeswecan.model.ProbMatrixGenerator;
import yeswecan.model.RateMatrix;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimHKY extends SimModel {
    
    private Tree tree;
    private TsTvRatioAdvanced kappa;
    private BaseFrequencies freqs;
    private Random rand;
    private AlignmentBuilder siteStates;
    private int length;
    private boolean printSubCounts;
    
    
    public SimHKY(Tree tree, double kappaValue, double[] baseFrequencyValues, 
            Random rand, int length, boolean printSubCounts){
      
        this.tree = tree;
        this.kappa = new TsTvRatioAdvanced(kappaValue);
        this.freqs = new BaseFrequencies(baseFrequencyValues);
        this.rand = rand;
        this.length = length;
        this.printSubCounts = printSubCounts;
    }
    
    public Alignment simulate(){
        Alignment[] sites = new Alignment[this.length];
        
        if (this.printSubCounts){
            System.out.println(super.subCountHeader);
        }
        
        RateMatrix Q = new RateMatrix(this.kappa, this.freqs);
        ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(Q);
        Node root = tree.getRoot();
        
        for (int iSite = 0; iSite < this.length; iSite++) {
            int rootState = SimModel.draw(this.freqs.get(), this.rand.nextDouble());
            
            this.siteStates = new AlignmentBuilder(this.tree.getExternalNodeCount()); // an 'alignment' for a single site, which will be populated with states by downTree.

            SubCount count = new SubCount();
            
            downTree(this.tree, Pgen, root, rootState, count, this.siteStates, rand);
            
            if (this.printSubCounts){
                System.out.println(iSite + Constants.OUTPUT_DELIMITER + "0" + Constants.OUTPUT_DELIMITER + iSite%3 + Constants.OUTPUT_DELIMITER + count.count); // hard coded 0 represents the alignment partition. Since there are no partitiions with HKY, every site is in partition 0
            }
            
            sites[iSite] = siteStates.generateAlignment(new Nucleotides());

        } // for iSite
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        return allSites;
    
    }// simulate

    
}
