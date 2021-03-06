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
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.RateMatrix;
import yeswecan.model.submodels.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.States;

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
    
    
    public SimHKY(Tree tree, Random rand, HKYModel hky, int length, boolean printSubCounts){
      
        this.tree = tree;
        this.kappa = hky.getKappa();
        this.freqs = hky.getPi();
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
            double r = this.rand.nextDouble();
            int rootState = States.draw(this.freqs.get(), r);
            
            this.siteStates = new AlignmentBuilder(this.tree.getExternalNodeCount()); // an 'alignment' for a single site, which will be populated with states by downTree.

            SubCount count = new SubCount();
            downTree(this.tree, Pgen, root, rootState, count, this.siteStates, rand);
            
            if (this.printSubCounts){
                System.out.println(Constants.DEL + iSite + Constants.DEL + "0" + Constants.DEL + iSite%3 + Constants.DEL + count.count); // hard coded 0 represents the alignment partition. Since there are no partitiions with HKY, every site is in partition 0
            }
            
            sites[iSite] = siteStates.generateAlignment(new Nucleotides());

        } // for iSite
        Alignment allSites = new ConcatenatedAlignment(sites);
        return allSites;
    
    }// simulate

    
}
