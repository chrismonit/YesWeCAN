
package yeswecan.model;

//for testing:
import pal.tree.ReadTree;
import pal.alignment.AlignmentReaders;
import pal.alignment.*;
import java.io.FileReader;
import pal.datatype.Nucleotides;
import yeswecan.phylo.AdvancedAlignment;

import pal.tree.Node;
import pal.tree.Tree;
import yeswecan.phylo.States;
import java.lang.Math;

/**
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * STILL NEEDS THOROUGH TESTING
 * 
 */
public class LogLikelihoodCalculator {
    //testing this requires a lot of annoying IO code.
    //Might be better to have a separate testing class for testing LC
    
    //real fields
    
    private AdvancedAlignment alignment;
    private Tree tree;
    //prvate double[] rateMatrix;
    private int site;
    
    //test fields
    
    double[] probMatrix = 
        new double[]{ 1.0, 0.0, 0.0, 0.0,
                      0.0, 1.0, 0.0, 0.0,              
                      0.0, 0.0, 1.0, 0.0,              
                      0.0, 0.0, 0.0, 1.0              
                    };
    
    double[] freqs = new double[]{ 0.25, 0.25, 0.25, 0.25 }; // A is only permissible state
    
    // to be removed
    public LogLikelihoodCalculator(AdvancedAlignment alignment, Tree tree, int site){
        this.alignment = alignment;
        this.tree = tree;
        //this.rateMatrix = rateMatrix
        this.site = site;
    }
    
    public LogLikelihoodCalculator(){
    
    }

    
    
    public double calculateSiteLogLikelihood(AdvancedAlignment alignment, Tree tree, int site, RateMatrix Q){ 
  
        double sum = 0.0;
        double[] rootConditionals = downTree( tree.getRoot(), alignment, tree, site, Q );
        for (int iRootState = 0; iRootState < States.NT_STATES; iRootState++) {
            sum +=  freqs[iRootState] * rootConditionals[iRootState];
            
        }
        return sum;
    }
    
    private double[] downTree(Node parent, AdvancedAlignment alignment, Tree tree, int site, RateMatrix Q){
        //Felsenstein's Pruning Algorithm. Inspired by (virtually copied from) R. Goldstein's downTree method in Play.java
        double[] parentConditionals = new double[States.NT_STATES]; //number of nucleotide states
        
        if (parent.isLeaf()){ // 'parent' is terminal node, i.e. has no children.
        
            String taxonName = parent.getIdentifier().getName();
            int state = alignment.getStateBySequenceName(taxonName, site);
            System.out.println("taxon Name = " + taxonName + "\t" + "state = " + state);
        
            if (state >= 0 && state < States.NT_STATES) //the observed state is recognised as a nucleotide
                parentConditionals[state] = 1.0;
            else  
                for (int i = 0; i < parentConditionals.length; i++) 
                    parentConditionals[i] = 1.0; //observed state is not recognised as nucleotide (may be gap). Treated as missing data. All conditional probabilities = 1.0.
            
        }//if leaf
        else{
            //parentConditionals dependent on children conditionals
            for (int i = 0; i < parentConditionals.length; i++) parentConditionals[i] = 1.0; // set all to one at first, will multiply by probabilities for each child, below

            for (int iChild = 0; iChild < parent.getChildCount(); iChild++){ //for every child of this parent node
                Node child = parent.getChild(iChild);
                
                double branchLength = child.getBranchLength();
                ProbMatrixGenerator pGenerator = new ProbMatrixGenerator(Q);
                double[][] P_t = pGenerator.getP(branchLength).getData();
                
                double[] childConditionals = downTree(child, alignment, tree, site, Q);
                
                for (int iParentState = 0; iParentState < States.NT_STATES; iParentState++) { //iParentState == x_i in Yang (2006) equation 4.4
                    double nodeConditionalLikelihood = 0.0; //prob of observing data below this node, if the state at this node were iParentState (i.e. 'conditional' on state being iParentState). This is L_i(x_i) in Yang (2006) eq 4.4
                    for (int jChildState = 0; jChildState < States.NT_STATES; jChildState++) {
                        double p_ij = P_t[iParentState][jChildState];
                        nodeConditionalLikelihood += p_ij * childConditionals[jChildState];
                    }// jChildState
                    parentConditionals[iParentState] *= nodeConditionalLikelihood;
                }// iParentState
                
            }//for iChild
            
        }//else (if not leaf)

        return parentConditionals;
    
    }//downTree()
    
    
    
    
    public static void main(String[] args){
//        System.out.println("Hello World, this is Likelihood Calculator");
//
//        String testTreePath = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/test.tre";
//        String testAlignPath = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/test.fasta";
//        
//        ReadTree testTree;
//        Alignment align;
//        SimpleAlignment simple;
//        AdvancedAlignment advanced;
//
//        try{
//            testTree = new ReadTree(testTreePath);
//
//            align = AlignmentReaders.readFastaSequences(new FileReader(testAlignPath), new Nucleotides());
//            simple = new SimpleAlignment(align);
//            advanced = new AdvancedAlignment(simple);
//            
//        }
//        catch(Exception e){ testTree=null; simple = null; advanced = null; System.out.println("Error reading alignment in LC main"); e.printStackTrace(); System.exit(1);  }
//
//        
////        String[] taxa = new String[]{ "human", "chimp", "gorilla", "orangutan" };
////        for (String s : taxa){
////            System.out.println("taxon: " + s + "\t" + "state: " + advanced.getStateBySequenceName(s, 0));
////        }
//        
//        LogLikelihoodCalculator LC = new LogLikelihoodCalculator(advanced, testTree, 0);
//        
//        Node root = testTree.getRoot();
//        
//        double lnL = Math.log( LC.calculateSiteLogLikelihood() );
//        
//        System.out.println( "Site lnL = " + lnL );
//        
//        System.out.println("End of test");
    }//main
    
    
}//class