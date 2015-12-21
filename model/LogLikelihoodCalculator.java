
package yeswecan.model;

import yeswecan.phylo.AdvancedAlignment;

import pal.tree.Node;
import pal.tree.Tree;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 */
public class LogLikelihoodCalculator {


    public LogLikelihoodCalculator(){}
    
    public static double calculateSiteLogLikelihood(AdvancedAlignment alignment, Tree tree, int site, ProbMatrixGenerator pGenerator){ 
        
        //System.out.println("root:" + tree.getRoot().toString());
        
        double sum = 0.0;
        double[] rootConditionals = downTree( tree.getRoot(), alignment, tree, site, pGenerator );
        //System.out.println("now");
        //ArrayPrinter.print(ReorderFrequencies.alphaToPaml(pGenerator.getQ().getBaseFrequencies().get()), ",");
        for (int iRootState = 0; iRootState < States.NT_STATES; iRootState++) {
            sum +=  pGenerator.getQ().getBaseFrequencies().get()[iRootState] * rootConditionals[iRootState]; //TODO surely we can make this more efficient
        }
        return Math.log(sum);
    }
    
    private static double[] downTree(Node parent, AdvancedAlignment alignment, Tree tree, int site, ProbMatrixGenerator pGenerator){
        //Felsenstein's Pruning Algorithm. Inspired by (virtually copied from) R. Goldstein's downTree method in Play.java
        double[] parentConditionals = new double[States.NT_STATES]; //number of nucleotide states
        //System.out.println("node id: " + parent.getIdentifier());
        
        if (parent.isLeaf()){ // 'parent' is terminal node, i.e. has no children.
        
            String taxonName = parent.getIdentifier().getName();
            int state = alignment.getStateBySequenceName(taxonName, site);
            //System.out.println("site: " + site + "\t" + "taxonName: " + taxonName + "\t" + "state: " + state + "\t");
        
            if (state >= 0 && state < States.NT_STATES){ //the observed state is recognised as a nucleotide
                //System.out.println("state: " + state);
                parentConditionals[state] = 1.0;
            }
            else  {
                for (int i = 0; i < parentConditionals.length; i++) {
                    parentConditionals[i] = 1.0; //observed state is not recognised as nucleotide (may be gap). Treated as missing data. All conditional probabilities = 1.0.
                    //System.out.println("ELSE state: " + state);
                }
            }
        }//if leaf
        else{
            //parentConditionals dependent on children conditionals
            for (int i = 0; i < parentConditionals.length; i++) 
                parentConditionals[i] = 1.0; // set all to one at first, will multiply by probabilities for each child, below
            //System.out.println(parent.getNumber() + "\t" + parent.getChildCount() + "\t" + parent.getBranchLength());
            
            
            for (int iChild = 0; iChild < parent.getChildCount(); iChild++){ //for every child of this parent node
                Node child = parent.getChild(iChild);
                
                double t = child.getBranchLength();
                
                double[][] P_t = pGenerator.getP(t).getData(); //don't need double[][], could use RealMatrix getElement method
                
                //MatrixPrinter.PrintMatrix(P_t, "P t="+t);
                
                double[] childConditionals = downTree(child, alignment, tree, site, pGenerator);
                
                for (int iParentState = 0; iParentState < States.NT_STATES; iParentState++) { //iParentState == x_i in Yang (2006) equation 4.4
                    double nodeConditionalLikelihood = 0.0; //prob of observing data below this node, if the state at this node were iParentState (i.e. 'conditional' on state being iParentState). This is L_i(x_i) in Yang (2006) eq 4.4
                    for (int jChildState = 0; jChildState < States.NT_STATES; jChildState++) {
                        double p_ij = P_t[iParentState][jChildState];
                        nodeConditionalLikelihood += p_ij * childConditionals[jChildState];
                    }// jChildState
                    parentConditionals[iParentState] *= nodeConditionalLikelihood;
                }// iParentState
            }//for iChild
            
//     
        }//else (if not leaf)
        
        //ArrayPrinter.print(ReorderFrequencies.alphaToPaml(parentConditionals), ",");
        //System.out.println(parent.toString());
        //System.out.println("");
        
        return parentConditionals;
    
    }//downTree()

}//class