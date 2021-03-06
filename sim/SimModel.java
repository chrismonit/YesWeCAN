/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;

import java.util.Random;
import org.apache.commons.math3.linear.RealMatrix;
import pal.alignment.Alignment;
import pal.alignment.AlignmentBuilder;
import pal.tree.Node;
import pal.tree.Tree;
import yeswecan.Constants;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public abstract class SimModel {
    
    abstract Alignment simulate();
    
    protected String subCountHeader = 
            "COUNT"+Constants.DEL+"site"+Constants.DEL+"part"+Constants.DEL+"type"+Constants.DEL+"count";
    
    protected class SubCount{
        int count = 0;
    }
    
    public static void downTree(Tree tree, ProbMatrixGenerator Pgen, 
            Node parentNode, int parentState, SubCount count, AlignmentBuilder siteStates, Random rand){
        
        if (parentNode.isLeaf()) {
            siteStates.addSequence(new int[]{parentState}, parentNode.getIdentifier().getName());
        }
        else{
            for (int iChild = 0; iChild < parentNode.getChildCount(); iChild++) {
                Node childNode = parentNode.getChild(iChild);
                double branchLength = childNode.getBranchLength();
                RealMatrix P = Pgen.getP( branchLength);

                double[] transProbDistribution = P.getRow(parentState);
                int childState = States.draw(transProbDistribution, rand.nextDouble());
                
                if (childState != parentState) {
                    count.count += 1;
                }
                
                downTree(tree, Pgen, childNode, childState, count, siteStates, rand);
            }// for iChild
        }// else (not leaf)
        
    }// downTree
    
    
}
