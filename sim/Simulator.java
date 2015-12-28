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
import pal.alignment.ConcatenatedAlignment;
import pal.datatype.Nucleotides;
import pal.tree.Node;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import yeswecan.model.ProbMatrixFactory;
import yeswecan.model.ProbMatrixGenerator;
import yeswecan.model.RateMatrix;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.FastaWriter;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Simulator {
    
    // first making a simulator for HKY which can be extended
    
    public static void main(String[] args){
        // read in command args etc
        String treePath = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/simulator_tests/basic.tre";
        double kappa = 2.0;
        double[] baseFrequencies = new double[]{.1,.2,.3,.4};
        double branchScaling = 1.0;
        
        String a = "0,1,1,1,2";
        String b = "2,2,0,0,0";
        String c = "3,3,3,3,0";
        String lengths = "50,50,50,50,50";
        
        Simulator sim = new Simulator(treePath, kappa, baseFrequencies, branchScaling);
    }
    
    Random rand = new Random();
    private AlignmentBuilder siteStates;

    
    public Simulator(String treePath, double kappaValue, double[] baseFrequencyValues, double branchScalingValue){
        Tree tree = loadTree(treePath);
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(kappaValue);
        BaseFrequencies freqs = new BaseFrequencies(baseFrequencyValues);
        BranchScaling scaling = new BranchScaling(branchScalingValue);
        int numberSites = 1000000;
        
        Alignment[] sites = new Alignment[numberSites];
        
        for (int iSite = 0; iSite < numberSites; iSite++) {
            
            RateMatrix Q = new RateMatrix(kappa, freqs);
            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(Q);

            Node root = tree.getRoot();

            double r = rand.nextDouble();
            int rootState = draw(freqs.get(), r);

            siteStates = new AlignmentBuilder(10);

            downTree(tree, Pgen, scaling, root, rootState);

            sites[iSite] = siteStates.generateAlignment(new Nucleotides());
        }
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        System.out.println(new FastaWriter().fastaString(allSites));
        //System.out.println(allSites.toString());
    }
    
    
    
    
    public void downTree(Tree tree, ProbMatrixGenerator Pgen, BranchScaling scaling, 
            Node parentNode, int parentState){
        
        if (parentNode.isLeaf()) {
            siteStates.addSequence(new int[]{parentState}, parentNode.getIdentifier().getName());
        }
        else{
            for (int iChild = 0; iChild < parentNode.getChildCount(); iChild++) {
                Node childNode = parentNode.getChild(iChild);
                double rawBranchLength = childNode.getBranchLength();
                RealMatrix P = Pgen.getP( scaling.get() * rawBranchLength);
                double[] transProbDistribution = P.getRow(parentState);
                
                int childState = draw(transProbDistribution, rand.nextDouble());
                downTree(tree, Pgen, scaling, childNode, childState);
            }// for iChild
        }// else (not leaf)
        
    }// downTree
    
    
    
    // draw a state (0/1/2/3 for nt states) from a probability distribution
    public int draw(double[] distribution, double randomUniform){
        // NB elements of distribution MUST sum to one!
        // randomUniform must be between 0 and 1
        double current = 0.0;
        for (int i = 0; i < distribution.length; i++) {
            if (randomUniform > current && randomUniform < current + distribution[i] ) {
                return i;
            }
            else{
                current += distribution[i];
            }
        } //for 
        return distribution.length-1; 
    }
    
    
    
    
    public Tree loadTree(String treePath){
        try{
            return new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println("Failed to read in tree for Simulator");
            e.printStackTrace();
            System.exit(1);
        }
        return null;
    }
    
    
    
}
