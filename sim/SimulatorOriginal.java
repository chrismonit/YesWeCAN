/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
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
import swmutsel.model.parameters.Omega;
import yeswecan.cli.CommandArgs;
import yeswecan.model.matrices.CodonAwareMatrix;
import yeswecan.model.ProbMatrixFactory;
import yeswecan.model.ProbMatrixGenerator;
import yeswecan.model.matrices.RateMatrix;
import yeswecan.model.ratioscaling.RatioScaler;
import yeswecan.model.ratioscaling.RatioScalerFactory;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.FastaWriter;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimulatorOriginal {
    
    

    private Tree tree;
    private GeneticStructure genStruct;
    private TsTvRatioAdvanced kappa;
    private BaseFrequencies freqs;
    private Omega[] geneOmegas;
    private BranchScaling scaling;
            
    private Random rand = new Random(123456789);
    private AlignmentBuilder siteStates;

    
    public SimulatorOriginal(Tree tree, GeneticStructure genStruct, double kappaValue, double[] baseFrequencyValues, double[] omegaValues, double branchScalingValue){
        
      
        this.tree = tree;
        this.genStruct = genStruct;
        this.kappa = new TsTvRatioAdvanced(kappaValue);
        this.freqs = new BaseFrequencies(baseFrequencyValues);
        this.scaling = new BranchScaling(branchScalingValue);
        
        // geneOmegas[0] will be for 'no gene' regions, ie where omega=1
        this.geneOmegas = new Omega[omegaValues.length+1];
        this.geneOmegas[0] = new Omega(1.0); // neutral evolution
        for (int i = 1; i < this.geneOmegas.length; i++) { // NB starting at i=1
            this.geneOmegas[i] = new Omega(omegaValues[i-1]);
        }
    }
    
    public Alignment simulate(){
       
        Alignment[] sites = new Alignment[this.genStruct.getTotalLength()];
        
        System.out.println(genStruct.toString());
        
        System.out.println("site\tpart\ttype\tcount");
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            // determine evolutionary process (from site type and partition)
            
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
            Omega aOmega = this.geneOmegas[genes[0]]; // genes[0] is the gene present in frame A
            Omega bOmega = this.geneOmegas[genes[1]]; 
            Omega cOmega = this.geneOmegas[genes[2]]; 
                        
            // make model
            
            RatioScaler ratioScaler = RatioScalerFactory.getRatioScaler();
            CodonAwareMatrix canQ = new CodonAwareMatrix(this.kappa, this.freqs, ratioScaler, siteType, aOmega, bOmega, cOmega, this.scaling);            
            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(canQ);
            
            // simulate according to process
            Node root = tree.getRoot();

            double r = rand.nextDouble();
            int rootState = draw(freqs.get(), r);

            siteStates = new AlignmentBuilder(this.tree.getExternalNodeCount()); // an 'alignment' for a single site, which will be populated with states by downTree.
            
            SubCount count = new SubCount();
            
            downTree(tree, Pgen, root, rootState, count);
            
            System.out.println(iSite + "\t" + genStruct.getPartitionIndex(iSite) + "\t" + iSite%3 + "\t" + count.count);
            
            // add this newly simulated site to the total set
            sites[iSite] = siteStates.generateAlignment(new Nucleotides());
        }
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        //////System.out.println(new FastaWriter().fastaString(allSites));
        return allSites;
    }
    
    
    
    public void downTree(Tree tree, ProbMatrixGenerator Pgen, 
            Node parentNode, int parentState, SubCount count){
        
        if (parentNode.isLeaf()) {
            siteStates.addSequence(new int[]{parentState}, parentNode.getIdentifier().getName());
        }
        else{
            for (int iChild = 0; iChild < parentNode.getChildCount(); iChild++) {
                Node childNode = parentNode.getChild(iChild);
                double branchLength = childNode.getBranchLength();
                RealMatrix P = Pgen.getP( branchLength);
                double[] transProbDistribution = P.getRow(parentState);
                
                int childState = draw(transProbDistribution, rand.nextDouble());
                
                if (childState != parentState) {
                    count.count += 1;
                }
                
                downTree(tree, Pgen, childNode, childState, count);
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
    
    private class SubCount{
        
        int count = 0;
    }
    
    
}
