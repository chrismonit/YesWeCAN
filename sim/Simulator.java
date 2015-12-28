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
import swmutsel.model.parameters.Omega;
import yeswecan.model.CodonAwareMatrix;
import yeswecan.model.ProbMatrixFactory;
import yeswecan.model.ProbMatrixGenerator;
import yeswecan.model.RateMatrix;
import yeswecan.model.RatioScaler;
import yeswecan.model.RatioScalerFactory;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.FastaWriter;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Simulator {
    
    
    public static void main(String[] args){
        // read in command args etc
        String treePath = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/simulator_tests/basic.tre";
        double kappa = 2.0;
        double[] baseFrequencies = new double[]{.1,.2,.3,.4};
        double branchScaling = 1.0;
        
        double[] omegas = new double[]{0.1, 0.2, 0.3};
        String a = "0,1,1,1,2";
        String b = "2,2,0,0,0";
        String c = "3,3,3,3,0";
        String lengths = "5,5,5,5,5";
        
        GeneticStructure genStruct = new GeneticStructure(a,b,c,lengths,",");
        
        Simulator sim = new Simulator(treePath, genStruct, kappa, baseFrequencies, omegas, branchScaling);
        Alignment aln = sim.simulate();
        System.out.println(new FastaWriter().fastaString(aln));
    }
    
    private Tree tree;
    private GeneticStructure genStruct;
    private TsTvRatioAdvanced kappa;
    private BaseFrequencies freqs;
    private Omega[] geneOmegas;
    private BranchScaling scaling;
            
    private Random rand = new Random();
    private AlignmentBuilder siteStates;

    
    public Simulator(String treePath, GeneticStructure genStruct, double kappaValue, double[] baseFrequencyValues, double[] omegaValues, double branchScalingValue){
        
        this.tree = loadTree(treePath);
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
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            // determine process
            
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
            Omega aOmega = this.geneOmegas[genes[0]]; // genes[0] is the gene present in frame A
            Omega bOmega = this.geneOmegas[genes[1]]; 
            Omega cOmega = this.geneOmegas[genes[2]]; 
                        
            // make model
            

            RatioScaler ratioScaler = RatioScalerFactory.getRatioScaler();
            CodonAwareMatrix canQ = new CodonAwareMatrix(this.kappa, this.freqs, ratioScaler, siteType, aOmega, bOmega, cOmega);            
            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(canQ);

            // simulate according to process
            Node root = tree.getRoot();

            double r = rand.nextDouble();
            int rootState = draw(freqs.get(), r);

            siteStates = new AlignmentBuilder(this.tree.getExternalNodeCount()); // an 'alignment' for a single site, which will be populated with states by downTree.

            downTree(tree, Pgen, scaling, root, rootState);

            // add this newly simulated site to the total set
            sites[iSite] = siteStates.generateAlignment(new Nucleotides());
        }
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        //////System.out.println(new FastaWriter().fastaString(allSites));
        return allSites;
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
