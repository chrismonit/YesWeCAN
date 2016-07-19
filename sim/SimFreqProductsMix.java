/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;

import java.util.ArrayList;
import java.util.Random;
import pal.alignment.Alignment;
import pal.alignment.AlignmentBuilder;
import pal.alignment.ConcatenatedAlignment;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.datatype.Nucleotides;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
import yeswecan.io.CommandArgs;
import yeswecan.model.matrices.CodonAwareMatrix;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.codonawareness.RatioScaler;
import yeswecan.model.codonawareness.RatioScalerFactory;
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.model.submodels.CANModelMixture;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimFreqProductsMix extends SimModel { // TODO this is just SimCANMixture class copied with the name changed!!!
    
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModelFrequenciesMix canMix;
    private Random rand;
    private AlignmentBuilder siteStates;
    private boolean printSubCounts;
    private CommandArgs comArgs;
    private CodonFrequencies[] codonFrequenciesArray;
    private TsTvRatioAdvanced kappa;
    private CodonTable universalTable;
    
    
    public SimFreqProductsMix(Tree tree, Random rand, 
            GeneticStructure genStruct, CommandArgs comArgs, boolean printSubCounts){
      
        this.tree = tree;
        this.genStruct = genStruct;
        this.rand = rand;
        
        this.printSubCounts = printSubCounts;

        this.comArgs = comArgs;
        
        
        this.kappa = new TsTvRatioAdvanced(this.comArgs.kappa());
        
        this.codonFrequenciesArray = new CodonFrequencies[this.genStruct.getNumberOfGenes()+1]; // plus 1 for no gene case
        codonFrequenciesArray[0] = new CodonFrequencies(); // default constructor has all freq = 1/64 for no gene case

        CodonFrequencies geneFrequencies = new CodonFrequencies(comArgs.getCodonFrequencyPath()); // only allowing one set of frequencies for all genes for now at least
        for (int i = 0; i < this.genStruct.getNumberOfGenes(); i++) {
            this.codonFrequenciesArray[i+1] = geneFrequencies; // use reference to same instance for all genes
        }
            
        this.universalTable = CodonTableFactory.createUniversalTranslator();
        System.out.println("MODEL"+Constants.DEL+"site"+Constants.DEL+"A"+Constants.DEL+"B"+Constants.DEL+"C"); // just a header for the class and w rows

    }

    
    
    
    public Alignment simulate(){
        Alignment[] sites = new Alignment[this.genStruct.getTotalLength()];
                
        Node root = tree.getRoot();
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            // determine evolutionary process (from site type and partition)
        
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
        
            int aFrameSiteClass = States.draw(this.canMix.getGeneProbabilities(genes[0]).get(), rand.nextDouble());
            int bFrameSiteClass = States.draw(this.canMix.getGeneProbabilities(genes[1]).get(), rand.nextDouble());
            int cFrameSiteClass = States.draw(this.canMix.getGeneProbabilities(genes[2]).get(), rand.nextDouble());
                        
            Omega aOmega = this.canMix.getGeneAndSiteClassOmega(genes[0], aFrameSiteClass); // genes[0] is the gene present in frame A
            Omega bOmega = this.canMix.getGeneAndSiteClassOmega(genes[1], bFrameSiteClass);  
            Omega cOmega = this.canMix.getGeneAndSiteClassOmega(genes[2], cFrameSiteClass);  
            
            Omega[] omegas = new Omega[]{ aOmega, bOmega, cOmega };
            
            System.out.println("class"+Constants.DEL+iSite+Constants.DEL+aFrameSiteClass+Constants.DEL+bFrameSiteClass+Constants.DEL+cFrameSiteClass);
            System.out.println("w"+Constants.DEL+iSite+Constants.DEL+aOmega.get()+Constants.DEL+bOmega.get()+Constants.DEL+cOmega.get());
            
            CANMatrixFreqProducts canQ = new CANMatrixFreqProducts(this.kappa, siteType, omegas, this.codonFrequenciesArray, this.universalTable);
            
            ProbMatrixGenerator Pgen = ProbMatrixFactory.getPGenerator(canQ);
            
            // simulate according to process
            
            int rootState = States.draw(canQ.getBaseFrequencies().get(), rand.nextDouble());

            this.siteStates = new AlignmentBuilder(this.tree.getExternalNodeCount()); // an 'alignment' for a single site, which will be populated with states by downTree.
            
            SubCount count = new SubCount();
            
            SimModel.downTree(tree, Pgen, root, rootState, count, this.siteStates, this.rand);
            
            if (this.printSubCounts){
                System.out.println(super.subCountHeader);
                System.out.println("count"+Constants.DEL +iSite + Constants.DEL + genStruct.getPartitionIndex(iSite) + Constants.DEL + iSite%3 + Constants.DEL + count.count); 
            }
            System.out.println("");

            // add this newly simulated site to the total set
            sites[iSite] = siteStates.generateAlignment(new Nucleotides());
                    
        }// for iSite
        
        Alignment allSites = new ConcatenatedAlignment(sites);
        return allSites;
    }
     
}
