package yeswecan.sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import pal.alignment.Alignment;
import pal.alignment.AlignmentBuilder;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.datatype.Codons;
import pal.datatype.Nucleotides;
import pal.tree.Node;
import pal.tree.ReadTree;
import pal.tree.SimpleNode;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class FrequencySimulator {
    
    private Tree tree;
    private Random rand;
    private GeneticStructure genStruct;
    private List<Omega> omegas;
    private TsTvRatioAdvanced kappa;
    private List<CodonFrequencies> codonFrequencies;
    private CodonTable codonTable;
    
    private int nTERMINAL_SITES_IGNORED = 2;
    
    
    public FrequencySimulator(Tree tree, String alignmentDestinationPath,
        Random rand, GeneticStructure genStruct, TsTvRatioAdvanced kappa,
        List<Omega> omegas, List<CodonFrequencies> codonFrequencies){
    
        this.tree = tree;
        this.rand = rand;
        this.genStruct = genStruct;
        this.omegas = omegas;
        this.kappa = kappa;
        this.codonFrequencies = codonFrequencies;
        this.codonTable = CodonTableFactory.createUniversalTranslator();
    }
    
    // TODO make this private after testing
    public double computeRate(int[] quintStates, int j, int site){
        // r_ijl = k * w_A * w_B * w_C * π_A * π_B * π_C

        int[] genes = this.genStruct.getGenes(site);
        double product = 1.0; // multiplicative identity
                
        product *= this.kappa.getKappaIfTransition(quintStates[2], j);
        for (int iFrame = 0; iFrame < 3; iFrame++) {
            int[] codonI_array = getCodon(quintStates, quintStates[2], iFrame, site%3);
            int codonI = Codons.getCodonIndexFromNucleotideStates(codonI_array);
            
            int[] codonJ_array = getCodon(quintStates, j, iFrame, site%3);

            int codonJ = Codons.getCodonIndexFromNucleotideStates(codonJ_array);

            if (!this.codonTable.isSynonymous(codonI, codonJ)) { 
                product *= this.omegas.get(genes[iFrame]).get();
            }
            
            CodonFrequencies geneCodonFreq = this.codonFrequencies.get(genes[iFrame]);
            int[] mappedToPaml = ReorderFrequencies.alphaToPaml(codonJ_array); // expecting codons will be ordered TCAG in codonFrequencies instances
            double pi_J = geneCodonFreq.getFrequency(mappedToPaml); 
            product *= pi_J;
            
        }// iFrame
        return product;
    }
    
    
    
    public double computeSumRates(int[] sequence){
        double sum = 0.0;
        
        for (int iSite = 2; iSite < sequence.length-2; iSite++) { // can't include first and last 2 nuceltodides, because we're working with quints
            int[] quint = new int[]{ sequence[iSite-2], sequence[iSite-1], sequence[iSite], sequence[iSite+1], sequence[iSite+2] };
            
            for (int jMutation : States.getMutationStates(sequence[iSite], 4)){ // only considering iStopCodon != j
                sum += computeRate(quint, jMutation, iSite);
            }
        }
        
        return sum;
    }
    
    public int[] substitute(int[] sequence){
        int[] newSequence = new int[sequence.length];
        System.arraycopy(sequence, 0, newSequence, 0, newSequence.length);
        
        double R = computeSumRates(newSequence);
        // Compute mutation probabilities and make substitution
        /*
            The method States.draw takes a distribution array and returns the index of the element chosen
            If I make a probabilities array, containing probs for all iStopCodon-> transitions at each site, 
            I will lose the information about the exact site and mutatino which is chosen by States.draw.
            I therefore have a Hashtable, with probabilities[] indices as keys and (site, mutation) pairs as values      
        */
        // chose a siteAndMutationIndex proportional to r_ijl
        Hashtable<Integer, List<Integer>> mutationStatesTable = new Hashtable(); // keys are indices in probabilities array, values are (site, siteAndMutationIndex) pairs
        int probabilitiesLength = (newSequence.length - 2*nTERMINAL_SITES_IGNORED) * (States.NT_STATES-1); // -4 because don't include first and last 2 bases; *3 because 3 siteAndMutationIndex states for each site
        double[] probabilities = new double[probabilitiesLength]; 
        for (int iSite = nTERMINAL_SITES_IGNORED; iSite < (newSequence.length - nTERMINAL_SITES_IGNORED); iSite++) {

            int[] quint = new int[]{ newSequence[iSite-2], newSequence[iSite-1], newSequence[iSite], newSequence[iSite+1], newSequence[iSite+2] };

            int[] mutationStates = States.getMutationStates(newSequence[iSite], 4); // all j != iStopCodon
            for (int iMutationState = 0; iMutationState < mutationStates.length; iMutationState++) {
                //System.out.println("Mutation "+mutationStates[iMutationState]);
                int probsIndex = ((iSite-nTERMINAL_SITES_IGNORED)*(States.NT_STATES-1))+iMutationState; // == (iGene-2) * 3 + iMutationState
                probabilities[probsIndex] = computeRate(quint, mutationStates[iMutationState], iSite) / R; // remember R = \sum_{\ell} \sum_{j \neq iStopCodon} r_{ij\ell}

                List<Integer> siteAndMutation = new ArrayList<Integer>();
                siteAndMutation.add(iSite);
                siteAndMutation.add(mutationStates[iMutationState]);
                mutationStatesTable.put(probsIndex, siteAndMutation);

            }// for mutation

        }// for iGene

        // sanity check
        double sumProbs = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            sumProbs += probabilities[i];
        }
        if (sumProbs < 1. - 1e-10 || sumProbs > 1. + 1e-10){
            throw new RuntimeException("Warning! probabilities don't sum to 1! Sum = "+sumProbs);
        }

        int siteAndMutationIndex = States.draw(probabilities, this.rand.nextDouble());
        List<Integer> siteMutationPair = mutationStatesTable.get(siteAndMutationIndex);
        int site = siteMutationPair.get(0);
        int mutationState = siteMutationPair.get(1);

        newSequence[site] = mutationState; // mutation is now a substitution (ie accepted)
        //ArrayPrinter.print(evolvingSequence, ",");
            
        return newSequence;
    }// substitute        
    
    
    
    
    
    
    // branchPosition is current position along the branch as me move along it, branch length is the total length
    public int[] evolveBranch(int[] startSequence, double branchLength){
        int[] evolvingSequence = new int[startSequence.length];
        System.arraycopy(startSequence, 0, evolvingSequence, 0, evolvingSequence.length);
            //System.out.println("branch length "+branchLength);
        double branchPosition = 0.0;
        while (branchPosition < branchLength){     
            
           
            //System.out.println("branch position: "+branchPosition);

            double R = computeSumRates(evolvingSequence);
            //System.out.println("R "+R);
            ExponentialDistribution expDist = new ExponentialDistribution(1./R); // expecting the MEAN of the distribtuion, which is inverse of rate parameter (1/lambda)
            double deltaT = expDist.sample();
            branchPosition += deltaT;
            //System.out.println("deltaT "+deltaT);
            //System.out.println("branchPosition "+branchPosition);
            if (branchPosition >= branchLength){
                break; // we're very close to the end of the branch, so don't evolve further
            }
            
            evolvingSequence = substitute(evolvingSequence);
////            
////            // Compute mutation probabilities
////            
////            /*
////                The method States.draw takes a distribution array and returns the index of the element chosen
////                If I make a probabilities array, containing probs for all iStopCodon-> transitions at each site, 
////                I will lose the information about the exact site and mutatino which is chosen by States.draw.
////                I therefore have a Hashtable, with probabilities[] indices as keys and (site, mutation) pairs as values      
////            */
////            // chose a siteAndMutationIndex proportional to r_ijl
////            Hashtable<Integer, List<Integer>> mutationStatesTable = new Hashtable(); // keys are indices in probabilities array, values are (site, siteAndMutationIndex) pairs
////            int probabilitiesLength = (evolvingSequence.length - 2*nTERMINAL_SITES_IGNORED) * (States.NT_STATES-1); // -4 because don't include first and last 2 bases; *3 because 3 siteAndMutationIndex states for each site
////            double[] probabilities = new double[probabilitiesLength]; 
////            for (int iSite = nTERMINAL_SITES_IGNORED; iSite < (evolvingSequence.length - nTERMINAL_SITES_IGNORED); iSite++) {
////                
////                int[] quint = new int[]{ evolvingSequence[iSite-2], evolvingSequence[iSite-1], evolvingSequence[iSite], evolvingSequence[iSite+1], evolvingSequence[iSite+2] };
////                
////                int[] mutationStates = States.getMutationStates(evolvingSequence[iSite], 4); // all j != iStopCodon
////                for (int iMutationState = 0; iMutationState < mutationStates.length; iMutationState++) {
////                    //System.out.println("Mutation "+mutationStates[iMutationState]);
////                    int probsIndex = ((iSite-nTERMINAL_SITES_IGNORED)*(States.NT_STATES-1))+iMutationState; // == (iGene-2) * 3 + iMutationState
////                    probabilities[probsIndex] = computeRate(quint, mutationStates[iMutationState], iSite) / R; // remember R = \sum_{\ell} \sum_{j \neq iStopCodon} r_{ij\ell}
////                    
////                    List<Integer> siteAndMutation = new ArrayList<Integer>();
////                    siteAndMutation.add(iSite);
////                    siteAndMutation.add(mutationStates[iMutationState]);
////                    mutationStatesTable.put(probsIndex, siteAndMutation);
////                    
////                }// for mutation
////                
////            }// for iGene
////            
////            //ArrayPrinter.print(probabilities, ",");
////
////            // sanity check
////            double sumProbs = 0.0;
////            for (int i = 0; i < probabilities.length; i++) {
////                sumProbs += probabilities[i];
////            }
////            if (sumProbs < 1. - 1e-10 || sumProbs > 1. + 1e-10){
////                throw new RuntimeException("Warning! probabilities don't sum to 1! Sum = "+sumProbs);
////            }
////
////            
////            int siteAndMutationIndex = States.draw(probabilities, this.rand.nextDouble());
////            List<Integer> siteMutationPair = mutationStatesTable.get(siteAndMutationIndex);
////            int site = siteMutationPair.get(0);
////            int mutationState = siteMutationPair.get(1);
////                        
////            evolvingSequence[site] = mutationState; // mutation is now a substitution (ie accepted)
////            //ArrayPrinter.print(evolvingSequence, ",");
        }        
        
        return evolvingSequence;
    }
    
    
    private int[] bytesToInts(byte[] bytes){
        int[] ints = new int[bytes.length];
        for (int i = 0; i < ints.length; i++) {
            ints[i] = (int)bytes[i];
        }
        return ints;
    }
    
    private byte[] intsToBytes(int[] ints){
        byte[] bytes = new byte[ints.length];
        for (int i = 0; i < ints.length; i++) {
            bytes[i] = (byte)ints[i];
        }
        return bytes;
    }
    // UNTESTED
    private void downTree(Node parent, AlignmentBuilder alnBuilder){
        int[] parentSequence = bytesToInts(parent.getSequence());
        
        if (parent.isLeaf()){
            alnBuilder.addSequence(parentSequence, parent.getIdentifier().getName());
        }else{
            for (int iChild = 0; iChild < parent.getChildCount(); iChild++) {
                Node child = parent.getChild(iChild);
                double branchLength = child.getBranchLength();
                
                int[] childSequence = evolveBranch(parentSequence, branchLength);
                child.setSequence(intsToBytes(childSequence));
                
                downTree(child, alnBuilder);
            }
        }
        
    }

    /**
     * We can generate a random sequence and then evolve it for a while according to this model
     * prior to starting the simulation properly, so that the actual starting sequence 
     * is consistent with the process we are describing
     */
    
    public int[] getRandomSequence(int sequenceLength){
    
        int[] sequence = new int[sequenceLength];
        
        for (int i = 0; i < sequence.length; i++) {
            sequence[i] = this.rand.nextInt(States.NT_STATES);
        }
        return sequence;
    }
    

    
    // assumes the gene sequence is in the frame of interest
    public boolean containsStop(int[] geneSequence){
        
        //System.out.println("containsStop "+ArrayPrinter.toString(geneSequence, ","));
        
        int[] stopCodons = codonTable.getTerminatorIndexes(); // will return 48, 50, 56
        int upperBound = geneSequence.length-2;
        
        // NB this will just ignore terminal nucletodies which don't belong to a whole codon
        for (int iSite = 0; iSite < upperBound; iSite+=3) {
            int[] codonArray = new int[]{ geneSequence[iSite], geneSequence[iSite+1], geneSequence[iSite+2]};
            int codon = Codons.getCodonIndexFromNucleotideStates(codonArray);
            //System.out.println("containsStop: codon = "+ArrayPrinter.toString(codonArray, ","));
            
            for (int iStopCodon = 0; iStopCodon < stopCodons.length; iStopCodon++) {
                if (codon == stopCodons[iStopCodon]){
                    //System.out.println("codon "+ArrayPrinter.toString(codonArray, ",")+" gene_site "+iSite);
                    return true;
                }
            }
            
        }// for iSite
        return false;
    }
    
    // crude way of dealing with stop codons in random sequences
    // just go through a sequence looking for stop codons and replacing the first T with C
    // nb this can't introduce stops into another frame because C never occurs in a stop codon
    
    int[][] indexAdjustments = new int[][]{
        { 3, 5, 4 },
        { 4, 3, 5 },
        { 5, 4, 3 }
    };
    
    
    public int[] changeStop(int[] geneSequence){
        Random rand = new Random();
        int[] newSequence = new int[geneSequence.length];
        System.arraycopy(geneSequence, 0, newSequence, 0, newSequence.length);
        
        int[] stopCodons = codonTable.getTerminatorIndexes(); // will return 48, 50, 56
        
        int[] adjustments = indexAdjustments[geneSequence.length%3];
        
        for (int iFrame = 0; iFrame < 3; iFrame++) {

            // NB this will just ignore terminal nucletodies which don't belong to a whole codon
            for (int iSite = iFrame; iSite <= (geneSequence.length - adjustments[iFrame]); iSite+=3) {

                int[] codonArray = new int[]{ newSequence[iSite], newSequence[iSite+1], newSequence[iSite+2]};
                int codon = Codons.getCodonIndexFromNucleotideStates(codonArray);

                for (int iStopCodon = 0; iStopCodon < stopCodons.length; iStopCodon++) {
                    if (codon == stopCodons[iStopCodon]){
                        //newSequence[iSite] = States.getMutationStates(evolvingSequence[iSite], States.NT_STATES)[rand.nextInt(States.NT_STATES)]; // choses one of the mutation states at random
                        newSequence[iSite] = 1; // replace with C
                    }
                }

            }// for iSite

        }// jFrame
        return newSequence;
        
    }
    
    /*
    We want to make sure the starting (root) sequence used is viable
    We define viable as simply not containing any stop codons in coding frames
    This method scans a sequence for stop codons
    */
    public boolean sequenceAcceptable(int[] sequence){
        // initialise arrays. Each ArrayList<Integer> represents the sites comprising a gene
        ArrayList<ArrayList<Integer>> genes = new ArrayList<ArrayList<Integer>>();
        for (int iGene = 0; iGene < genStruct.getNumberOfGenes()+1; iGene++) { // let's include noncoding region as a gene
            genes.add( new ArrayList<Integer>() );
        } 
                
        // go through sites and assign sites from the sequence to lists of sites representing specific genes
        for (int iSite = 0; iSite < sequence.length; iSite++) {
            for (int iFrame = 0; iFrame < 3; iFrame++) {
                
                int gene = this.genStruct.getGenes(iSite)[iFrame]; // find which gene we're in, given this site and frame
                genes.get(gene).add( sequence[iSite] );

                //System.out.println("site, frame, gene: " + iSite + "\t" + iFrame + "\t" + gene);
                
            } // iFrame
        }// iSite
        
        // go through genes and check for stop codons
        
        for (int iGene = 1; iGene < genStruct.getNumberOfGenes()+1; iGene++) { // don't include noncoding regions (gene 0)
            // first need to convert Integer[] to int[]
            Integer[] integerGeneSeq = new Integer[genes.get(iGene).size()];
            integerGeneSeq = genes.get(iGene).toArray(integerGeneSeq);
            int[] geneSeq = new int[integerGeneSeq.length];
            for (int i = 0; i < geneSeq.length; i++) {
                geneSeq[i] = integerGeneSeq[i].intValue();
            }
            
            if (containsStop(geneSeq)){
                return false;
            }
            
        }// for iGene
        return true;
    }// sequenceAcceptable
    
    public Alignment simulate(int[] rootSequence){
        AlignmentBuilder alnBuilder = new AlignmentBuilder(this.tree.getExternalNodeCount());
        
        Node root = new SimpleNode(this.tree.getRoot());
        root.setSequence(intsToBytes(rootSequence));
        downTree(root, alnBuilder);
        
        return alnBuilder.generateAlignment(new Nucleotides());
    }
    
     /*Matrix A is accessed using the site type and the frame you want the terminators for.
     * Each row off matrix B contains the co-ordinates of the nt positions for the terminators you want, relative to position
     * Using the two in combination, we get the characters in the parentSequence which
     * correspond to the terminators we want for each frame
     */

    private static final int[][] A = {  //i = site type (alpha, beta, gamma), jIndex = frame (a, b, c)
        { 2, 0, 1 },
        { 1, 2, 0 },
        { 0, 1, 2 }
    };

    private static final int[][] B = {   //represent nt positions in CODON in the quintuplet, relative to 'position' variable
        { -2, -1, 0 }, //1st terminators in quint
        { -1, 0, 1 }, //2nd terminators in quint
        { 0, 1, 2 }  //3rd terminators in quint
    };
    
    
    /*
        provide array of length 5 with nuceltodie states, with site of interest in the centre (index=2)
        centralState is the nuceltide state at quint[2]. 
            If the quint is not mutated: centralState == quint[2]
            If quint is mutated: centralState == jIndex
        siteType is siteIndex%3 (and obviously not quint[2]%3)
    */
    
    private static int[] getCodon(int[] quint, int centralState, int frame, int siteType){
        int[] quintCopy = new int[quint.length]; // we create a copy so we can change the central state if needed
        System.arraycopy(quint, 0, quintCopy, 0, quint.length);
        quintCopy[2] = centralState; // if centralState == quint[2], then quintCopy == quint
        
        int[] codon = new int[3];
        codon[0] = quintCopy[ 2 + B[ A[siteType][frame] ][0] ];
        codon[1] = quintCopy[ 2 + B[ A[siteType][frame] ][1] ];
        codon[2] = quintCopy[ 2 + B[ A[siteType][frame] ][2] ];
        return codon;
    }
    
    
    public static Tree loadTree(String treePath){  
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
