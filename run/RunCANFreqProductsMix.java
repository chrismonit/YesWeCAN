/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.run;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.tree.Tree;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Probabilities;
import yeswecan.Constants;
import yeswecan.io.CommandArgs;
import yeswecan.model.functions.CANFunctionFreqProductsMix;
import yeswecan.model.parameters.OmegaNegative;
import yeswecan.model.parameters.OmegaPositive;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import static yeswecan.run.RunCANMixture.makeMixture;
import static yeswecan.run.RunCANMixture.numberSiteClasses;

/**
 *
 * @author cmonit1
 */
public class RunCANFreqProductsMix extends RunModel {
    
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    private GeneticStructure genStruct;
    
    private CodonTable codonTable;
    private CodonFrequencies[] codonFrequenciesArray;
    
    private int numSiteClasses;
    
    private int model;
    
    public RunCANFreqProductsMix(AdvancedAlignment alignment, Tree tree, CommandArgs input, int model){
        this.comArgs = input;
        this.alignment = alignment;
        this.tree = tree;
        
        this.genStruct = new GeneticStructure(
                this.comArgs.aFrame(),
                this.comArgs.bFrame(),
                this.comArgs.cFrame(),
                this.comArgs.lengths()
        );
    
        this.codonFrequenciesArray = new CodonFrequencies[this.genStruct.getNumberOfGenes()+1];

        this.codonFrequenciesArray[0] = new CodonFrequencies(); // 1/64 for frames where there is no gene

        // use same codonFrequencies reference for all coding genes. Could make gene-specific if desired
        CodonFrequencies codonFrequencies = new CodonFrequencies(input.getCodonFrequencyPath());        
        System.out.println(Constants.CODON_FREQ_PATH+Constants.DEL+input.getCodonFrequencyPath());
        for (int iGene = 0; iGene < this.genStruct.getNumberOfGenes(); iGene++) {
            this.codonFrequenciesArray[iGene+1] = codonFrequencies;
        }
        
        this.codonTable = CodonTableFactory.createUniversalTranslator();
        
        if (alignment != null) { // in Simulate class, I make instances of this class but pass in null alignment and tree
            if (this.alignment.getLength() != this.genStruct.getTotalLength()){
                throw new RuntimeException("Sum of partition lengths (-l) and alignment length differ! Please check that partition lengths are correct.");
            }
        }
                
        this.numSiteClasses = numberSiteClasses(model);
        this.model = model;
        System.out.println("test "+1e6);
    }// constructor
    
    
    @Override
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "model", "lnL", "kappa", "sc");
        for (int iGene = 0; iGene < this.comArgs.getGeneNumber(); iGene++) {
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                columns.add(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.OMEGA_STRING + Integer.toString(jClass));
            }
            
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                columns.add(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.PROB_STRING + Integer.toString(jClass));
            }
        }
        return columns.toArray(new String[columns.size()]);
    }
    
    
    
    public static CANModelFrequenciesMix makeCAN(CommandArgs comArgs, int mixtureModel, 
            int complexModelIdentifier, int numberOfSiteClasses){

        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(comArgs.kappa());
        if (comArgs.fix().contains(Constants.FIX_KAPPA)) {
            kappa.setOptimisable(false);
        }
        
        BranchScaling scaling = new BranchScaling(comArgs.scaling());
        if (comArgs.fix().contains(Constants.FIX_SCALING)){
            scaling.setOptimisable(false);
        }
        
        List<Omega> omegas = new ArrayList<Omega>();
        List<Probabilities> probs = new ArrayList<Probabilities>();
        
        // neutral (for noncoding frames)
        
        Probabilities neutralProbs;
       
        OmegaNegative neutralW_0 = new OmegaNegative(1.0); // should never be used
        neutralW_0.setOptimisable(false);
        omegas.add(neutralW_0);
        
        Omega neutralW_1 = new Omega(1.0);
        neutralW_1.setOptimisable(false);
        omegas.add(neutralW_1);
        
        if (mixtureModel == complexModelIdentifier){
            OmegaPositive neutralW_2 = new OmegaPositive(1.0); 
            neutralW_2.setOptimisable(false);
            omegas.add(neutralW_2);
            
            neutralProbs = new Probabilities(new double[]{ 0.0, 1.0, 0.0 }); // all density on w_1. w_1 == 1.0
        }else{
            neutralProbs = new Probabilities(new double[]{ 0.0, 1.0 });
        }
        
        neutralProbs.setOptimisable(false);
        probs.add(neutralProbs);
        
        
        // for coding frames
                
        for (int iGene = 0; iGene < comArgs.getGeneNumber(); iGene++) {
       
            OmegaNegative geneW_0 = new OmegaNegative(comArgs.omega0()[iGene]);
            // fix if needs fixing
            
            if (comArgs.fix().contains(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.OMEGA_STRING + Constants.SITE_CLASS_0)) //+1 for zero based
               geneW_0.setOptimisable(false);
        
            omegas.add(geneW_0);

            Omega geneW_1 = new Omega(1.0);
            geneW_1.setOptimisable(false); // w_1 always fixed to 1
            omegas.add(geneW_1);
            
            Probabilities geneProbs;

            if (mixtureModel == complexModelIdentifier){
                OmegaPositive geneW_2 = new OmegaPositive(comArgs.omega2()[iGene]); 
                if (comArgs.fix().contains(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.OMEGA_STRING + Constants.SITE_CLASS_2)) //+1 for zero based
                   geneW_2.setOptimisable(false);
                omegas.add(geneW_2);

                geneProbs = new Probabilities(new double[]{ comArgs.prob0()[iGene], comArgs.prob1()[iGene], comArgs.prob2()[iGene]});
            }else{
                geneProbs = new Probabilities(new double[]{ comArgs.prob0()[iGene], comArgs.prob1()[iGene] });
            }
            
            if (comArgs.fix().contains(Integer.toString(iGene+1)+Constants.WITIHIN_FIELD_SEPARATOR+Constants.PROB_STRING)) // zero based
                geneProbs.setOptimisable(false);
                
            probs.add(geneProbs);
           
        } // for iGene
        
        return new CANModelFrequenciesMix(
                kappa, scaling, omegas, probs, numberOfSiteClasses
        );
    }

    
    private double[] getValueArray(CANModelFrequenciesMix canMix){
        List<Double> resultList = new ArrayList<Double>();
        resultList.add((double)this.model);
        resultList.add(canMix.getLnL()); 


        resultList.add(canMix.getKappa().get());

        resultList.add(canMix.getScaling().get());
        
        for (int iGene = 1; iGene < this.comArgs.getGeneNumber()+1; iGene++) {
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                resultList.add( canMix.getGeneAndSiteClassOmega(iGene, jClass).get() );
            }
            
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                resultList.add( canMix.getProbability(iGene, jClass) );
            }
        }
        
        double[] resultArray = new double[resultList.size()];
        for (int i = 0; i < resultList.size(); i++) {
            resultArray[i] = resultList.get(i).doubleValue();
        }
        return resultArray;

    }

    @Override
    public  double[] getInitialValues(){ // NB first element does not contain lnL
        return getValueArray( makeCAN(this.comArgs, this.model, Constants.CODON_FREQ_MIX2_IDENTIFIER, numberSiteClasses(this.model)) );
    }
    
   @Override
    public double[] calculate(){
                
        CANModelFrequenciesMix can = makeCAN(this.comArgs, this.model, Constants.CODON_FREQ_MIX2_IDENTIFIER, numberSiteClasses(this.model));
        double[] optimisableParams = Mapper.getOptimisable(can.getParameters()); // map parameters to optimisation space, so FunctionHKY.value canMix use them
        CANFunctionFreqProductsMix calculator = 
                new CANFunctionFreqProductsMix(this.alignment, this.tree, this.genStruct, 
                        can, this.codonFrequenciesArray, this.codonTable, this.numSiteClasses
                );
        
        double[] resultArray = getInitialValues();
        resultArray[1] = calculator.value(optimisableParams);
        return resultArray;
    }
    
    
    @Override
    public double[] fit(){
        
        CANModelFrequenciesMix can = makeCAN(this.comArgs, this.model, Constants.CODON_FREQ_MIX2_IDENTIFIER, numberSiteClasses(this.model));
        
        CANFunctionFreqProductsMix optFunction = 
                new CANFunctionFreqProductsMix(this.alignment, this.tree, this.genStruct, 
                        can, this.codonFrequenciesArray, this.codonTable, this.numSiteClasses
                );

        Optimise opt = new Optimise();
        CANModelFrequenciesMix result = (CANModelFrequenciesMix)opt.optNMS(optFunction, can);
        
        double[] mles = getValueArray(result);        
        return mles;
        
    }
    
    // TODO refactor. Identical method in SimFreqsMix
    public static int numberSiteClasses(int mixtureModel){ // based on near identical method in RunCANMixture
        int numSiteClasses = -1;

        if (mixtureModel == Constants.CODON_FREQ_MIX2_IDENTIFIER)
            numSiteClasses = Constants.NUM_M2_SITE_CLASSES;
        else
            numSiteClasses = Constants.NUM_M1_SITE_CLASSES;
        return numSiteClasses;
    }

    
}
