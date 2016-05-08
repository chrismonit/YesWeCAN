/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.sim;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import pal.tree.Tree;
import yeswecan.Constants;
import yeswecan.cli.CommandArgs;
import yeswecan.model.submodels.CANModelMixture;
import yeswecan.phylo.GeneticStructure;
import yeswecan.run.RunCANMixture;

/**
 *
 * @author cmonit1
 */
public class SimFreqsMix extends SimFreqs {
        
        protected CANModelMixture canMix;
        protected CommandArgs comArgs;
        protected int numSiteClasses;
        
        public SimFreqsMix(Tree tree, Random rand, GeneticStructure genStruct, 
            CommandArgs comArgs){
        
            super(tree, rand, genStruct, comArgs); 
            /* creates instance of FrequencySimulator rather than FrequencySimulatorMix
            Creates list of omegas based on -w argument, rather than -w1 etc  
            
            */          
            
            this.canMix = RunCANMixture.makeMixture(comArgs, comArgs.getModel(), Constants.CODON_FREQ_MIX2_IDENTIFIER);
            this.omegas = canMix.getOmegas(); // probably redundant, for completeness
            System.out.println("can mix num site classes "+canMix.numSiteClasses);

            
            System.exit(1);
            this.comArgs = comArgs;
            
            this.numSiteClasses = numberSiteClasses(comArgs.getModel());

            System.out.println("num site classes "+numSiteClasses);
            
            this.simulator = new FrequencySimulatorMix(
                    tree, rand, genStruct, kappa, this.canMix.getOmegas(),
                    this.codonFrequencies, this.canMix.getProbabilities(),
                    numSiteClasses
            );
            
                       
            
            // override header method
            // override get params method
            
            // print true stite classes thing
            
        }
        
        public void printSiteClasses(){
            //TODO?
        }
        
        
        @Override
        public String[] getHeader(){
            ArrayList<String> columns = new ArrayList<String>();
            Collections.addAll(columns, "model", "kappa", "nuNumRepeats", "nuNumSubs", "equiBranchLength", "allowStops");
            
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
        
    @Override
    public double[] getSimParameters(){
        List<Double> resultList = new ArrayList<Double>();
        
        resultList.add((double)this.comArgs.getModel());
        resultList.add(this.kappa.get());
        resultList.add((double)this.nRepeats);
        resultList.add((double)this.nSubs);
        resultList.add(this.equiBranchLength);
        resultList.add((this.allowStops) ? 1.0 : 0.0);
        
        System.out.println("size "+omegas.size());
        
        for (int iGene = 1; iGene < this.comArgs.getGeneNumber()+1; iGene++) {
            System.out.println("iGene "+iGene);
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                resultList.add( canMix.getOmega(iGene, jClass).get() );
                System.out.println("jClass "+ jClass +" " + canMix.getOmega(iGene, jClass));
            }
            
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                resultList.add( canMix.getProbability(iGene, jClass) );
            }
            
            System.out.println("");
        }
        
        double[] resultArray = new double[resultList.size()];
        for (int i = 0; i < resultList.size(); i++) {
            resultArray[i] = resultList.get(i).doubleValue();
        }
        
        for (swmutsel.model.parameters.Omega w : canMix.getOmegas()){
            System.out.println(w);
            
        }
        
        //System.exit(1);
        return resultArray;
    }
        
        
    public static int numberSiteClasses(int mixtureModel){ // based on near identical method in RunCANMixture
        int numSiteClasses = -1;

        if (mixtureModel == Constants.CODON_FREQ_MIX2_IDENTIFIER)
            numSiteClasses = Constants.NUM_M2_SITE_CLASSES;
        else
            numSiteClasses = Constants.NUM_M1_SITE_CLASSES;
        return numSiteClasses;
    }
        
                
}