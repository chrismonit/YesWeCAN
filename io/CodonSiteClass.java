/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public abstract class CodonSiteClass {

    protected abstract String getHeader();
    
    protected abstract List<String> getGeneRows();
    
    
    public void print(){
        System.out.println( getHeader() );
        for (String row : getGeneRows()){
            System.out.println(row);
        }
    }
    
    public void write(String filePath){
    
        try{
            FileWriter fileWriter = new FileWriter(filePath);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            
            bufferedWriter.write( getHeader() );
            bufferedWriter.newLine();
            for (String row : getGeneRows()){
                bufferedWriter.write(row);
                bufferedWriter.newLine();
            }

            bufferedWriter.close();
        
        }catch (IOException e){
            e.printStackTrace();
        }
    }
    
    
    // static utility methods
    
    public static int[] intArray(List<Integer> intList){
        int[] array = new int[intList.size()];
        for (int i = 0; i < intList.size(); i++) {
            array[i] = intList.get(i);
        }
        return array;
    }
    
    public static boolean sitesSequential(int[] codon){
        return (codon[1] == codon[0]+1 && codon[2] == codon[0]+2);
    }
    
    public static int[] getGeneContiguousSites(int gene, GeneticStructure genStruct){
        List<Integer> geneSites = new ArrayList<Integer>();
        
        for (int iSite = 0; iSite < genStruct.getTotalLength(); iSite++) {
            int partition = genStruct.getPartitionIndex(iSite);
            
            if (genStruct.containsGene(partition, gene)) {
                geneSites.add(iSite);
            }
        }//for iSite
        
        return intArray(geneSites);
    }
    
    // if array length is not multiple of three, will just ignore the last 1 or 2 elements
    public static int[][] getCodons(int[] sites){
        
        int nCodons = sites.length/3; // integer division. Only considers number of complete triplets
        
        int[][] codons = new int[nCodons][3];
        
        for (int iSite = 0; iSite < sites.length-2; iSite+=3) {
            int codonIndex = iSite/3;
            codons[codonIndex] = new int[]{ sites[iSite], sites[iSite+1], sites[iSite+2] };
        }
        return codons;
    }
    
}
