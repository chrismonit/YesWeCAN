/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.phylo;

import java.util.List;
import yeswecan.utils.ArrayPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SequenceProfile {
    
    private SiteProfile[] siteProfiles;
    private GeneticStructure genStruct;
    public static final int nTERMINAL_SITES = 2;
    public static final int nBOUNDARY_SITES = 2;
    
    public SequenceProfile(GeneticStructure genStruct){
        this.genStruct = genStruct;
        this.siteProfiles = new SiteProfile[(genStruct.getTotalLength() - nTERMINAL_SITES*2 )]; // ignoring first and last two positions in sequence
        init();
    }
    
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
    
    // suitable only for contiguous qunituplets (ie no good for partition boundaries)
    private static int[] getCodonSiteIndices(int[] quintIndices, int frame, int siteType){ // siteType is the site type of the central site
        int[] codon = new int[3];
        codon[0] = quintIndices[ 2 + B[ A[siteType][frame] ][0] ];
        codon[1] = quintIndices[ 2 + B[ A[siteType][frame] ][1] ];
        codon[2] = quintIndices[ 2 + B[ A[siteType][frame] ][2] ];
        return codon;
    }
    
    // what if there is no other occurrance of this gene?
    // I guess I have to make sure my genes start and end with complete codons
    private int[] getNext(int currentPartition, int nSiteIndices, int gene){
        int[] toReturn = new int[nSiteIndices];
        
        int partitionOfInterest = -1;
        for (int iPartition = (currentPartition+1); iPartition < genStruct.getNumberOfPartitions(); iPartition++) {
            if (genStruct.containsGene(iPartition, gene)){
                partitionOfInterest = iPartition;
                break;
            }
        }
        
        if (partitionOfInterest == -1){
            // freak out
        }
        
        int partitionSites = genStruct.getPartitionStart(partitionOfInterest);
        for (int i = 0; i < toReturn.length; i++) {
            toReturn[i] = partitionSites;
            partitionSites++;
        }
        
        return toReturn;
    }
    
    // please rename this
    private int[] getCompletedCodon(int siteIndex, int frame, int nSiteIndices){
        
        int[] codon;
        
        if (siteIndex%3 == 0){
            switch (frame){
                case 0: codon = new int[]{ siteIndex, siteIndex+1, getNext(  ) };
                    break;
            }
        
        }
        
    }
    
    private void init(){
        
        // deal with first partition, where first bases are invariant
        
        // could have a method for doing this loop, which you provide as arguments the start and end indices
        for (int iSite = nTERMINAL_SITES; iSite < (genStruct.getPartitionLength(0)-nBOUNDARY_SITES); iSite++) {
            int[] quint = new int[]{ iSite-2, iSite-1, iSite, iSite+1, iSite+2 };
            int[] aCodon = getCodonSiteIndices(quint, 0, iSite%3);
            int[] bCodon = getCodonSiteIndices(quint, 1, iSite%3);
            int[] cCodon = getCodonSiteIndices(quint, 2, iSite%3);
            
            siteProfiles[iSite] = new SiteProfile(iSite, aCodon, bCodon, cCodon);
        }
        // deal with sites at end of first partition
        // penultimate
        
        
        
        for (int iPartition = 1; iPartition < (genStruct.getNumberOfPartitions() - 1); iPartition++) {
            
            // deal with first two sites
            
            
            
            
            // deal with middle sites
            
            // deal with last two sites
                        
            
        }
        
        // deal with last partition, where last bases are invariant

    
    }
    
    
    public static void main(String[] args){
        int[] indices = new int[]{ 0,1,2,3,4,5,6,7,8,9,10,11,12 };
        int[] quint = {0,1,2,3,4};
        
        int[] codon = getCodonSiteIndices(quint, 0, 1);
        ArrayPrinter.print(codon, "\t");
        
        int a = 1;
        int b = 2;
        
        System.out.println("a "+a);
        System.out.println("b "+b);
        
        int c = a;
        
        System.out.println("c "+c);
        c++;
        System.out.println("c "+c);
        System.out.println("a "+a);
    }
    
}
