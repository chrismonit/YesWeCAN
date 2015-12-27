/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.phylo;

import yeswecan.Constants;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class GeneticStructure {
    
    private Partition[] partitions;
   
    // version of constructore which takes strings, as if from command args
    public GeneticStructure(String aLayout, String bLayout, String cLayout,
        String partitionLengths, String delimiter){
        
        // i) gene positions
        String[][] genePositionsString =
        {
            aLayout.split(delimiter),
            bLayout.split(delimiter),
            cLayout.split(delimiter),
        };
        
        int[][] genePositions = new int[genePositionsString.length][genePositionsString[0].length];
        for (int i = 0; i < genePositions.length; i++) {
            for (int j = 0; j < genePositions[0].length; j++) {
                genePositions[i][j] = Integer.parseInt(genePositionsString[i][j]);
            }
        }
        
        // ii) partition lengths
        String[] partitionLengthsString = partitionLengths.split(delimiter);
        int[] lengths = new int[partitionLengthsString.length];
        for (int i = 0; i < lengths.length; i++) {
            lengths[i] = Integer.parseInt(partitionLengthsString[i]);
        }
        
        init(genePositions, lengths);
    }
        
    public GeneticStructure(int[][] genePositions, int[] partitionLengths){
        init(genePositions, partitionLengths);
    }
    
    private void init(int[][] genePositions, int[] partitionLengths){
       int numPartitions = genePositions[0].length;
        this.partitions = new Partition[numPartitions];
        
        int partitionFirstSite = 0;
        for (int iPartition = 0; iPartition < numPartitions; iPartition++) {
            int[] genes = new int[3]; // for each of three frames
            
            for (int jFrame = 0; jFrame < 3; jFrame++) {
                genes[jFrame] = genePositions[jFrame][iPartition];
            }
            
            int partitionLastSite = (partitionFirstSite+partitionLengths[iPartition])-1;
            this.partitions[iPartition] = new Partition(genes, partitionFirstSite, partitionLastSite);
            partitionFirstSite = (partitionFirstSite+partitionLengths[iPartition]); // for the next iteration
        }
    }
    
    @Override
    public String toString(){
        String d = "\t"; // column delimiter
        StringBuilder builder = new StringBuilder(Constants.LAYOUT + d + "Sites");
        for (int i = 0; i < this.partitions.length; i++) {
            builder.append(d);
            String first = Integer.toString(this.partitions[i].firstSite+1);// +1 to correct for zero-based
            String last = Integer.toString(this.partitions[i].lastSite+1);
            builder.append(first+"-"+last);
        }
        
        for (int i = 0; i < 3; i++) {
            builder.append(System.lineSeparator());
            builder.append(Constants.LAYOUT + d);
            builder.append(Constants.FRAMES[i] + d);
            for (int j = 0; j < this.partitions.length; j++) {
                
                builder.append(partitions[j].genes[i] + d);
            }
            //builder.append(System.lineSeparator());
        }
        
        return builder.toString();
    }
    
    
    private class Partition{
        
        private int[] genes; // will be of length 3, one for each frame. Integer value therin represents the gene
        private int firstSite, lastSite;
                
        public Partition(int[] genes, int firstSite, int lastSite){
            this.genes = genes;
            this.firstSite = firstSite;
            this.lastSite = lastSite;
        }
        
        public int getGene(int frame){
            return this.genes[frame];
        }
        
        public int getLength(){
            return (lastSite - firstSite)+1;
        }
        
        public int getFirstSite(){
            return firstSite;
        }
        
        public int getLastSite(){
            return lastSite;
        }
        
    } // class Partition
    
    public static void main(String[] args){
        int[][] genePositions =
        {   {0,1,1,1,2},
            {2,2,0,0,0},
            {3,3,3,3,0}
        };
        int[] lengths = { 10, 20, 30, 40, 50 };
        
        GeneticStructure structure = new GeneticStructure(genePositions, lengths);
        System.out.println(structure.toString());
    }
    
    
    
}// class GeneticStructure
