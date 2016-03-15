/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.phylo;

import java.util.List;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SequenceProfile {
    
    private SiteProfile[] siteProfiles;
    private GeneticStructure genStruct;
    public static final int nTERMINAL_IGNORE = 2;
    
    public SequenceProfile(GeneticStructure genStruct){
        this.genStruct = genStruct;
        this.siteProfiles = new SiteProfile[(genStruct.getTotalLength() - nTERMINAL_IGNORE*2 )]; // ignoring first and last two positions in sequence
        init();
    }
    
    private void init(){
        
        // deal with first partition, where first bases are invariant
        
        for (int iPartition = 1; iPartition < (genStruct.getNumberOfPartitions() - 1); iPartition++) {
            
            // deal with first two sites
            
            
            
            
            // deal with middle sites
            
            // deal with last two sites
                        
            
        }
        
        // deal with last partition, where last bases are invariant

    
    }
    
    
}
