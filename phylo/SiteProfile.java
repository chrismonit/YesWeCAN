/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.phylo;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SiteProfile {
    
    
    private int siteIndex;
    private int aCodonNumber; // not sure I need these
    private int bCodonNumber;
    private int cCodonNumber;
    
    private int siteType;
    private int[] aCodon;
    private int[] bCodon;
    private int[] cCodon;
    
    public SiteProfile(int siteIndex, int[] aCodon, int[] bCodon, int[] cCodon){
        this.siteIndex = siteIndex;
        this.aCodon = aCodon;
        this.bCodon = bCodon;
        this.cCodon = cCodon;
  
    }
    
    public int getSiteIndex(){ 
        return this.siteIndex; 
    }
    
    public int[] getACodon(){
        return this.aCodon;
    }
    
    public int[] getBCodon(){
        return this.aCodon;
    }
    
    public int[] getCCodon(){
        return this.aCodon;
    }
    
}// class
