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
import yeswecan.Constants;

/**
 *
 * @author cmonit1
 */
public abstract class SiteClassTable {
    
    protected List<List<String>> table;
    protected String majorDelimiter = Constants.MAJOR_DELIM;
    protected String minorDelimiter = ","; //Constants.DEL;
    
    public SiteClassTable(){
        this.table = new ArrayList<List<String>>();
    }
    
    // pointless?
    public SiteClassTable(List<List<String>> table){
        this.table = table;
    }
    
    
    public void addDelimiterColumn(String delimiter){
        for (List<String> stringList : this.table){
            stringList.add(delimiter);
        }        
    }
    // could have equivalent for adding rows
    
    public void write(String filePath){
        
        try{
            FileWriter fileWriter = new FileWriter(filePath);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            
            for (List<String> stringList : this.table){
                bufferedWriter.write(String.join(this.minorDelimiter, stringList));
            }
            
            bufferedWriter.close();
        
        }catch (IOException e){
            e.printStackTrace();
        }
        
    }
    
    public void print(String lineMarker){
        for (List<String> stringList : this.table){
            //System.out.println(lineMarker+this.minorDelimiter+String.join(this.minorDelimiter, stringList));
            System.out.println(String.join(this.minorDelimiter, stringList));
        }
    }
    
    protected abstract void makeTable();
    
    protected void prependHeader(int numberOfSiteClasses, int numberOfGenes){
        List<String> header = new ArrayList<String>();
        header.add(""); // blank (could name column SITES or something)
        header.add(this.majorDelimiter);
        for (int iGene = 1; iGene < numberOfGenes+1; iGene++) {
            
            for (int iSiteClass = 0; iSiteClass < numberOfSiteClasses; iSiteClass++) {
                header.add(Integer.toString(iGene));
            }
            header.add(this.majorDelimiter);
        }
        
        this.table.add(0, header);
    }
    
    
}
