/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 *
 * @author cmonit1
 */
public abstract class SiteClassTable {
    
    protected List<List<String>> table;
    
    public SiteClassTable(){
        
    }
    
    public void addDelimiterColumn(String delimiter){
        for (List<String> stringList : this.table){
            stringList.add(delimiter);
        }        
    }
    // could have equivalent for adding rows
    
    public void write(String filePath, String delimiter){
        
        try{
            FileWriter fileWriter = new FileWriter(filePath);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            
            for (List<String> stringList : this.table){
                bufferedWriter.write(String.join(delimiter, stringList));
            }
            
            bufferedWriter.close();
        
        }catch (IOException e){
            e.printStackTrace();
        }
        
    }
    
    
    public void print(String lineMarker, String delimiter){
        for (List<String> stringList : this.table){
            System.out.println(lineMarker+delimiter+String.join(delimiter, stringList));
        }
    }
    
}
