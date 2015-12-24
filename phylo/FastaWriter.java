/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.phylo;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import pal.alignment.Alignment;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class FastaWriter {
    
    String newline = "\n";
    
    public String fastaString(Alignment alignment){
        StringBuilder builder = new StringBuilder();
        for (int iSequence = 0; iSequence < alignment.getSequenceCount(); iSequence++) {
            builder.append(">"+alignment.getIdentifier(iSequence).getName());
            builder.append(newline);
            builder.append(alignment.getAlignedSequenceString(iSequence));
            builder.append(newline);
        }
        return builder.toString();
    }
    
    public void writeFasta(Alignment alignment, String outputPath){
        try (Writer writer = new BufferedWriter(new OutputStreamWriter(
              new FileOutputStream(outputPath), "utf-8"))) {
       writer.write( fastaString(alignment) );
        }
        catch(Exception e){
        }
        
    } // writeFasta
    
}
