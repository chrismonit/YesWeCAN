package yeswecan.phylo;

import pal.alignment.Alignment;
import pal.alignment.SimpleAlignment;
import pal.datatype.DataTypeTool;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * I find the PAL alignment objects to be less suited to my design plan.
 * This is therefore an extension of the PAL SimpleAlignment class,
 * which includes the kinds of methods I would like to use.
 * 
 * Should consider creating an alignment wrapper which is suited to dealing with overlapping frames
 * (e.g. conveniently returns only alpha sites, or sites in a given partition)
 */
public class AdvancedAlignment extends SimpleAlignment {

    public AdvancedAlignment(Alignment a){
        super(a);
    }
    
    public int getStateBySequenceName(String name, int site){ 
        //might be faster if I didn't assign all of these variables?
        int sequenceID = super.whichIdNumber(name);
        char stateAsChar = super.getData(sequenceID, site);
        int stateAsInt = DataTypeTool.getNucleotides().getState(stateAsChar); 
        
        return stateAsInt;
    }
    
    
    
}// class

