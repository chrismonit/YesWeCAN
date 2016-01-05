/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model;

import java.io.FileReader;
import java.util.Arrays;
import java.util.List;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.Nucleotides;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Parameter;
import yeswecan.Constants;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.ReorderFrequencies;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class HKYlnl {
    
    AdvancedAlignment alignment;
    Tree tree;
    
    public static void main(String[] args){
        
        new HKYlnl();
        
        
    }
    
    public HKYlnl(){
    
        loadData("/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/calculate/0.5kappa/hky/yang/sim/fastas/1.fasta",
                "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/calculate/0.5kappa/hky/yang/fit/tree.tre"
        );
        
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(0.5);
        double[] tcag = new double[]{0.10170,0.20036,0.30310,0.39484}; // TCAG
        BaseFrequencies freq = new BaseFrequencies(ReorderFrequencies.pamlToAlpha(tcag));
        
        List<Parameter> parameters = Arrays.asList(kappa, freq);
        
        double[] optSpaceParams = Mapper.getOptimisable(parameters);
        
        FunctionHKY calculator = new FunctionHKY(this.alignment, this.tree);
        System.out.println("lnl: " + calculator.value(optSpaceParams));
    }
    
    
    public void loadData(String alignmentPath, String treePath){
        try{
            this.alignment = new AdvancedAlignment(
                                new SimpleAlignment(
                                        AlignmentReaders.readFastaSequences(new FileReader(alignmentPath), new Nucleotides())));
            this.tree = new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println(Constants.ERROR_PREFIX + "Unable to load alignment or tree file(s)");
            e.printStackTrace();
            System.exit(1);
        }

    }//loadData
    
}
