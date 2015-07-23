/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.optim;
import java.util.Hashtable;
/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class ParamMap {
    
    /* have method which takes a double[] and can populate a useful hashtable
    But this class could take on a more complicated role of both assigning 
    'point' array elements to params AND housing fixed param values
    
    algorithm:
point
endlist = {}


count = 0 // counts how far along the point array we have moved
for param in range(N)
	if param.fixed == true
		for i in range(param.value.length)
			endlist.append(param.value[i])
	else
		for i in range(count, param.value.length)
			endlist.append(point[i])
		count += param.value.length



convert end list into something accessible


    */
    private Hashtable<String, Double> table;
    private String[] identifiers;
    
    public ParamMap(String[] identifiers){
        this(new double[identifiers.length], identifiers);
    }
    
    
    public ParamMap(double[] params, String[] identifiers){
        this.table = new Hashtable<String, Double>();
        this.identifiers = identifiers;
        
        if (params.length != identifiers.length){
            throw new RuntimeException("ParamMap: Lengths of params and identifiers do not match");
        }
        
        for (int i = 0; i < params.length; i++) {
            table.put(identifiers[i], params[i]);
        }
    
    }
    
    // takes point array as given to Function.value() and places values in datastructure
    public void update(double[] params){
        if (params.length != identifiers.length){
            throw new RuntimeException("ParamMap.update: incorrect length of params");
        }
        
        for (int i = 0; i < identifiers.length; i++) {
            table.put(identifiers[i], params[i]);
        }

    }

}//class
