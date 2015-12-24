
package yeswecan.utils;

/**
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * written originally for omegaABC/omega3 model, which is why the original methods only allowed 3 rows
 * 
 * Class with all static members. 

 * 
 */
public final class  MatrixPrinter {
    
    private MatrixPrinter(){} //no instantiation
    
    public static double roundDouble(double toRound, int decPlaces){
        int tenMultiple = (int)Math.pow(10.0, (double)decPlaces);
        return (double)Math.round( toRound * tenMultiple )  / tenMultiple;
    }
    
    public static void Print3RowMatrix(int[][] matrix2D, String message){
        
        StringBuilder builderRow0 = new StringBuilder("");
        StringBuilder builderRow1 = new StringBuilder("");
        StringBuilder builderRow2 = new StringBuilder("");

        for (int j = 0; j < matrix2D[0].length; j++) {   
            builderRow0.append(matrix2D[0][j] + "\t");
            builderRow1.append(matrix2D[1][j] + "\t");
            builderRow2.append(matrix2D[2][j] + "\t");
        }

        System.out.println(message);
        System.out.println(builderRow0.toString());
        System.out.println(builderRow1.toString());
        System.out.println(builderRow2.toString());
        System.out.println("------------");
 
    }
    
    
    public static void Print3RowMatrix(double[][] matrix2D, String message){
        
        StringBuilder builderRow0 = new StringBuilder("");
        StringBuilder builderRow1 = new StringBuilder("");
        StringBuilder builderRow2 = new StringBuilder("");

        for (int j = 0; j < matrix2D[0].length; j++) {   
            builderRow0.append(matrix2D[0][j] + "\t");
            builderRow1.append(matrix2D[1][j] + "\t");
            builderRow2.append(matrix2D[2][j] + "\t");
        }

        System.out.println(message);
        System.out.println(builderRow0.toString());
        System.out.println(builderRow1.toString());
        System.out.println(builderRow2.toString());
        System.out.println("------------");
 
    }
    
    public static void PrintMatrix(double[][] matrix2D, String message){
        StringBuilder[] builders = new StringBuilder[ matrix2D.length ];
        
        for (int iRow = 0; iRow < matrix2D.length; iRow++) {
            builders[iRow] = new StringBuilder("");
        }
        
        for (int iRow = 0; iRow < matrix2D.length; iRow++) {
            for (int jColumn = 0; jColumn < matrix2D[0].length; jColumn++) {
                builders[iRow].append( matrix2D[iRow][jColumn] + "\t" );
            }
        }
        
        System.out.println(message);
        for (int iRow = 0; iRow < matrix2D.length; iRow++) {
            System.out.println( builders[iRow].toString() );
        }
        System.out.println("------------");
        
    }
    
    public static void PrintMatrix(double[][] matrix2D, String message, String round){
        int sigFigures = 4;
        
        StringBuilder[] builders = new StringBuilder[ matrix2D.length ];
        
        for (int iRow = 0; iRow < matrix2D.length; iRow++) {
            builders[iRow] = new StringBuilder("");
        }
        
        for (int iRow = 0; iRow < matrix2D.length; iRow++) {
            for (int jColumn = 0; jColumn < matrix2D[0].length; jColumn++) {
                builders[iRow].append( roundDouble(matrix2D[iRow][jColumn], sigFigures) + "\t" );
            }
        }
        
        System.out.println(message);
        for (int iRow = 0; iRow < matrix2D.length; iRow++) {
            System.out.println( builders[iRow].toString() );
        }
        System.out.println("------------");
        
    }
    
    
    
}//class
