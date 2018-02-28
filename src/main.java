import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class main {
    public static void main(String[] args){
//        runMatrixProject();
        try {
            runEigenProject();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void runEigenProject() throws Exception {
        //1) Find mean vector and covariance matrix
        Vector[] vectors = initializeSingleSetOfVectorsFromFile("data/p2data.txt");
        Matrix covarianceMatrix = Matrix.toCovariance(vectors);
        System.out.printf("The covariance matrix is: \n%s\n", covarianceMatrix);

        //2) Determine the trace of the covariance matrix
        double trace = covarianceMatrix.trace();
        System.out.printf("The trace of the covariance matrix is %f\n", trace);


        //The determinant of the covariance matrix
        double determinant = covarianceMatrix.determinate();
        System.out.printf("\nThe determinant is %f\n", determinant);
    }

    private static void runMatrixProject(){
        Vector[][] vectors = initializeVectorsFromFile("data/data.txt");
        System.out.println("Class One:");
        Matrix classOneCovariance = Matrix.toCovariance(vectors[0]);
        System.out.println("");
        System.out.println("Class One Covariance:");
        System.out.println(classOneCovariance);
        System.out.println("Class Two:");
        Matrix classTwoCovariance = Matrix.toCovariance(vectors[1]);
        System.out.println("");
        System.out.println("Class Two Covariance:");
        System.out.println(classTwoCovariance);

        //Next, test inverses
        Matrix identityMatrix = Matrix.createIdentityMatrix(classOneCovariance.numRows());
        System.out.println("");
        System.out.println("Class one inverse matrix:");
        Matrix classOneInverse = classOneCovariance.copy();
        try {
            classOneInverse.inverse(identityMatrix);
            //Remove first two columns to remove identity matrix
            for(int i = 0; i < identityMatrix.numColumns(); i++) {
                classOneInverse.removeFirstColumn();
            }
            System.out.println(classOneInverse);
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("");
        System.out.println("Class two inverse matrix:");
        Matrix classTwoInverse = classTwoCovariance.copy();
        try {
            classTwoInverse.inverse(identityMatrix);
            //Remove first two columns to remove identity matrix
            for(int i = 0; i < identityMatrix.numColumns(); i++) {
                classTwoInverse.removeFirstColumn();
            }
            System.out.println(classTwoInverse);
        } catch (Exception e) {
            e.printStackTrace();
        }

        //Next, find determinants
        System.out.println("Class one determinant");
        double classOneDeterminate = 0.0;
        try {
            classOneDeterminate = classOneCovariance.determinate();
            System.out.println(Double.toString(classOneDeterminate));
        } catch (InversionException e) {
            e.printStackTrace();
        }
        System.out.println("Class Two determinant");
        double classTwoDeterminate = 0.0;
        try {
            classTwoDeterminate = classTwoCovariance.determinate();
            System.out.println(Double.toString(classTwoDeterminate));
        } catch (InversionException e) {
            e.printStackTrace();
        }

        //Create Discriminant objects
        Vector classOneMeanVector = Vector.mean(vectors[0]);
        Discriminate classOneDiscriminant = new Discriminate(classOneMeanVector, classOneInverse, classOneDeterminate);
        Vector classTwoMeanVector = Vector.mean(vectors[1]);
        Discriminate classTwoDiscriminant = new Discriminate(classTwoMeanVector, classTwoInverse, classTwoDeterminate);
        Discriminate[] discriminates = {classOneDiscriminant, classTwoDiscriminant};

        //Problem 6
        try {
            double g1m1Value = classOneDiscriminant.classify(classOneMeanVector.getX(), classOneMeanVector.getY());
            double g2m1Value = classTwoDiscriminant.classify(classOneMeanVector.getX(), classOneMeanVector.getY());
            double g1m2Value = classOneDiscriminant.classify(classTwoMeanVector.getX(), classTwoMeanVector.getY());
            double g2m2Value = classTwoDiscriminant.classify(classTwoMeanVector.getX(), classTwoMeanVector.getY());
            System.out.format("%ng1 of mu1 =%f%ng2 of mu1 is %f%n%ng1 of mu2 =%f%ng2 of mu2 is %f%n",g1m1Value,g2m1Value,g1m2Value,g2m2Value);
            System.out.format("Therefore, mu1 should be in g1, and mu2 should be in g2%n");
        } catch (Exception e) {
            e.printStackTrace();
        }

        //Problem 7
        try {
            findMisclassifiedPoints(vectors, discriminates);
        } catch (Exception e) {
            e.printStackTrace();
        }

        //Estimate Boundary Contour (Problem 8)
        try {
            estimateBoundaryContour("data/out/boundaries.csv", discriminates);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void findMisclassifiedPoints(Vector[][] vectors, Discriminate[] discriminates) throws Exception {
        ArrayList<Vector> classOneMisclassifiedVectors = new ArrayList<>(); //Vectors that should be in class one
        ArrayList<Vector> classTwoMisclassifiedVectors = new ArrayList<>(); //Vectors that should be in class two
        //Create Top line of table
        System.out.format("--------------------------------------------------------------%n");
        System.out.format("%20s%20s%20s%n", "Misclassified", "g1(x)", "g2(x)");
        System.out.format("--------------------------------------------------------------%n");
        for(int i = 0; i < vectors.length; i++){
            for(int j = 0; j < vectors[i].length; j++){
                if(i == 0){ //class one vectors
                    Vector vector = vectors[i][j];
                    double g1Value = discriminates[0].classify(vector.getX(), vector.getY());
                    double g2Value = discriminates[1].classify(vector.getX(), vector.getY());
                    if(g1Value >= g2Value) {} //Vector is correctly classified
                    else{
                        classOneMisclassifiedVectors.add(vector);
                        printTableLine(vector, g1Value, g2Value);
                    }
                }
                if(i == 1){ //Class two vectors
                    Vector vector = vectors[i][j];
                    double g1Value = discriminates[0].classify(vector.getX(), vector.getY());
                    double g2Value = discriminates[1].classify(vector.getX(), vector.getY());
                    if(g2Value >= g1Value) {} //Vector is correctly classified
                    else{
                        classTwoMisclassifiedVectors.add(vector);
                        printTableLine(vector, g1Value, g2Value);
                    }
                }
            }
        }
        //Create bottom line of table
        System.out.format("--------------------------------------------------------------%n");

        //Create summary table
        createSummaryTable(classOneMisclassifiedVectors, classTwoMisclassifiedVectors);
    }

    private static void printTableLine(Vector vector, double g1Value, double g2Value){
        String format = "<%+f.5, %+f.5>         %+-17.5f%+f%n";
        System.out.format(format, vector.getX(), vector.getY(), g1Value, g2Value);
    }

    private static void createSummaryTable(ArrayList<Vector> classOneErrors, ArrayList<Vector> classTwoErrors){
        //Create beginning of table
        System.out.format("--------------------------------------------------------------%n");
        System.out.format("%20s%20s%20s%n", "", "Class 1", "Class 2");
        System.out.format("--------------------------------------------------------------%n");
        System.out.format("%-20s%15d%20d%n", "Correctly Identified: ", 110-classOneErrors.size(), 110-classTwoErrors.size());
        System.out.format("%-17s%13d%20d%n", "Incorrectly Identified: ", classTwoErrors.size(), classOneErrors.size());
        System.out.format("--------------------------------------------------------------%n");
    }

    //This method computes the boundaries points in cartesian axes, and outputs them to a .csv file for external use
    private static void estimateBoundaryContour(String outputFilePath, Discriminate[] discriminates) throws Exception {
        final double SIGMA = 0.01;

        PrintWriter out = null;
        try {
            out = new PrintWriter(new FileWriter(outputFilePath));
        } catch (IOException e) {
            e.printStackTrace();
        }
        out.printf("boundaryX,boundaryY%n");
        for(double x = -5.0; x < 4; x+= 0.01){
            for(double y = -4.0; y<8.0; y+= 0.01){
                double g1 = discriminates[0].classify(x,y);
                double g2 = discriminates[1].classify(x,y);
                if(Math.abs(g1 - g2) < SIGMA) {
                    out.printf("%f,%f%n", x, y);
                }
            }
        }
        out.close();
    }

    private static Vector[][] initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }

    private static Vector[] initializeSingleSetOfVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertP2data(filepath);
    }
}
