import java.util.ArrayList;

public class main {
    public static void main(String[] args){
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

        //Problem 7
        try {
            findMisclassifiedPoints(vectors, discriminates);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static void findMisclassifiedPoints(Vector[][] vectors, Discriminate[] discriminates) throws Exception {
        ArrayList<Vector> misclassifiedVectors = new ArrayList<Vector>();
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
                        misclassifiedVectors.add(vector);
                        printTableLine(vector, g1Value, g2Value);
                    }
                }
                if(i == 1){ //Class two vectors
                    Vector vector = vectors[i][j];
                    double g1Value = discriminates[0].classify(vector.getX(), vector.getY());
                    double g2Value = discriminates[1].classify(vector.getX(), vector.getY());
                    if(g2Value >= g1Value) {} //Vector is correctly classified
                    else{
                        misclassifiedVectors.add(vector);
                        printTableLine(vector, g1Value, g2Value);
                    }
                }
            }
        }
        //Create bottom line of table
        System.out.format("--------------------------------------------------------------%n");

        //Create summary table
        createSummaryTable(misclassifiedVectors);
    }

    private static void printTableLine(Vector vector, double g1Value, double g2Value){
        String format = "<%+f.5, %+f.5>         %+-17.5f%+f%n";
        System.out.format(format, vector.getX(), vector.getY(), g1Value, g2Value);
    }

    private static void createSummaryTable(ArrayList<Vector> misclassifiedVectors){

    }

    private static Vector[][] initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }
}
