public class LinearSystems {
    public static void main(String[] args){
        Matrix linearMatrix = initializeVectorsFromFile("data/LinearSystems.txt");
        Matrix solutionMatrix = linearMatrix.copy();
        try {
            solutionMatrix.inverse();
            System.out.println("Solution to the linear equations:");
            System.out.println(solutionMatrix);
        } catch (InversionException e) {
            e.printStackTrace();
        }
        double determinant = 0.0;
        try {
            Matrix determinantMatrix = linearMatrix.copy();
            determinantMatrix.removeLastColumn();
            determinant = determinantMatrix.determinate();
            System.out.println("Determinant is: "+Double.toString(determinant));
        } catch (InversionException e) {
            e.printStackTrace();
        }

        Matrix inverseMatrix = linearMatrix.copy();
        inverseMatrix.removeLastColumn();
        Matrix identityMatrix = Matrix.createIdentityMatrix(inverseMatrix.numRows());
        try {
            inverseMatrix.inverse(identityMatrix);
            for(int i = 0; i < identityMatrix.numColumns(); i++){
                inverseMatrix.removeFirstColumn();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Inverse matrix: ");
        System.out.println(inverseMatrix);

        double inverseDeterminant = 0.0;
        try {
            inverseDeterminant = inverseMatrix.determinate();
        } catch (InversionException e) {
            e.printStackTrace();
        }

        double productOfDeterminants = inverseDeterminant * determinant;
        System.out.println("Product of determinants: ");
        System.out.println(productOfDeterminants);

//        //Next, test inverses
//        Matrix identityMatrix = Matrix.createIdentityMatrix(classOneCovariance.numRows());
//        System.out.println("");
//        System.out.println("Class one inverse matrix:");
//        Matrix classOneInverse = classOneCovariance.copy();
//        try {
//            classOneInverse.inverse(identityMatrix);
//            classOneInverse.removeFirstColumn();
//            classOneInverse.removeFirstColumn();
//            System.out.println(classOneInverse);
//        } catch (Exception e) {
//            e.printStackTrace();
//            fail();
//        }
    }

    private static Matrix initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertLinear(filepath);
    }
}
