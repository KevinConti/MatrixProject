public class LinearSystems {
    public static void main(String[] args){
        Matrix linearMatrix = initializeVectorsFromFile("data/LinearSystems.txt");
        Matrix solutionMatrix = linearMatrix.copy();
        try {
            solutionMatrix.inverse();
            System.out.println("Solution to the linear equations:");
            System.out.println(solutionMatrix);
            System.out.println("");
        } catch (InversionException e) {
            e.printStackTrace();
        }
        double determinant = 0.0;
        try {
            Matrix determinantMatrix = linearMatrix.copy();
            determinantMatrix.removeLastColumn();
            determinant = determinantMatrix.determinate();
            System.out.println("Determinant is: ");
            System.out.println(determinant);
            System.out.println("");
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
        System.out.println("");

        double inverseDeterminant = 0.0;
        try {
            inverseDeterminant = inverseMatrix.determinate();
        } catch (InversionException e) {
            e.printStackTrace();
        }
        System.out.format("Inverse Determinate: %n%.5f", inverseDeterminant);
        System.out.println("");

        double productOfDeterminants = inverseDeterminant * determinant;
        System.out.format("Product of determinants: %n%.5f%n", productOfDeterminants);
        System.out.println("");

        try {
            double determinateOfInverse = inverseMatrix.determinate();
            double inverseOfDeterminate = 1/determinant;
            System.out.format("The inverse of the determinant is %.2f%nThe determinate of the inverse is %.2f%n", inverseOfDeterminate, determinateOfInverse);
        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            Matrix conditionNumberMatrix = linearMatrix.copy();
            conditionNumberMatrix.removeLastColumn();
            double conditionNumber = conditionNumberMatrix.conditionNumber();
            System.out.format("%nThe condition number for this matrix is %f%n", conditionNumber);
        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            Matrix squareMatrix = linearMatrix.copy();
            squareMatrix.removeLastColumn();
            System.out.println(Matrix.multiply(squareMatrix, inverseMatrix));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static Matrix initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertLinear(filepath);
    }
}
