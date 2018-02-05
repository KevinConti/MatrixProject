import java.util.Arrays;

public class Matrix {

    double[][] matrix;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public double getValueAt(int row, int column){
        return matrix[row][column];
    }

    public void setMatrix(double[][] matrix) {
        this.matrix = matrix;
    }

    public void setMatrix(int row, int column, double value){
        matrix[row][column] = value;
    }

    public int numColumns(){
        return this.matrix[0].length;
    }

    public int numRows(){
        return this.matrix.length;
    }

    //Class Methods
    public static Matrix add(Matrix matrixOne, Matrix matrixTwo){
        //Check if matrices are appropriately sized (both must be MxN, returns a MxN matrix
        //Add the two together, store the result in the a new matrix that is returned
        int numberOfRows = matrixOne.getMatrix().length;
        int numberOfColumns = matrixOne.getMatrix()[0].length;
        double[][] addedMatrixArray = new double[numberOfRows][numberOfColumns];
        for(int row = 0; row < numberOfRows; row++){
            for(int column = 0; column < numberOfColumns; column++){
                addedMatrixArray[row][column] = matrixOne.getMatrix()[row][column] + matrixTwo.getMatrix()[row][column];
            }
        }
        return new Matrix(addedMatrixArray);
    }

    //Warning, this is a destructive method.
    public static void divideByScalarDestructive(Matrix matrix, double scalar){
        //enumerate through each entry in matrix table, divide by the scalar and set that new value into the matrix
        int numRows = matrix.getMatrix().length;
        int numColumns = matrix.getMatrix()[0].length;
        for (int i = 0; i < numRows; i++){
            for (int j = 0; j < numColumns; j++){
                double matrixValue = matrix.getMatrix()[i][j];
                matrix.setMatrix(i, j, matrixValue/scalar);
            }
        }
    }

    public static Matrix createAugmentedMatrix(Matrix squareMatrix, Matrix coefficientMatrix) throws Exception {
        //Determine what length the matrix will be / verify appropriate length
        if (squareMatrix.numColumns() == squareMatrix.numRows() && squareMatrix.numRows() == coefficientMatrix.numRows()){
            //create a new matrix with the appropriate size
            Matrix augmentedMatrix = new Matrix(new double[squareMatrix.numRows()][squareMatrix.numColumns() + 1]);

            //Set matrix values
            //For each row:
            for (int i = 0; i < augmentedMatrix.numRows(); i++){

                for (int j = 0; j < augmentedMatrix.numColumns(); j++){
                    //First n columns (where n = squareMatrix.numColumns) will be square matrix's values
                    if (j < augmentedMatrix.numColumns() - 1){
                        double value = squareMatrix.getValueAt(i,j);
                        augmentedMatrix.setMatrix(i, j, value);
                    }
                    else {
                        //Last column will be coefficientMatrix's value
                        augmentedMatrix.setMatrix(i,j, coefficientMatrix.getValueAt(i, 0));
                    }

                }

            }
            return augmentedMatrix;
        }
        else throw new Exception("Input matrices not formatted properly");
    }

    public static Matrix parseVector(Vector vector){

        double x = vector.getX();
        double y = vector.getY();

        double[][] table = new double[1][2];
        table[0][0] = x;
        table[0][1] = y;

        return new Matrix(table);
    }

    public static Matrix toCovariance(Vector[] vectors){
        //calculate mean of vectors
        System.out.println("");
        Vector meanVector = Vector.mean(vectors);
        System.out.println("Mean Vector:");
        System.out.println(meanVector);
        //Subtract mean
        Vector[] subtractedVectors = Vector.subtract(vectors, meanVector);
        //convert to matrices
        Matrix[] matrices = Matrix.parseVectors(subtractedVectors);
        System.out.println("");
        System.out.println("Matrices after subtracting mean");
        System.out.println(Arrays.toString(matrices));
        //multiply by inverses
        System.out.println("");
        Matrix[] twoByTwoMatrices = multiplyByTranspose(matrices);
        System.out.println("After multiplying by inverse");
        System.out.println(Arrays.toString(twoByTwoMatrices));
        //Calculate matrix mean
        System.out.println("");
        System.out.println("Covariance Matrix:");
        Matrix matrixMean = matrixMean(twoByTwoMatrices);
        System.out.println(matrixMean);
        //return answer]
        return matrixMean;
    }

    private static Matrix[] multiplyByTranspose(Matrix[] matrices){
        Matrix[] twoByTwoMatrices = new Matrix[matrices.length];

        //for each matrix in matrices:
        for (int i = 0; i < matrices.length; i++) {
            Matrix matrix = matrices[i];
            //Determine the inverse
            Matrix transposedMatrix = Matrix.transpose(matrix);
            //Multiply the inverse times matrix
            try {
                //Store two by two result in twoByTwoMatrices
                Matrix twoByTwoMatrix = Matrix.multiply(transposedMatrix, matrix);
                twoByTwoMatrices[i] = twoByTwoMatrix;
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        //return result
        return twoByTwoMatrices;
    }

    //NOTE: This method currently only works to invert a 1x2 matrix into a 2x1
    public static Matrix transpose(Matrix matrix){
        double[][] transposedTable = new double[2][1];
        transposedTable[0][0] = matrix.getMatrix()[0][0];
        transposedTable[1][0] = matrix.getMatrix()[0][1];

        return new Matrix(transposedTable);
    }

    public static Matrix multiply(Matrix matrixOne, Matrix matrixTwo) throws Exception {
        Matrix multipliedMatrix;

        //Assuming a matrixOne is Mx*y, and matrixTwo is Mi*j, then they can be multiplied if y==i
        int matrixOneNumberOfColumns = matrixOne.getMatrix()[0].length;
        int matrixTwoNumberOfRows = matrixTwo.getMatrix().length;
        if(matrixOneNumberOfColumns == matrixTwoNumberOfRows) {
            int numberOfMultiplicationsPerCell = matrixOneNumberOfColumns;

            //Assuming a matrixOne is Mx*y, and matrixTwo is Mi*j, then the resulting matrix will be Mx*j
            int matrixOneNumberOfRows = matrixOne.getMatrix().length;
            int matrixTwoNumberOfColumns = matrixTwo.getMatrix()[0].length;

            int numberOfColumnsForAnswer = matrixOneNumberOfRows;
            int numberOfRowsForAnswer = matrixTwoNumberOfColumns;
            multipliedMatrix = new Matrix(new double[numberOfRowsForAnswer][numberOfColumnsForAnswer]);

            //Multiplication begins. Formula is answerMatrix[x][y] = sum(matrixOne[x][i] + matrixTwo[i][y])
            //Where i = numberOfMultiplicationsPerCell
            int x,y;
            for(x = 0; x < numberOfRowsForAnswer; x++) {
                for(y = 0; y < numberOfColumnsForAnswer; y++) {
                    //for each cell...
                    double sum = 0.0;
                    for (int i = 0; i < numberOfMultiplicationsPerCell; i++) {
                        sum += matrixOne.getMatrix()[x][i] * matrixTwo.getMatrix()[i][y];
                    }

                    multipliedMatrix.setMatrix(x, y, sum);
                }
            }

        }
        else {
            throw new Exception("Matrices cannot be multiplied");
        }
        return multipliedMatrix;
    }

    //WARNING: This method requires that matrix sizes are uniform in the given array
    public static Matrix matrixMean(Matrix[] matrixArray){
        int numOfRowsInMatrices = matrixArray[0].getMatrix().length;
        int numOfColumnsInMatrices = matrixArray[0].getMatrix()[0].length;
        Matrix matrixMean = new Matrix(new double[numOfRowsInMatrices][numOfColumnsInMatrices]);
        try {
            int i;
            for (i = 0; i < matrixArray.length; i++){
                matrixMean = Matrix.add(matrixMean, matrixArray[i]);
            }
            Matrix.divideByScalarDestructive(matrixMean, i);
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println("Some error in matrixMean method");
        }

        return matrixMean;
    }

    public static Matrix[] parseVectors(Vector[] vectors){
        Matrix[] convertedMatrices = new Matrix[vectors.length];
        for (int i = 0; i < vectors.length; i++){
            convertedMatrices[i] = Matrix.parseVector(vectors[i]);
        }

        return convertedMatrices;
    }

    @Override
    public String toString() {
        String result = "";
        for(int i = 0; i < matrix.length; i++){
            String rowString = "[";
            for(int j = 0; j < matrix[0].length; j++){
                rowString += Double.toString(matrix[i][j]);
                if(j != matrix[0].length - 1) {
                    rowString += " ";
                }
            }
            rowString += "]\n";
            result += rowString;
        }
        return result;
    }
}
