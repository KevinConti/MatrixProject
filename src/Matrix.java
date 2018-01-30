import java.util.Arrays;

public class Matrix {

    double[][] matrix;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(double[][] matrix) {
        this.matrix = matrix;
    }

    public void setMatrix(int row, int column, double value){
        matrix[row][column] = value;
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
        Vector meanVector = Vector.mean(vectors);
        //Subtract mean
        Vector[] subtractedVectors = Vector.subtract(vectors, meanVector);
        //convert to matrices
        Matrix[] matrices = Matrix.parseVectors(subtractedVectors);
        //multiply by inverses
        Matrix[] twoByTwoMatrices = multiplyByInverse(matrices);
        System.out.println(Arrays.toString(twoByTwoMatrices));
        //Calculate matrix mean
        //return answer


        return matrices[0];
    }

    private static Matrix[] multiplyByInverse(Matrix[] matrices){
        Matrix[] twoByTwoMatrices = new Matrix[matrices.length];

        //for each matrix in matrices:
        for (int i = 0; i < matrices.length; i++) {
            Matrix matrix = matrices[i];
            //Determine the inverse
            Matrix inverseMatrix = Matrix.inverse(matrix);
            //Multiply the inverse times matrix
            try {
                //Store two by two result in twoByTwoMatrices
                Matrix twoByTwoMatrix = Matrix.multiply(inverseMatrix, matrix);
                twoByTwoMatrices[i] = twoByTwoMatrix;
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        //return result
        return twoByTwoMatrices;
    }

    //NOTE: This method currently only works to invert a 1x2 matrix into a 2x1
    public static Matrix inverse(Matrix matrix){
        double[][] invertedTable = new double[2][1];
        invertedTable[0][0] = matrix.getMatrix()[0][0];
        invertedTable[1][0] = matrix.getMatrix()[0][1];

        return new Matrix(invertedTable);
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
