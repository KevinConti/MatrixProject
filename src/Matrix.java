import java.util.ArrayList;
import java.util.Arrays;
import java.lang.Math;

public class Matrix {

    double[][] matrix;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
    }

    public Matrix(int numRows, int numColumns){
        this.matrix = new double[numRows][numColumns];
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

    public Matrix copy(){
        Matrix copy = new Matrix(this.numRows(), this.numColumns());
        for (int rowIndex = 0; rowIndex < copy.numRows(); rowIndex++){
            for (int columnIndex = 0; columnIndex < copy.numColumns(); columnIndex++){
                copy.setMatrix(rowIndex,columnIndex, this.getValueAt(rowIndex, columnIndex));
            }
        }
        return copy;
    }

    public void removeFirstColumn(){
        Matrix tempMatrix = new Matrix(this.numRows(), this.numColumns() - 1);
        for (int i = 0; i < this.numRows(); i++){
            for (int x = 0; x < this.numColumns(); x++){
                if (x != 0){
                    tempMatrix.setMatrix(i,x-1, this.getValueAt(i,x));
                }
            }
        }
        this.setMatrix(tempMatrix.getMatrix());
    }

    public void removeLastColumn(){
        Matrix tempMatrix = new Matrix(this.numRows(), this.numColumns() - 1);
        for (int i = 0; i < this.numRows(); i++){
            for (int x = 0; x < this.numColumns(); x++){
                if (x != this.numColumns()-1){
                    tempMatrix.setMatrix(i,x, this.getValueAt(i,x));
                }
            }
        }
        this.setMatrix(tempMatrix.getMatrix());
    }

    public int largestAbsoluteValueIndex(int currentColumn){
        int index = -1;
        double highestValue = 0.0;
        //i does not start at 0 as we do not check if the highest value is in a row we don't care about
        //i is doing two things here, it marks the earliest index we need to check, and it refers to the
        //total number of rows in the matrix
        for (int i = currentColumn; i < this.numRows(); i++){
            double currentValue = Math.abs(this.getValueAt(i, currentColumn));
            if (currentValue > highestValue){
                highestValue = currentValue;
                index = i;
            }
        }
        return index;
    }

    public double[] rowSums(){
        double[] sums = new double[this.numRows()];
        for(int i = 0; i < numRows(); i++){
            double currentRowValue = 0.0;
            for(int j = 0; j < numColumns(); j++){
                double absCellValue = Math.abs(this.getValueAt(i,j));
                currentRowValue += absCellValue;
            }
            sums[i] = currentRowValue;
        }
        return sums;
    }

    public double conditionNumber(){
        Matrix copyMatrix = this.copy();
        Matrix identity = createIdentityMatrix(this.numRows());
        try {
            copyMatrix.inverse(identity);
        } catch (Exception e) {
            e.printStackTrace();
        }
        for(int i = 0; i < identity.numColumns(); i++){ //Purge identity matrix from copyMatrix
            copyMatrix.removeFirstColumn();
        }

        double[] inverseRowSums = copyMatrix.rowSums();
        double[] thisRowSums = this.rowSums();

        double inverseMax = Matrix.maximum(inverseRowSums);
        double thisMax = Matrix.maximum(thisRowSums);

        return thisMax * inverseMax;
    }

    public void inverse() throws InversionException {
        //'this' is an augmented coefficient matrix
        for(int j = 0; j < this.numRows(); j++){
            //2a: Find largest absolute value in column
            int p = this.largestAbsoluteValueIndex(j);
            //2b: Check for flag, exit if answer cannot be found
            if (this.getValueAt(p,j) == 0.0){
                throw new InversionException("No answer found");
            }
            //2c: Pivot if necessary
            else if (p > j) {
                this.setMatrix(this.pivot(p, j).getMatrix());
            }
            //2d: Divide row j by the pivot value
            double scalar = this.getValueAt(j, j);
            Matrix.divideByScalarDestructive(this, scalar, j);

            //2e: For each i != j, do row reduction
            Matrix tempMatrix = new Matrix(this.numRows(), this.numColumns());
            for (int y = 0; y < this.numRows(); y++){
                if (y == j){
                    for (int x = 0; x < tempMatrix.numColumns(); x++){
                        tempMatrix.setMatrix(y, x, this.getValueAt(y, x));
                    }
                }
                else {
                    for (int x = 0; x < this.numColumns(); x++){
                        double cellValue = this.getValueAt(y, x) - this.getValueAt(y, j) * this.getValueAt(j, x);
                        tempMatrix.setMatrix(y, x, cellValue);
                    }
                }
            }
            this.setMatrix(tempMatrix.getMatrix());
        }
    }

    public double determinate() throws InversionException {
        Matrix gaussianMatrix = this.copy();
        for(int j = 0; j < gaussianMatrix.numRows() - 1; j++){
            //2a: Find largest absolute value in column
            int p = gaussianMatrix.largestAbsoluteValueIndex(j);
            //2b: Check for flag, exit if answer cannot be found
            if (gaussianMatrix.getValueAt(p,j) == 0.0){
                throw new InversionException("No answer found");
            }
            //2c: Pivot if necessary
            else if (p > j) {
                gaussianMatrix.setMatrix(gaussianMatrix.pivot(p, j).getMatrix());
            }

            //2e: For each i > j, do row reduction
            Matrix tempMatrix = new Matrix(gaussianMatrix.numRows(), gaussianMatrix.numColumns());
            for (int y = 0; y < gaussianMatrix.numRows(); y++){
                if (y <= j){
                    for (int x = 0; x < tempMatrix.numColumns(); x++){
                        tempMatrix.setMatrix(y, x, gaussianMatrix.getValueAt(y, x));
                    }
                }
                else {
                    for (int x = 0; x < gaussianMatrix.numColumns(); x++){
                        double cellValue = gaussianMatrix.getValueAt(y, x) - gaussianMatrix.getValueAt(y, j) / gaussianMatrix.getValueAt(j, j) * gaussianMatrix.getValueAt(j,x);
                        tempMatrix.setMatrix(y, x, cellValue);
                    }
                }
            }
            gaussianMatrix.setMatrix(tempMatrix.getMatrix());
        }
        //calculate along diagonal
        double diagonalValue = 0.0;
        for (int i= 0; i < gaussianMatrix.numColumns(); i++){
            if (i == 0){
                diagonalValue = gaussianMatrix.getValueAt(i, i);
            } else {
                diagonalValue *= gaussianMatrix.getValueAt(i, i);
            }
        }
        //Negative 1 for delta
        if(gaussianMatrix.numRows() % 2 != 0){
            diagonalValue *= -1;
        }
        return diagonalValue;
    }

    public Matrix pivot(int pivotIndex, int rowToPivotInto){
        Matrix pivotedMatrix = new Matrix(this.numRows(), this.numColumns());
        //p=1, j=0
        for (int i = 0; i < pivotedMatrix.numRows(); i++){
            for(int x = 0; x < pivotedMatrix.numColumns(); x++){
                if(i == pivotIndex){
                    //copy rowToPivotInto's values
                    pivotedMatrix.setMatrix(i, x, this.getValueAt(rowToPivotInto, x));
                } else if (i == rowToPivotInto){
                    //copy pivotIndex row's values
                    pivotedMatrix.setMatrix(i, x, this.getValueAt(pivotIndex, x));
                } else {
                    //Copy row
                    pivotedMatrix.setMatrix(i, x, this.getValueAt(i, x));
                }
            }
        }
        return pivotedMatrix;
    }

    //Call this on a square (nxn) matrix and send the coefficient matrix as a parameter
    public void inverse(Matrix coefficientMatrix) throws Exception {
        //Create augmented matrix
        Matrix augmentedMatrix = Matrix.createAugmentedMatrix(this, coefficientMatrix);
        augmentedMatrix.inverse();
        this.setMatrix(augmentedMatrix.getMatrix());
    }

    //Determines the trace of a matrix (the sum of it's diagonal entries)
    public double trace() throws Exception {
        double sumOfDiagonals = 0.0;

        //Matrix must be square
        if (this.numColumns() != this.numRows()){
            throw new Exception("Matrix must be square");
        }

        //Enumerate through diagonals, add to sum
        for(int i = 0; i < this.numRows(); i++){
            double currentValue = this.getValueAt(i,i);
            sumOfDiagonals += currentValue;
        }

        return sumOfDiagonals;
    }

    //This method applies leverrier's method to a copy of the matrix that it is called upon,
    //Returns: A coefficient matrix with zeros in the upper and lower triangle (only has values down the middle)
    //these values represent the coefficients for the characteristic equation, starting with the highest degree
    //Throws: Exception for various reasons, such as the inability to create an identity matrix for some reason
    public Matrix leverriersMethod() throws Exception {
        Matrix characteristicPolynomial = Matrix.createIdentityMatrix(this.numRows());
        //Create a copy of this matrix that we can utilize
        Matrix A = this.copy();
        //Set Bn = A, and a(n) = -trace(Bn)
        Matrix Bn = A.copy();
        double an = -1 * Bn.trace();

        //For k=n-1 down to 1, compute:
        int n = this.numRows();
        //assign value to characteristic polynomial (off by one issue)
        int count = 0;
        characteristicPolynomial.setMatrix(count, count, an);
        count++;
        for(int k = n-1; k > 0; k--){
            //Bk = A(Bk+1 + a(k+1)*I)
            Matrix identityMatrix = Matrix.createIdentityMatrix(this.numRows());
            Matrix scaledMatrix = identityMatrix.multiplyByScalar(an);
            Matrix combinedMatrix = Matrix.add(Bn, scaledMatrix);
            Bn = Matrix.multiply(A, combinedMatrix);

            //ak = -1*(trace(Bk))/(n-k+1)
            double trace = Bn.trace();
            an = -1 * trace / (n - k + 1);

            //Assign value to characteristic polynomial
            characteristicPolynomial.setMatrix(count, count, an);
            count++;
        }
        return characteristicPolynomial;
    }

    //Purpose: Applies the power method to a Matrix
    //Returns: An array of Matrices where Matrix[0] = 1x1 matrix containing the estimated eigenvalue
    // and Matrix[1] = the associated eigenvector
    //Parameters:
    //  sigma: How close the residual error must be to zero in order for the function to return
    //  maxIterations: how many loops the program will run through before returning it's best result (if sigma is not met)
    //This method is NON-Destructive to the original array
    public Matrix[] powerMethod(double sigma, int maxIterations) throws Exception {
        Matrix[] answers = new Matrix[2];

        Matrix y = new Matrix(this.numRows(), 1);
        for(int i = 0; i < y.numRows(); i++){
            if(i%2 == 0) {
                y.setMatrix(i, 0, 0.0);
            } else {
                y.setMatrix(i, 0, 1.0);
            }
        }
        Matrix x = Matrix.multiply(this, y);
        Matrix r = new Matrix(1,1);
        Matrix mu;
        int k = 0;
        do {
            Matrix xDoubleBar = new Matrix(new double[][]{
                    {x.maximumAbsoluteValue()}
            }); //||x||
            y = Matrix.divide(x, xDoubleBar);
            x = Matrix.multiply(this, y); //x=Ay
            Matrix temp1 = Matrix.multiply(Matrix.transpose(y), x);
            Matrix temp2 = Matrix.multiply(Matrix.transpose(y), y);
            mu = Matrix.divide(temp1, temp2);
            temp1 = Matrix.multiply(y, mu);
            r = Matrix.subtract(temp1, x);

            k++;
        } while(r.maximumAbsoluteValue() > sigma && k < maxIterations);
        System.out.println(k);
        answers[0] = mu;
        answers[1] = y;
        return answers;
    }

    //Purpose: Given a matrix, apply jacobi's method to it (non-destructive)
    //Returns: Array of matrices where:
      //Matrix[0] = A permutation of A which has zeros everywhere except the diagonals,
        //and the diagonals are the eigenvalues
      //Matrix[1] = P, whose columns are the corresponding eigenvectors
    //Parameters:
      //sigma: The maximum magnitude allowed above the diagonal. Smaller values lead to more accurate results
    public Matrix[] jacobi(double sigma) throws Exception {
        if(this.numColumns() != this.numRows()) {
            return null;
        }
        Matrix A = this.copy();
        Matrix P = Matrix.createIdentityMatrix(A.numRows());
        Matrix R;
        Vector largestIndex = A.largestMagnitudeAboveDiagonal();
        int p = (int) largestIndex.getX();
        int q = (int) largestIndex.getY();
        double Apq = A.getValueAt(p, q);

        while(Math.abs(Apq) > sigma){
            double arctanNumerator = 2 * Apq;
            double arctanDenominator = A.getValueAt(p,p) - A.getValueAt(q,q);
            double arctanAnswer = Math.atan(arctanNumerator/arctanDenominator);
            double angle = 0.5 * arctanAnswer;

            R = Matrix.createIdentityMatrix(A.numRows());
            R.setMatrix(p, p, Math.cos(angle));
            R.setMatrix(q, q, Math.cos(angle));
            R.setMatrix(p, q, -1*Math.sin(angle));
            R.setMatrix(q, p, Math.sin(angle));

            P = Matrix.multiply(P, R);
            Matrix RTranspose = Matrix.transpose(R);
            A = Matrix.multiply(Matrix.multiply(RTranspose, A), R);

            largestIndex = A.largestMagnitudeAboveDiagonal();
            p = (int) largestIndex.getX();
            q = (int) largestIndex.getY();
            Apq = A.getValueAt(p, q);
        }
        Matrix[] results = new Matrix[2];
        results[0] = A;
        results[1] = P;
        return results;
    }

    //Class Methods
    public static Matrix createIdentityMatrix(int numberOfRowsAndColumns){
        Matrix identityMatrix = new Matrix(numberOfRowsAndColumns, numberOfRowsAndColumns);
        for(int i = 0; i < identityMatrix.numRows(); i++){
            for (int j = 0; j < identityMatrix.numColumns(); j++){
                if (i == j){
                    identityMatrix.setMatrix(i, j, 1.0);
                } else {
                    identityMatrix.setMatrix(i, j, 0.0);
                }
            }
        }
        return identityMatrix;
    }

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

    public static Matrix subtract(Matrix matrixOne, Matrix matrixTwo){
        //divide MatrixTwo by -1
        Matrix.divideByScalarDestructive(matrixTwo, -1.0);
        //Add the two together
        return Matrix.add(matrixOne, matrixTwo);
    }
    //Warning, this is a destructive method.
    public static void divideByScalarDestructive(Matrix matrix, double scalar){
        //enumerate through each entry in matrix table, divide by the scalar and set that new value into the matrix
        int numRows = matrix.getMatrix().length;
        int numColumns = matrix.getMatrix()[0].length;
        for (int i = 0; i < numRows; i++){
            for (int j = 0; j < numColumns; j++){
                double matrixValue = matrix.getMatrix()[i][j];
                matrix.setMatrix(i, j, matrixValue / scalar);
            }
        }
    }

    public static Matrix divide(Matrix numerator, Matrix denominator) throws Exception {
        Matrix inverted = new Matrix(denominator.numRows(), denominator.numColumns());
        //Set all values in 'inverted' = 1/x, where x = the value in that cell in 'denominator'
        for(int i = 0; i < inverted.numRows(); i++){
            for(int j = 0; j < inverted.numColumns(); j++){
                inverted.setMatrix(i,j, 1.0/denominator.getValueAt(i,j));
            }
        }
        return Matrix.multiply(numerator, inverted);
    }

    //Warning, this is a destructive method.
    //This method is used for Gauss-Jordan elimination
    public static void divideByScalarDestructive(Matrix matrix, double scalar, int rowIndex){
        for (int columnIndex = 0; columnIndex < matrix.numColumns(); columnIndex++){
            double matrixValue = matrix.getValueAt(rowIndex, columnIndex);
            matrix.setMatrix(rowIndex, columnIndex, matrixValue / scalar);
        }
    }

    public Matrix multiplyByScalar(double scalar){
        Matrix multipliedMatrix = this.copy();
        for (int i = 0; i < multipliedMatrix.numRows(); i++){
            for(int j = 0; j < multipliedMatrix.numColumns(); j++){
                multipliedMatrix.setMatrix(i,j,multipliedMatrix.getValueAt(i,j) * scalar);
            }
        }
        return multipliedMatrix;
    }

    public static Matrix createAugmentedMatrix(Matrix squareMatrix, Matrix coefficientMatrix) throws Exception {
        //Determine what length the matrix will be / verify appropriate length
        if (squareMatrix.numColumns() == squareMatrix.numRows() && squareMatrix.numRows() == coefficientMatrix.numRows()){
            //create a new matrix with the appropriate size
            Matrix augmentedMatrix = new Matrix(new double[squareMatrix.numRows()][squareMatrix.numColumns() + coefficientMatrix.numColumns()]);

            //Set matrix values
            //For each row:
            for (int i = 0; i < augmentedMatrix.numRows(); i++){

                for (int j = 0; j < augmentedMatrix.numColumns(); j++){
                    //First n columns (where n = squareMatrix.numColumns) will be square matrix's values
                    if (j < squareMatrix.numColumns()){
                        double value = squareMatrix.getValueAt(i,j);
                        augmentedMatrix.setMatrix(i, j, value);
                    }
                    else {
                        //Last columns will be coefficientMatrix's value
                        augmentedMatrix.setMatrix(i,j, coefficientMatrix.getValueAt(i, j - squareMatrix.numColumns()));
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
        //multiply by inverses
        System.out.println("");
        Matrix[] twoByTwoMatrices = multiplyByTranspose(matrices);
        //Calculate matrix mean
        Matrix matrixMean = matrixMean(twoByTwoMatrices);
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

    public static Matrix transpose(Matrix originalMatrix){
        Matrix transposedMatrix = new Matrix(originalMatrix.numColumns(), originalMatrix.numRows());
        for (int i = 0; i < transposedMatrix.numRows(); i++){
            for (int j = 0; j < transposedMatrix.numColumns(); j++){
                transposedMatrix.setMatrix(i, j, originalMatrix.getValueAt(j, i));
            }
        }

        return transposedMatrix;

//        old method that only worked on a 1x2 matrix
//        double[][] transposedTable = new double[2][1];
//        transposedTable[0][0] = matrix.getMatrix()[0][0];
//        transposedTable[1][0] = matrix.getMatrix()[0][1];
//        return new Matrix(transposedTable);
    }

    public static Matrix multiply(Matrix matrixOne, Matrix matrixTwo) throws Exception {
        Matrix multipliedMatrix;

        //Assuming a matrixOne is Mx*y, and matrixTwo is Mi*j, then they can be multiplied if y==i
        //In other words, number of columns in matrixOne must be equal to the number of rows in matrixTwo
        if(matrixOne.numColumns() == matrixTwo.numRows()) {
            int numberOfMultiplicationsPerCell = matrixOne.numColumns();

            //Assuming a matrixOne is Mx*y, and matrixTwo is Mi*j, then the resulting matrix will be Mx*j
            int numberOfColumnsForAnswer = matrixTwo.numColumns();
            int numberOfRowsForAnswer = matrixOne.numRows();
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
    //Takes an array of doubles and returns the largest value
    public static double maximum(double[] sums){
        double largest = -99999999;
        for(double sum: sums){
            if (sum > largest){
                largest = sum;
            }
        }
        return largest;
    }
    //Returns the largest value in this matrix
    public double maximumAbsoluteValue(){
        double max = -9999999;
        for(int i = 0; i < this.numRows(); i++){
            for(int j = 0; j < this.numColumns(); j++){
                double currentValue = this.getValueAt(i, j);
                if(currentValue < 0){
                    currentValue *= -1;
                }
                if (currentValue > max){
                    max = currentValue;
                }
            }
        }
        return max;
    }

    public Vector largestMagnitudeAboveDiagonal(){
        ArrayList<Vector> indices = new ArrayList<>();
        Vector result;

        for(int currentRowIndex = 0; currentRowIndex < this.numRows(); currentRowIndex++){
            for(int currentColumnIndex = 0; currentColumnIndex < this.numColumns(); currentColumnIndex++){
                if(currentColumnIndex > currentRowIndex) { //Above the diagonal
                    indices.add(new Vector(currentRowIndex, currentColumnIndex)); //Add the current value
                }
            }
        }

        result = findLargestMagnitude(this, indices);

        return result;
    }

    private Vector findLargestMagnitude(Matrix matrix, ArrayList<Vector> indices){
        Vector largestVector = null;
        double largestMagnitude = -1;
        for(Vector currentVector: indices){
            double magnitude = matrix.getValueAt((int) currentVector.getX(), (int) currentVector.getY());
            if (magnitude < 0) {magnitude *= -1;}
            if (magnitude > largestMagnitude){
                largestMagnitude = magnitude;
                largestVector = currentVector;
            }
        }
        return largestVector;
    }
}
