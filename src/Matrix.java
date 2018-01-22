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

    @Override
    public String toString() {
        StringBuilder fullArray = new StringBuilder();
        for(int row = 0; row < matrix.length; row++){
            fullArray.append(Arrays.toString(matrix[row]));
        }
        return fullArray.toString();
    }
}
