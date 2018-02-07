public class Discriminate {
    private Vector meanVector;
    private Matrix inverseMatrix;
    private double determinate;

    public Discriminate() {
        super();
    }

    public Discriminate(Vector meanVector, Matrix inverseMatrix, double determinate) {
        this.meanVector = meanVector;
        this.inverseMatrix = inverseMatrix;
        this.determinate = determinate;
    }

    public Vector getMeanVector() {
        return meanVector;
    }

    public void setMeanVector(Vector meanVector) {
        this.meanVector = meanVector;
    }

    public Matrix getInverseMatrix() {
        return inverseMatrix;
    }

    public void setInverseMatrix(Matrix inverseMatrix) {
        this.inverseMatrix = inverseMatrix;
    }

    public double getDeterminate() {
        return determinate;
    }

    public void setDeterminate(double determinate) {
        this.determinate = determinate;
    }

    public double classify(double x, double y) throws Exception {
        double result = 0.0;
        //Create x-mu, y-mu 1x2 matrix
        Matrix oneByTwoMeanMatrix = createOneByTwoMatrix(x, y);
        //Multiply meanMatrix by inverseMatrix
        Matrix tempMatrixOne = Matrix.multiply(oneByTwoMeanMatrix, inverseMatrix);
        //Create x-mu, y-mu 2x1 matrix
        Matrix twoByOneMeanMatrix = createTwoByOneMatrix(x, y);
        //multiply resultMatrix by 2x1
        Matrix tempMatrixTwo = Matrix.multiply(tempMatrixOne, twoByOneMeanMatrix);
        //Convert the single cell matrix to a double
        result = tempMatrixTwo.getValueAt(0,0);
        //Multiply by -0.5
        result *= -0.5;
        //Subtract by 1/2*ln(determinate)
        result -= .5 * Math.log(this.determinate);
        //???Add ln(1/2)
        return result;
    }
    private Matrix createOneByTwoMatrix(double x, double y){
        //Creates a matrix of the following form:
        //[x-mu y-mu2]
        double xMinusMuOne = x - this.getMeanVector().getX();
        double yMinusMuTwo = y - this.getMeanVector().getY();

        Matrix oneByTwo = new Matrix(1,2);
        oneByTwo.setMatrix(0,0, xMinusMuOne);
        oneByTwo.setMatrix(0,1, yMinusMuTwo);

        return oneByTwo;
    }

    private Matrix createTwoByOneMatrix(double x, double y){
        //Creates a matrix of the following form:
        //[x-mu ]
        //[y-mu2]
        double xMinusMuOne = x - this.getMeanVector().getX();
        double yMinusMuTwo = y - this.getMeanVector().getY();

        Matrix twoByOne = new Matrix(2,1);
        twoByOne.setMatrix(0,0, xMinusMuOne);
        twoByOne.setMatrix(1,0, yMinusMuTwo);

        return twoByOne;
    }
}
