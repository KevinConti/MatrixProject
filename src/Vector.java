import java.util.ArrayList;

public class Vector {
    private double x;
    private double y;

    //Returns a Vector with x and y values which are the mean of the array of Vectors sent in params
    public static Vector mean(Vector[] vectors){

        double xSum = 0.0;
        double ySum = 0.0;

        double xMean = 0.0;
        double yMean = 0.0;

        for(int i = 0; i < vectors.length ; i++){
            Vector vectorAtCurrentIndex = vectors[i];
            xSum = xSum + vectorAtCurrentIndex.getX();
            ySum = ySum + vectorAtCurrentIndex.getY();
        }

        xMean = xSum / vectors.length;
        yMean = ySum / vectors.length;

        return new Vector(xMean, yMean);
    }

    public static Vector subtract(Vector minuend, Vector subtrahend){
        double differenceX = minuend.getX() - subtrahend.getX();
        double differenceY = minuend.getY() - subtrahend.getY();
        return new Vector(differenceX, differenceY);
    }

    //The below method subtracts all the vectors in the vector array by the vector called "subtrahend". Useful for subtracting the mean
    public static Vector[] subtract(Vector[] minuends, Vector subtrahend){
        Vector[] subtractedVectors = new Vector[minuends.length];
        for(int i = 0; i < minuends.length; i++){
            subtractedVectors[i] = Vector.subtract(minuends[i], subtrahend);
        }

        return subtractedVectors;
    }

    public Vector(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    @Override
    public String toString() {
        return "Vector{" +
                "x=" + x +
                ", y=" + y +
                '}';
    }
}
