import java.util.ArrayList;

public class Vector {
    private double x;
    private double y;

    //Returns a Vector with x and y values which are the mean of the ArrayList of Vectors sent in params
    public static Vector add(ArrayList<Vector> vectors){

        double xSum = 0.0;
        double ySum = 0.0;

        double xMean = 0.0;
        double yMean = 0.0;

        for(int i = 0; i < vectors.size() ; i++){
            Vector vectorAtCurrentIndex = vectors.get(i);
            xSum = xSum + vectorAtCurrentIndex.getX();
            ySum = ySum + vectorAtCurrentIndex.getY();
        }

        xMean = xSum / vectors.size();
        yMean = ySum / vectors.size();

        return new Vector(xMean, yMean);
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
