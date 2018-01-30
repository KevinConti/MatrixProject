public class main {
    public static void main(String[] args){
        Vector[][] vectors = initializeVectorsFromFile("data/data.txt");
        Matrix classOneCovariance = Matrix.toCovariance(vectors[0]);
        Matrix classTwoCovariance = Matrix.toCovariance(vectors[1]);
    }

    private static Vector[][] initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }
}
