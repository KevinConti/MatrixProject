public class main {
    public static void main(String[] args){
        Vector[][] vectors = initializeVectorsFromFile("data/data.txt");
        Matrix classOneCovariance = Matrix.toCovariance(vectors[0]);
        System.out.println("");
        System.out.println("Class One Covariance:");
        System.out.println(classOneCovariance);
        Matrix classTwoCovariance = Matrix.toCovariance(vectors[1]);
        System.out.println("");
        System.out.println("Class Two Covariance:");
        System.out.println(classTwoCovariance);
    }

    private static Vector[][] initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }
}
