public class main {
    public static void main(String[] args){
        Vector[][] vectors = initializeVectorsFromFile("data/data.txt");
    }

    private static Vector[][] initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }
}
