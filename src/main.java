import java.util.ArrayList;

public class main {
    public static void main(String[] args){
        InputParser ip = new InputParser();
        String filepath = "data/data.txt";
        Vector[][] vectors = ip.convertFile(filepath);

        for (Vector vector : vectors[0]){
            System.out.println(vector);
        }
    }
}
