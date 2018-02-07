import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import static javax.script.ScriptEngine.FILENAME;

//The specific implementation which takes the input from Dr. Tag's file and converts it into Vectors
//This is a specific implementation and switch cases must be made if a different form of input is to be implemented
//The return is a 2d array of vectors, where each row is a class of vectors
public class InputParser {
    public InputParser() {
        super();
    }

    public Vector[][] convertFile(String filename) {
        Vector[][] allVectors;
        ArrayList<Vector> classOneVectors = new ArrayList<>();
        ArrayList<Vector> classTwoVectors = new ArrayList<>();

        final String FILENAME = filename;

        try (BufferedReader br = new BufferedReader(new FileReader(FILENAME))) {

            String currentLine;
            ArrayList<String> values = new ArrayList<String>();

            while ((currentLine = br.readLine()) != null) {
                String[] currentLineValues = currentLine.split("\t");
                for (int i = 0; i < currentLineValues.length; i++) {
                    values.add(currentLineValues[i]);
                }
            }

            for (int i = 0; i < values.size(); i += 4) {
                double classOneXValue = Double.parseDouble(values.get(i));
                double classOneYValue = Double.parseDouble(values.get(i + 1));
                double classTwoXValue = Double.parseDouble(values.get(i + 2));
                double classTwoYValue = Double.parseDouble(values.get(i + 3));

                Vector classOneVector = new Vector(classOneXValue, classOneYValue);
                classOneVectors.add(classOneVector);
                Vector classTwoVector = new Vector(classTwoXValue, classTwoYValue);
                classTwoVectors.add(classTwoVector);
            }

            allVectors = new Vector[2][classOneVectors.size()];
            for (int j = 0; j < classOneVectors.size(); j++) {
                allVectors[0][j] = classOneVectors.get(j);
                allVectors[1][j] = classTwoVectors.get(j);
            }


        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("Error in convertFile");
            allVectors = new Vector[1][1];
        }

        return allVectors;
    }

    public Matrix convertLinear(String filename) {
        Matrix linearMatrix = null;
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {

            String currentLine;
            ArrayList<String> values = new ArrayList<String>();
            boolean firstPass = true;
            int numColumns = 0;

            while ((currentLine = br.readLine()) != null) {
                String[] currentLineValues = currentLine.split("\t");
                values.addAll(Arrays.asList(currentLineValues));
                if (firstPass){
                    numColumns = values.size();
                    firstPass = false;
                }
            }
            int count = 0;
            linearMatrix = new Matrix(values.size() / numColumns, numColumns);
            for(int i = 0; i < linearMatrix.numRows(); i++){
                for(int j = 0; j < linearMatrix.numColumns(); j++){
                    linearMatrix.setMatrix(i,j,Double.parseDouble(values.get(count)));
                    count++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return linearMatrix;
    }
}
