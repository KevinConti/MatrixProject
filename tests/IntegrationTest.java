import org.junit.Test;
import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class IntegrationTest {
    @Test
    public void givenSampleDataTheProgramRunsAsExpected(){
        Vector[][] vectors = initializeVectorsFromFile("data/fake_data.txt");
        System.out.println("Class One:");
        Matrix classOneCovariance = Matrix.toCovariance(vectors[0]);
        System.out.println("Class One Covariance:");
        System.out.println(classOneCovariance);
        System.out.println("Class Two:");
        Matrix classTwoCovariance = Matrix.toCovariance(vectors[1]);
        System.out.println("Class Two Covariance:");
        System.out.println(classTwoCovariance);

        //Next, test inverses
        Matrix identityMatrix = Matrix.createIdentityMatrix(classOneCovariance.numRows());
        System.out.println("");
        System.out.println("Class one inverse matrix:");
        try {
            System.out.println(classOneCovariance.copy().inverse(identityMatrix));
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("");
        System.out.println("Class two inverse matrix:");
        try {
            System.out.println(classTwoCovariance.copy().inverse(identityMatrix));
        } catch (Exception e) {
            e.printStackTrace();
        }

        //Next, find determinants
        System.out.println("Class one determinant");
        try {
            System.out.println(Double.toString(classOneCovariance.determinate()));
        } catch (InversionException e) {
            e.printStackTrace();
        }
        System.out.println("Class Two determinant");
        try {
            System.out.println(Double.toString(classTwoCovariance.determinate()));
        } catch (InversionException e) {
            e.printStackTrace();
        }

    }

    private Vector[][] initializeVectorsFromFile(String filepath) {
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }

}
