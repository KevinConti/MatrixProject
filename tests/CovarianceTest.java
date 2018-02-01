import org.junit.Test;
import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class CovarianceTest {
    @Test
    public void givenSampleDataTheToCovarianceMethodReturnsTheProperResult(){
        Vector[][] vectors = initializeVectorsFromFile("data/fake_data.txt");
        Matrix classOneCovariance = Matrix.toCovariance(vectors[0]);
        System.out.println("");
        System.out.println("Class One Covariance:");
        System.out.println(classOneCovariance);
        Matrix classTwoCovariance = Matrix.toCovariance(vectors[1]);
        System.out.println("");
        System.out.println("Class Two Covariance:");
        System.out.println(classTwoCovariance);
    }

    private Vector[][] initializeVectorsFromFile(String filepath) {
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }
}
