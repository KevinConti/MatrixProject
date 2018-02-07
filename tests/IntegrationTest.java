import org.junit.Test;
import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

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
        Matrix classOneInverse = classOneCovariance.copy();
        try {
            classOneInverse.inverse(identityMatrix);
            classOneInverse.removeFirstColumn();
            classOneInverse.removeFirstColumn();
            System.out.println(classOneInverse);
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        System.out.println("");
        System.out.println("Class two inverse matrix:");
        Matrix classTwoInverse = classTwoCovariance.copy();
        try {
            classTwoInverse.inverse(identityMatrix);
            //Remove first two columns to eliminate identity matrix
            classTwoInverse.removeFirstColumn();
            classTwoInverse.removeFirstColumn();
            System.out.println(classTwoInverse);
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        //Next, find determinants
        System.out.println("Class one determinant");
        double classOneDeterminate = 0.0;
        try {
            classOneDeterminate = classOneCovariance.determinate();
            System.out.println(Double.toString(classOneDeterminate));
        } catch (InversionException e) {
            e.printStackTrace();
            fail();
        }
        System.out.println("Class Two determinant");
        double classTwoDeterminate = 0.0;
        try {
            classTwoDeterminate = classTwoCovariance.determinate();
            System.out.println(Double.toString(classTwoDeterminate));
        } catch (InversionException e) {
            e.printStackTrace();
            fail();
        }

    }

    private Vector[][] initializeVectorsFromFile(String filepath) {
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }

}
