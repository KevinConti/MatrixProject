public class main {
    public static void main(String[] args){
        Vector[][] vectors = initializeVectorsFromFile("data/data.txt");
        System.out.println("Class One:");
        Matrix classOneCovariance = Matrix.toCovariance(vectors[0]);
        System.out.println("");
        System.out.println("Class One Covariance:");
        System.out.println(classOneCovariance);
        System.out.println("Class Two:");
        Matrix classTwoCovariance = Matrix.toCovariance(vectors[1]);
        System.out.println("");
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

//        Vector classOneMeanVector = Vector.mean(vectors[0]);
//        Matrix classOneInverse = classOneCovariance.copy();
//        try {
//            classOneCovariance.inverse();
//        } catch (InversionException e) {
//            e.printStackTrace();
//        }
//        Discriminate classOneDiscriminate = new Discriminate();
//        try {
//            classOneDiscriminate = new Discriminate(classOneMeanVector, classOneInverse, classOneCovariance.copy().determinate());
//        } catch (InversionException e) {
//            e.printStackTrace();
//
//        }
//        try {
//            classOneDiscriminate.classify(2.125675004181818, 3.1825644208181822);
//        } catch (Exception e) {
//            e.printStackTrace();
//        }

    }

    private static Vector[][] initializeVectorsFromFile(String filepath){
        InputParser ip = new InputParser();
        return ip.convertFile(filepath);
    }
}
