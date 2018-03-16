public class ProblemTwo {
    public static void main(String[] args){
        Matrix originalProblem = new Matrix(new double[][]{
                {0, 0, 0, 0, 0, 756/30},
                {1, 0, 0, 0, 0, 2733.0/30.0},
                {0, 1, 0, 0, 0, -4903.0/30.0},
                {0, 0, 1, 0, 0, 1689.0/30},
                {0, 0, 0, 1, 0, 139.0/30.0},
                {0, 0, 0, 0, 1, 1}
        });
        runIteration(originalProblem);

        Matrix newIteration = new Matrix(new double[][]{
                {0, 0, 0, 0, 108.0/30.0},
                {1, 0, 0, 0, 375.0/30.0},
                {0, 1, 0, 0, -754.0/30.0},
                {0, 0, 1, 0, 349.0/30.0},
                {0, 0, 0, 1, 1.0}
        });
        runIteration(newIteration);

        newIteration = new Matrix((new double[][]{
                {0, 0, 0, -12.0/30.0},
                {1, 0, 0, -43.0/30.0},
                {0, 1, 0, 79.0/30.0},
                {0, 0, 1, -1}
        }));
        runIteration(newIteration);

        newIteration = new Matrix(new double[][]{
                {0, 0, 4.0/15.0},
                {1, 0, 17.0/15.0},
                {0, 1, 15.0/15.0}
        });
        runIteration(newIteration);

        

    }

    private static void runIteration(Matrix matrix){
        try {
            Matrix solution = matrix.leverriersMethod();
            System.out.println("Coefficient Matrix:");
            System.out.println(solution);
            Matrix otherSolution[] = solution.powerMethod(.00001, 100000);
            System.out.println("Largest Eigenvalue");
            System.out.println(otherSolution[0]);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
