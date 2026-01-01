package heparansulfate;

import java.util.Random;

/**
 * Projection of a vector {@code y} onto a polyhedral set ({@code Ax <= b}) via coordinate
 * descent. Solves the optimization problem: minimize {@code ||x - y||^2} with respect to 
 * {@code x} subject to the constraints {@code Ax <= b}.
 */
public class ProjOnPolyHSet {
    /**
     * Convergence certificate. Set to {@code false} if the coordinate descent 
     * process failed to converge within the allowed iterations.
     */
    boolean success = false;
    /**
     * The resulting optimal vector {@code x} representing the projection of {@code y} 
     * onto the polyhedral set {@code Ax <= b}.
     */
    public double[] optimum = null;

    /**
     * Projects vector {@code y} onto a polyhedral set defined by the linear 
     * inequalities {@code Ax <= b}.
     * 
     * @param y The vector to be projected.
     * @param A The matrix defining the inequality constraints {@code Ax <= b}.
     * @param b The vector of targets in the inequality constraints {@code Ax <= b}.
     * @param rand Random number generator used to seed the initial coordinate descent state.
     */
    public ProjOnPolyHSet(double[] y, double[][] A, double[] b, Random rand) {
        double[][] AT = MatrixOp.transpose(A);
        double[][] P = MatrixOp.multMat(A, AT);
        for (int i = 0; i < P.length; i++) {
            for (int j = 0; j < P[i].length; j++) {
                P[i][j] *= 0.5;
            }
        }
        double[] t = MatrixOp.multMatVec(A, y);
        for (int j = 0; j < t.length; j++) {
            t[j] = b[j] - t[j];
        }
        double[] s = new double[P.length];
        for (int j = 0; j < s.length; j++) {
            s[j] = rand.nextDouble();
        }
        CoordDescNNQP cd = new CoordDescNNQP(P, t, s);
        success = cd.success;
        optimum = MatrixOp.multMatVec(AT, cd.optimum);
        for (int j = 0; j < optimum.length; j++) {
            optimum[j] = y[j] - 0.5 * optimum[j];
        }
    }
}