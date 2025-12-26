package heparansulfate;

import java.util.Random;

/**
 * projection of a vector y onto a polyhedral set ({@code Ax <= b}) via coordinate
 * descent: minimize w.r.t. x {@code ||x-y||^2} subject to {@code Ax <= b}
 */
public class ProjOnPolyHSet {
    /**
     * certificate (false if coordinate descent was too slow)
     */
    boolean success = false;
    /**
     * projection of y onto polyhedral set {@code Ax <= b}
     */
    public double[] optimum = null;

    /**
     * projection of a vector y onto a polyhedral set ({@code Ax <= b}) via coordinate
     * descent: minimize w.r.t. x {@code ||x-y||^2} subject to {@code Ax <= b}
     * @param y vector to project
     * @param A matrix in inequality constraints {@code Ax <= b}
     * @param b vector in inequality constraints {@code Ax <= b}
     * @param rand used to seed the coordinate descent
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