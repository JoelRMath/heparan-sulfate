package heparansulfate;

/**
 * Solution to nonnegative quadratic programming:
 * minimize {@code (1/2)x^TQx + c^Tx} subject to {@code x >= 0}
 * via coordinate descent. Note that Q must be symmetric
 * and positive definite
 */
public class CoordDescNNQP {
    /**
     * Dimension of x (Q is n*n)
     */
    int n = 0;
    /**
     * Quadratic coefficients in the cost {@code (1/2)x^TQx + c^Tx}
     */
    double[][] Q = null;
    /**
     * Linear coefficients in the cost {@code (1/2)x^TQx + c^Tx}
     */
    double[] c = null;
    /**
     * Point being iterated
     */
    double[] x = null;
    /**
     * Keeps track of the number of iterations (1 iter == change of
     * all components x_1, x_2, ..., x_n)
     */
    public int iter = 0;
    /**
     * Maximum number of iterations (1e6); if it goes above,
     * matrix Q is likely ill-conditioned
     */
    int maxIter = 1000000;
    /**
     * Final value of the objective function {@code (1/2)x^TQx + c^Tx}
     */
    double finalCost = 0.;
    /**
     * Optimal value of {@code x}.
     */
    public double[] optimum = null;
    /**
     * Used to stop iterations
     */
    double eps = 1e-10;
    /**
     * Certificate
     */
    boolean success = false;

    /**
     * Solves: minimize {@code (1/2)x^TQx + c^Tx} subject to {@code x >= 0}.
     * Solution x is in this.optimum
     * @param Q coefficients of the quadratic terms in the cost function.
     * @param c coefficients of the linear terms in the cost function.
     * @param y starting point (any random point in the positive orthant).
     */
    public CoordDescNNQP(double[][] Q, double[] c, double[] y) {
        // initializations
        n = Q.length;
        this.Q = Q;
        this.c = c;
        x = new double[n];
        for (int j = 0; j < n; j++) {
            x[j] = y[j];
        }
        double alpha = 1.;
        iter = 0;
        double cost = getCost(x);
        boolean b = true;
        // coordinate descent method
        while (b) {
            iter++;
            alpha = 0.;
            for (int j = 0; j < n; j++) {
                double t = c[j];
                for (int k = 0; k < n; k++) {
                    t += (Q[j][k] * x[k]);
                }
                t /= Q[j][j];
                t *= -1.;
                t += x[j];
                if (t >= 0.) {
                    alpha += ((x[j] - t) * (x[j] - t));
                    x[j] = t;
                } else {
                    alpha += (x[j] * x[j]);
                    x[j] = 0.;
                }
            }
            alpha = Math.sqrt(alpha);
            alpha /= (double) x.length;
            cost = getCost(x);
            if (alpha < eps) {
                b = false;
                success = true;
            }
            if (b) {
                if (iter > maxIter) {
                    success = false;
                    b = false;
                }
            }
        }
        // solution and final cost
        optimum = new double[n];
        for (int j = 0; j < n; j++) {
            optimum[j] = x[j];
        }
        finalCost = cost;
    }

    /**
     * Computes the objective function value: {@code (1/2)z^TQz + c^Tz}
     * @param z vector z
     * @return {@code (1/2)z^TQz + c^Tz}
     */
    double getCost(double[] z) {
        double res = 0.;
        res = MatrixOp.geInnerProd(MatrixOp.multMatVec(Q, z), z);
        res /= 2.;
        res += MatrixOp.geInnerProd(c, z);
        return res;
    }

    /**
     * Returns the Euclidean norm of a vector: {@code z^T z}.
     * @param z Input vector.
     * @return {@code z^T z}.
     */
    double getNorm2(double[] z) {
        double res = 0.;
        for (int j = 0; j < z.length; j++) {
            res += (z[j] * z[j]);
        }
        return res;
    }
}