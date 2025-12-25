package heparansulfate;

/**
 * solution to nonnegative quadratic programming:
 * minimize (1/2)x^TQx + c^Tx subject to x >= 0
 * via coordinate descent. Note that Q must be symmetric
 * and positive definite
 */
public class CoordDescNNQP {
    /**
     * dimension of x (Q is n*n)
     */
    int n = 0;
    /**
     * quadratic coefficients in the cost (1/2)x^TQx + c^Tx
     */
    double[][] Q = null;
    /**
     * linear coefficients in the cost (1/2)x^TQx + c^Tx
     */
    double[] c = null;
    /**
     * point being iterated
     */
    double[] x = null;
    /**
     * keeps track of the number of iterations (1 iter == change of
     * all components x_1, x_2, ..., x_n)
     */
    public int iter = 0;
    /**
     * maximum number of iterations (1e6); if it goes above,
     * matrix Q is likely ill-conditioned
     */
    int maxIter = 1000000;
    /**
     * final value of the objective function (1/2)x^TQx + c^Tx
     */
    double finalCost = 0.;
    /**
     * solution
     */
    public double[] optimum = null;
    /**
     * used to stop iterations
     */
    double eps = 1e-10;
    /**
     * certificate
     */
    boolean success = false;

    /**
     * Solves: minimize (1/2)x^TQx + c^Tx subject to x >= 0.
     * Solution x is in this.optimum
     * @param Q coefficients of the quadratic terms in the cost function
     * @param c coefficients of the linear terms in the cost function
     * @param y starting point (any random point in the positive orthant)
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
     * computes (1/2)z^TQz + c^Tz
     * @param z
     * @return (1/2)z^TQz + c^Tz
     */
    double getCost(double[] z) {
        double res = 0.;
        res = MatrixOp.geInnerProd(MatrixOp.multMatVec(Q, z), z);
        res /= 2.;
        res += MatrixOp.geInnerProd(c, z);
        return res;
    }

    /**
     * returns z^Tz
     * @param z vector
     * @return z^Tz
     */
    double getNorm2(double[] z) {
        double res = 0.;
        for (int j = 0; j < z.length; j++) {
            res += (z[j] * z[j]);
        }
        return res;
    }
}