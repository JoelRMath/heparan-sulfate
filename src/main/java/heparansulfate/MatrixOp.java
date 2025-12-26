package heparansulfate;

/**
 * This class contains a few static methods for matrix and vector operations.
 */
public class MatrixOp {

    /**
     * Default constructor
     */
    public MatrixOp() {
    }

    /**
     * Prints out a matrix
     * @param mat matrix
     */
    public static void printMat(double[][] mat) {
        String s = null;
        System.out.println("*************");
        for (int i = 0; i < mat.length; i++) {
            s = "";
            for (int j = 0; j < mat[0].length; j++) {
                s += String.valueOf(mat[i][j]);
                if (j != mat[0].length - 1)
                    s += "\t";
            }
            System.out.println(s);
        }
        System.out.println("*************");
    }

    /**
     * Prints out a matrix, with "." for 0 and "x" for entries not 0
     * @param mat matrix
     * @param eps definition of 0
     */
    public static void printMatS(double[][] mat, double eps) {
        String s = null;
        System.out.println("*************");
        for (int i = 0; i < mat.length; i++) {
            s = "";
            for (int j = 0; j < mat[0].length; j++) {
                String e = "x";
                if (Math.abs(mat[i][j]) <= eps) {
                    e = ".";
                }
                s += e;
                if (j != mat[0].length - 1)
                    s += "\t";
            }
            System.out.println(s);
        }
        System.out.println("*************");
    }

    /**
     * Returns a String representation of a matrix
     * @param mat matrix
     * @return a String representation of matrix ’mat’
     */
    public static String getMatStringRep(double[][] mat) {
        String res = "*************\n";
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                res += String.valueOf(mat[i][j]);
                if (j != mat[0].length - 1)
                    res += "\t";
            }
            res += "\n";
        }
        res += "*************\n";
        return res;
    }

    /**
     * Returns a String representation of a vector
     * @param x vector
     * @return String representation of vector ’x’
     */
    public static String getVecStringRep(double[] x) {
        String res = "*************\n";
        for (int i = 0; i < x.length; i++) {
            res += String.valueOf(x[i]) + "\n";
        }
        res += "*************\n";
        return res;
    }

    /**
     * Prints out a vector
     * @param x vector
     */
    public static void printVec(double[] x) {
        System.out.println("*************");
        for (int i = 0; i < x.length; i++) {
            System.out.println(x[i]);
        }
        System.out.println("*************");
    }

    /**
     * Returns the transpose of a matrix
     * @param a matrix
     * @return transpose of a
     */
    public static double[][] transpose(double[][] a) {
        double[][] at = new double[a[0].length][a.length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                at[j][i] = a[i][j];
            }
        }
        return at;
    }

    /**
     * Returns the L2 norm square of a vector
     * @param x vector
     * @return L2 norm square of vector x
     */
    public static double getL2Norm(double[] x) {
        double res = 0.;
        for (int i = 0; i < x.length; i++) {
            res += (x[i] * x[i]);
        }
        return res;
    }

    /**
     * Returns c = ab, where a is a matrix and b a vector
     * @param a matrix
     * @param b vector
     * @return product ab
     */
    public static double[] multMatVec(double[][] a, double[] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < c.length; i++) {
            c[i] = 0.;
            for (int j = 0; j < a[i].length; j++) {
                c[i] += (a[i][j] * b[j]);
            }
        }
        return c;
    }

    /**
     * Returns the product of two matrices
     * @param a matrix
     * @param b matrix
     * @return c = ab
     */
    public static double[][] multMat(double[][] a, double[][] b) {
        double[][] c = new double[a.length][b[0].length];
        for (int i = 0; i < c.length; i++) {
            for (int j = 0; j < c[0].length; j++) {
                c[i][j] = 0.;
                for (int k = 0; k < a[0].length; k++)
                    c[i][j] += (a[i][k] * b[k][j]);
            }
        }
        return c;
    }

    /**
     * Performs the Cholesky factorization of a symmetric and positive
     * definite matrix A: A = LL^T, where L is lower triangular with positive
     * diagonal elements. This method returns L
     * @param A input matrix (symmetric and positive definite)
     * @return Cholesky factor L of A = LL^T
     */
    public static double[][] cholesky(double[][] A) {
        int m = A.length;
        double[][] l = new double[m][m];
        for (int k = 0; k < m; k++) {
            l[k][k] = A[k][k];
            for (int j = 0; j < k; j++) {
                l[k][k] -= (l[k][j] * l[k][j]);
            }
            l[k][k] = Math.sqrt(l[k][k]);
            for (int i = k + 1; i < m; i++) {
                l[i][k] = A[i][k];
                for (int j = 0; j < k; j++) {
                    l[i][k] -= (l[i][j] * l[k][j]);
                }
                l[i][k] /= l[k][k];
            }
        }
        return l;
    }

    /**
     * Solves the linear system ax = b, where a is lower triangular
     * @param a matrix
     * @param b vector
     * @return solution x to ax = b
     */
    public static double[] forwardElimination(double[][] a, double[] b) {
        int m = b.length;
        double[] x = new double[m];
        for (int i = 0; i < m; i++) {
            x[i] = b[i];
            for (int j = 0; j < i; j++) {
                x[i] -= (a[i][j] * x[j]);
            }
            x[i] /= a[i][i];
        }
        return x;
    }

    /**
     * Solves the system ax = b, where a is upper triangular
     * @param a matrix
     * @param b vector
     * @return solution x to ax = b
     */
    public static double[] backwardSubstitution(double[][] a, double[] b) {
        int m = b.length;
        double[] x = new double[m];
        for (int i = m - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < m; j++) {
                x[i] -= (a[i][j] * x[j]);
            }
            x[i] /= a[i][i];
        }
        return x;
    }

    /**
     * Solves the linear system ax = b, where matrix a is symmetric and
     * positive definite
     * @param a matrix
     * @param b vector
     * @return solution x to ax = b
     */
    public static double[] choleskySolve(double[][] a, double[] b) {
        double[][] l = cholesky(a);
        double[][] lt = transpose(l);
        double[] y = forwardElimination(l, b);
        double[] x = backwardSubstitution(lt, y);
        return x;
    }

    /**
     * Computes the inner product between two vectors
     * @param x vector of dim n
     * @param y vector of dim n
     * @return sum_{i=1}^n x_iy_i
     */
    public static double geInnerProd(double[] x, double[] y) {
        double res = 0.;
        for (int i = 0; i < x.length; i++) {
            res += (x[i] * y[i]);
        }
        return res;
    }
}