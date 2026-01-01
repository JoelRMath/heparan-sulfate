package heparansulfate;

/**
 * Provides static methods for common matrix and vector operations.
 * <p>
 * This class includes utilities for printing, matrix multiplication, transposition, 
 * and solving linear systems using Cholesky factorization.
 */
public class MatrixOp {

    /**
     * Default constructor.
     */
    public MatrixOp() {
    }

    /**
     * Prints a matrix to the standard output.
     * @param mat The matrix to print.
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
     * Prints a symbolic representation of a matrix.
     * <p>
     * Entries with absolute values less than or equal to {@code eps} are printed as ".", 
     * while other entries are printed as "x".
     * @param mat The matrix to print.
     * @param eps The threshold value below which an entry is considered zero.
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
     * Returns a {@code String} representation of a matrix.
     * @param mat The matrix to convert.
     * @return A tab-delimited string representation of {@code mat}.
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
     * Returns a {@code String} representation of a vector.
     * @param x The vector to convert.
     * @return A string representation of vector {@code x}.
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
     * Prints a vector to the standard output.
     * @param x The vector to print.
     */
    public static void printVec(double[] x) {
        System.out.println("*************");
        for (int i = 0; i < x.length; i++) {
            System.out.println(x[i]);
        }
        System.out.println("*************");
    }

    /**
     * Returns the transpose of the given matrix.
     * @param a The input matrix.
     * @return The transposed matrix {@code a^T}.
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
     * Computes the squared L2 norm of a vector.
     * @param x The input vector.
     * @return The value of {@code x^T x}.
     */
    public static double getL2Norm(double[] x) {
        double res = 0.;
        for (int i = 0; i < x.length; i++) {
            res += (x[i] * x[i]);
        }
        return res;
    }

    /**
     * Performs matrix-vector multiplication.
     * @param a The input matrix.
     * @param b The input vector.
     * @return The product vector {@code c = ab}.
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
     * Performs matrix-matrix multiplication.
     * @param a The first matrix.
     * @param b The second matrix.
     * @return The product matrix {@code c = ab}.
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
     * Performs the Cholesky factorization of a symmetric and positive-definite matrix {@code A}.
     * <p>
     * Decomposes {@code A} into {@code A = LL^T}, where {@code L} is a lower triangular matrix 
     * with positive diagonal elements.
     * 
     * @param A The input matrix (must be symmetric and positive definite).
     * @return The lower triangular Cholesky factor {@code L}.
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
     * Solves a linear system {@code ax = b} where the matrix {@code a} is lower triangular 
     * using forward elimination.
     * @param a A lower triangular matrix.
     * @param b The constant vector.
     * @return The solution vector {@code x}.
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
     * Solves a linear system {@code ax = b} where the matrix {@code a} is upper triangular 
     * using backward substitution.
     * @param a An upper triangular matrix.
     * @param b The constant vector.
     * @return The solution vector {@code x}.
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
     * Solves the linear system {@code ax = b} for a symmetric and positive-definite matrix {@code a}.
     * <p>
     * This method uses Cholesky factorization followed by forward elimination 
     * and backward substitution.
     * @param a Symmetric and positive-definite matrix.
     * @param b The constant vector.
     * @return The solution vector {@code x}.
     */
    public static double[] choleskySolve(double[][] a, double[] b) {
        double[][] l = cholesky(a);
        double[][] lt = transpose(l);
        double[] y = forwardElimination(l, b);
        double[] x = backwardSubstitution(lt, y);
        return x;
    }

    /**
     * Computes the inner product (dot product) between two vectors.
     * @param x First vector of dimension {@code n}.
     * @param y Second vector of dimension {@code n}.
     * @return The scalar value {@code sum_{i=1}^n x_i * y_i}.
     */
    public static double geInnerProd(double[] x, double[] y) {
        double res = 0.;
        for (int i = 0; i < x.length; i++) {
            res += (x[i] * y[i]);
        }
        return res;
    }
}