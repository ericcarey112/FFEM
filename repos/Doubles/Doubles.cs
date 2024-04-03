public class Doubles
{
    int rows { get; set; }
    public Doubles()
    {
        this.rows = 0;
    }

    // Methods

    /// <summary>
    /// Add []A + [] B
    /// </summary>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <returns></returns>
    public static double[] AddDoubles(double[] A, double[] B)
    {
        int rowsA = A.Length;
        int rowsB = B.Length;
        if (rowsA != rowsB)
        {
            throw new Exception("Vectors must be same length for addition");
        }
        double[] C = new double[rowsA];
        for (int i = 0; i < rowsA; i++)
        {
            C[i] = A[i] + B[i];
        }
        return C;
    }

    /// <summary>
    /// Subtract []A - []B
    /// </summary>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <returns></returns>
    public static double[] SubtractDoubles(double[] A, double[] B)
    {
        int rows = A.Length;
        double[] C = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            C[i] = A[i] - B[i];
        }
        return C;
    }

    /// <summary>
    /// Multiply []A * [] B
    /// </summary>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <returns></returns>
    public static double MultiplyDoubles(double[] A, double[] B)
    {
        int rowsA = A.Length;
        int rowsB = B.Length;
        if (rowsA != rowsB)
        {
            throw new Exception("Vectors must be same length for multiplication");
        }

        double C = 0.0;
        for (int i = 0; i < rowsA; i++)
        {
            C = C + A[i] + B[i];
        }
        return C;
    }

    /// <summary>
    /// Find the max of a double[]
    /// </summary>
    /// <param name="A"></param>
    /// <returns></returns>
    public static double Max(double[] A)
    {
        double max = A[0];
        for (int i = 0; i < A.Length;i++)
        {
            if (A[i] > max) { max = A[i]; }
        }
        return max;
    }

    /// <summary>
    /// Find the min of a double[]
    /// </summary>
    /// <param name="A"></param>
    /// <returns></returns>
    public static double Min(double[] A)
    {
        double min = A[0];
        for (int i = 0; i < A.Length; i++)
        {
            if (A[i] < min) { min = A[i]; }
        }
        return min;
    }

    /// <summary>
    /// Find the index of the minimum value from a double
    /// </summary>
    /// <param name="A"></param>
    /// <returns></returns>
    public static int MinIndex(double[] A)
    {
        double min = A[0];
        int index = 0;
        for(int i = 0; i < A.Length;i++)
        {
            if (A[i] < min)
            {
                min = A[i];
                index = i;
            }
        }
        return index;
    }

    /// <summary>
    /// Find the index of the maximum value from a double
    /// </summary>
    /// <param name="A"></param>
    /// <returns></returns>
    public static int MaxIndex(double[] A)
    {
        double max = A[0];
        int index = 0;
        for (int i = 0; i < A.Length; i++)
        {
            if (A[i] > max)
            {
                max = A[i];
                index = i;
            }
        }
        return index;
    }

    /// <summary>
    /// Multiply a double[] by a constant
    /// </summary>
    /// <param name="A"></param>
    /// <param name="n"></param>
    /// <returns></returns>
    public static double[] MultiplyVectorByConstant(double[] A, double n)
    {
        int rows = A.Length;
        for (int i = 0; i < rows; i++)
        {
            A[i] *= n;
        }
        return A;
    }

    /// <summary>
    /// Create a copy of an existing double[]
    /// </summary>
    /// <param name="A"></param>
    /// <returns></returns>
    public static double[] CreateCopyOfDouble(double[] A)
    {
        int rows = A.Length;
        double[] B = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            B[i] = A[i];
        }
        return B;
    }

    /// <summary>
    /// Display double[] to command line
    /// </summary>
    /// <param name="A"></param>
    public static void DisplayDouble(double[] A)
    {
        int rows = A.GetLength(0);
        for (int i = 0; i < rows; i++)
        {
            Console.WriteLine(A[i]);
        }
    }

    /// <summary>
    /// Return the specified row k of a 2D double[,] A
    /// </summary>
    /// <param name="A"></param>
    /// <param name="k"></param>
    /// <returns></returns>
    public static double[] GetRowFromDouble(double[,] A, int k)
    {
        int cols = A.GetLength(1);
        double[] B = new double[cols];
        for (int i = 0; i < cols; i++)
        {
            B[i] = A[k, i];
        }
        return B;
    }

    /// <summary>
    /// Return the specified col k of a 2D double[,] A
    /// </summary>
    /// <param name="A"></param>
    /// <param name="k"></param>
    /// <returns></returns>
    public static double[] GetColFromDouble(double[,] A, int k)
    {
        int rows = A.GetLength(0);
        double[] B = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            B[i] = A[i, i];
        }
        return B;
    }

    /// <summary>
    /// Makes a list of zeros [rows, cols] of doubles []
    /// </summary>
    /// <param name="rows"></param>
    /// <param name="cols"></param>
    /// <returns></returns>
    public static List<double[]> ListOfZeros(int rows, int cols)
    {
        List<double[]> A = new List<double[]>();
        for (int i = 0; i < rows; i++)
        {
            double[] B = new double[cols];
            for (int j = 0; j < cols; j++)
            {
                B[j] = 0.0;
            }
            A.Add(B);
        }
        return A;
    }

    /// <summary>
    /// Calculates the euclidean norm of a vector
    /// </summary>
    /// <returns></returns>
    public static double EuclideanNorm(double[] A)
    {
        double norm = 0.0;
        for (int i = 0; i < A.Length; i++)
        {
            norm += Math.Pow(Math.Abs(A[i]), 2);
        }
        norm = Math.Sqrt(norm);
        return norm;
    }
}