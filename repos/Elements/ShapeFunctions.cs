using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Elements
{
    /// <summary>
    /// Shape functions and derivatives for each element type
    /// </summary>
    public class ShapeFunctions
    {
        // Methods
        /// <summary>
        /// Evaluates shape function derivatives dN_dXI
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static Matrix EvaluateShapeFunctionDerivatives(Elements Element, double xi, double eta)
        {
            string type = Element.Type;
            switch (type)
            {
                // 3-Noded Triangles
                case "3NT":
                    Matrix dN_3NT = ShapeFunctionDerivatives_3NT(xi, eta);
                    return dN_3NT;

                // 4-Noded Triangles
                case "4NT":
                    Matrix dN_4NT = ShapeFunctionDerivatives_4NT(xi, eta);
                    return dN_4NT;

                // 6-Noded Triangles
                case "6NT":
                    Matrix dN_6NT = ShapeFunctionDerivatives_6NT(xi, eta);
                    return dN_6NT;

                // 4-Noded Quads
                case "4NQ":
                    Matrix dN_4NQ = ShapeFunctionDerivatives_4NQ(xi, eta);
                    return dN_4NQ;

                // 6-Noded Quads:
                case "6NQ":
                    Matrix dN_6NQ = ShapeFunctionDerivatives_6NQ(xi, eta);
                    return dN_6NQ;

                // 8-Noded Quads
                case "8NQ":
                    Matrix dN_8NQ = ShapeFunctionDerivatives_8NQ(xi, eta);
                    return dN_8NQ;

                // Default
                default:
                    throw new Exception("Invalid element type when calculating derivatives of shape functions");
            }
        }

        /// <summary>
        /// Shape Functions: 3-Noded Triangle
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static double[] ShapeFunctions_3NT(double xi, double eta)
        {
            double[] N = new double[3];

            double N1 = 1.0 - xi - eta;
            double N2 = xi;
            double N3 = eta;

            N[0] = N1; N[1] = N2; N[2] = N3;

            return N;
        }

        /// <summary>
        /// Derivatives of Shape Functions: 3-Noded Triangle
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static Matrix ShapeFunctionDerivatives_3NT(double xi, double eta)
        {
            Matrix dN = new(4, 6);

            double dN1_dxi = -1.0;
            double dN1_deta = -1.0;

            double dN2_dxi = 1.0;
            double dN2_deta = 0.0;

            double dN3_dxi = 0.0;
            double dN3_deta = 1.0;

            dN[0, 0] = dN1_dxi;  dN[0, 1] = 0.0; dN[0, 2] = dN2_dxi;  dN[0, 3] = 0.0; dN[0, 4] = dN3_dxi;  dN[0, 5] = 0.0;
            dN[1, 0] = dN1_deta; dN[1, 1] = 0.0; dN[1, 2] = dN2_deta; dN[1, 3] = 0.0; dN[1, 4] = dN3_deta; dN[1, 5] = 0.0;
            dN[2, 0] = 0.0; dN[2, 1] = dN1_dxi;  dN[2, 2] = 0.0; dN[2, 3] = dN2_dxi;  dN[2, 4] = 0.0; dN[2, 5] = dN3_dxi;
            dN[3, 0] = 0.0; dN[3, 1] = dN1_deta; dN[3, 2] = 0.0; dN[3, 3] = dN2_deta; dN[3, 4] = 0.0; dN[3, 5] = dN3_deta;

            return dN;
        }

        /// <summary>
        /// Shape Functions: 4-Noded Triangle
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static double[] ShapeFunctions_4NT(double xi, double eta)
        {
            double[] N = new double[4];

            double N1 = eta;
            double N2 = (1.0 - 2.0 * xi - eta) * (1.0 - xi - eta) / (1.0 - eta);
            double N3 = 4.0 * xi * (1.0 - xi - eta) / (1.0 - eta);
            double N4 = xi * (2.0 * xi + eta - 1.0) / (1.0 - eta);

            N[0] = N1; N[1] = N2; N[2] = N3; N[3] = N4;

            return N;
        }

        /// <summary>
        /// Derivatives of Shape Functions: 4-Noded Triangle
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static Matrix ShapeFunctionDerivatives_4NT(double xi, double eta)
        {
            Matrix dN = new(4, 8);

            double dN1_dxi = 0.0;
            double dN1_deta = 1.0;

            double dN2_dxi = -(eta + 2.0 * xi - 1.0) / (eta - 1.0) - (2.0 * (eta + xi - 1.0)) / (eta - 1.0);
            double dN2_deta = (eta + xi - 1) * (eta + 2 * xi - 1) / ((eta - 1) * (eta - 1.0)) - (eta + xi - 1) / (eta - 1) - (eta + 2 * xi - 1) / (eta - 1);

            double dN3_dxi = 4.0 * (eta + xi - 1.0) / (eta - 1.0) + 4.0 * xi / (eta - 1.0);
            double dN3_deta = 4.0 * xi / (eta - 1.0) - 4.0 * xi * (eta + xi - 1.0) / ((eta - 1.0) * (eta - 1.0));

            double dN4_dxi = -(eta + 2.0 * xi - 1.0) / (eta - 1.0) - (2.0 * xi) / (eta - 1.0);
            double dN4_deta = xi * (eta + 2.0 * xi - 1.0) / ((eta - 1.0) * (eta - 1.0)) - xi / (eta - 1.0);

            dN[0, 0] = dN1_dxi;  dN[0, 1] = 0.0; dN[0, 2] = dN2_dxi;  dN[0, 3] = 0.0; dN[0, 4] = dN3_dxi;  dN[0, 5] = 0.0; dN[0, 6] = dN4_dxi;  dN[0, 7] = 0.0;
            dN[1, 0] = dN1_deta; dN[1, 1] = 0.0; dN[1, 2] = dN2_deta; dN[1, 3] = 0.0; dN[1, 4] = dN3_deta; dN[1, 5] = 0.0; dN[1, 6] = dN4_deta; dN[1, 7] = 0.0;
            dN[2, 0] = 0.0; dN[2, 1] = dN1_dxi;  dN[2, 2] = 0.0; dN[2, 3] = dN2_dxi;  dN[2, 4] = 0.0; dN[2, 5] = dN3_dxi;  dN[2, 6] = 0.0; dN[2, 7] = dN4_dxi;
            dN[3, 0] = 0.0; dN[3, 1] = dN1_deta; dN[3, 2] = 0.0; dN[3, 3] = dN2_deta; dN[3, 4] = 0.0; dN[3, 5] = dN3_deta; dN[3, 6] = 0.0; dN[3, 7] = dN4_deta;

            return dN;
        }

        /// <summary>
        /// Shape Functions: 6-Noded Triangle
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static double[] ShapeFunctions_6NT(double xi, double eta)
        {
            double[] N = new double[6];

            double N1 = eta * (2.0 * eta - 1.0);
            double N2 = 4.0 * eta * (1.0 - xi - eta);
            double N3 = (1.0 - xi - eta) * (1.0 - 2.0 * xi - 2.0 * eta);
            double N4 = 4.0 * xi * (1.0 - xi - eta);
            double N5 = xi * (2.0 * xi - 1.0);
            double N6 = 4.0 * xi * eta;

            N[0] = N1; N[1] = N2; N[2] = N3;
            N[3] = N4; N[4] = N5; N[5] = N6;

            return N;
        }
        /// <summary>
        /// Derivatives of Shape Functions: 6-Noded Triangle
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static Matrix ShapeFunctionDerivatives_6NT(double xi, double eta)
        {
            Matrix dN = new(4, 12);

            double dN1_dxi = 0;
            double dN1_deta = 4 * eta - 1;

            double dN2_dxi = -4 * eta;
            double dN2_deta = 4 - 4 * xi - 8 * eta;

            double dN3_dxi = 4 * eta + 4 * xi - 3;
            double dN3_deta = 4 * eta + 4 * xi - 3;

            double dN4_dxi = 4 - 8 * xi - 4 * eta;
            double dN4_deta = -4 * xi;

            double dN5_dxi = 4 * xi - 1;
            double dN5_deta = 0;

            double dN6_dxi = 4 * eta;
            double dN6_deta = 4 * xi;

            dN[0, 0] = dN1_dxi;  dN[0, 1] = 0.0; dN[0, 2] = dN2_dxi; dN[0, 3] = 0.0;  dN[0, 4] = dN3_dxi; dN[0, 5] = 0.0;  dN[0, 6] = dN4_dxi;  dN[0, 7] = 0.0; dN[0, 8] = dN5_dxi;  dN[0, 9] = 0.0; dN[0, 10] = dN6_dxi;  dN[0, 11] = 0.0;
            dN[1, 0] = dN1_deta; dN[1, 1] = 0.0; dN[1, 2] = dN2_deta; dN[1, 3] = 0.0; dN[1, 4] = dN3_deta; dN[1, 5] = 0.0; dN[1, 6] = dN4_deta; dN[1, 7] = 0.0; dN[1, 8] = dN5_deta; dN[1, 9] = 0.0; dN[1, 10] = dN6_deta; dN[1, 11] = 0.0;
            dN[2, 0] = 0.0; dN[2, 1] = dN1_dxi;  dN[2, 2] = 0.0; dN[2, 3] = dN2_dxi;  dN[2, 4] = 0.0; dN[2, 5] = dN3_dxi;  dN[2, 6] = 0.0; dN[2, 7] = dN4_dxi;  dN[2, 8] = 0.0; dN[2, 9] = dN5_dxi;  dN[2, 10] = 0.0; dN[2, 11] = dN6_dxi;
            dN[3, 0] = 0.0; dN[3, 1] = dN1_deta; dN[3, 2] = 0.0; dN[3, 3] = dN2_deta; dN[3, 4] = 0.0; dN[3, 5] = dN3_deta; dN[3, 6] = 0.0; dN[3, 7] = dN4_deta; dN[3, 8] = 0.0; dN[3, 9] = dN5_deta; dN[3, 10] = 0.0; dN[3, 11] = dN6_deta;

            return dN;
        }

        /// <summary>
        /// Shape Functions: 4-Noded Quad
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static double[] ShapeFunctions_4NQ(double xi, double eta)
        {
            double[] N = new double[4];

            double N1 = 1.0 / 4.0 * (xi - 1.0) * (eta - 1.0);
            double N2 = -1.0 / 4.0 * (xi + 1.0) * (eta - 1.0);
            double N3 = 1.0 / 4.0 * (xi + 1.0) * (eta + 1.0);
            double N4 = -1.0 / 4.0 * (xi - 1.0) * (eta + 1.0);

            N[0] = N1; N[1] = N2; N[2] = N3; N[3] = N4;

            return N;
        }
        /// <summary>
        /// Derivatives of Shape Functions: 4-Noded Quad
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static Matrix ShapeFunctionDerivatives_4NQ(double xi, double eta)
        {
            Matrix dN = new(4, 8);

            double dN1_dxi  = 1.0 / 4.0 * (eta - 1.0);
            double dN1_deta = 1.0 / 4.0 * (xi - 1.0);

            double dN2_dxi  = -1.0 / 4.0 * (eta - 1.0);
            double dN2_deta = -1.0 / 4.0 * (xi + 1.0);

            double dN3_dxi  = 1.0 / 4.0 * (eta + 1.0);
            double dN3_deta = 1.0 / 4.0 * (xi + 1.0);

            double dN4_dxi  = -1.0 / 4.0 * (eta + 1.0);
            double dN4_deta = -1.0 / 4.0 * (xi - 1.0);

            dN[0,0] = dN1_dxi;  dN[0,1] = 0.0;      dN[0,2] = dN2_dxi;  dN[0,3] = 0.0;      dN[0,4] = dN3_dxi;  dN[0,5] = 0.0;      dN[0,6] = dN4_dxi;  dN[0,7] = 0.0;
            dN[1,0] = dN1_deta; dN[1,1] = 0.0;      dN[1,2] = dN2_deta; dN[1,3] = 0.0;      dN[1,4] = dN3_deta; dN[1,5] = 0.0;      dN[1,6] = dN4_deta; dN[1,7] = 0.0;
            dN[2,0] = 0.0;      dN[2,1] = dN1_dxi;  dN[2,2] = 0.0;      dN[2,3] = dN2_dxi;  dN[2,4] = 0.0;      dN[2,5] = dN3_dxi;  dN[2,6] = 0.0;      dN[2,7] = dN4_dxi;
            dN[3,0] = 0.0;      dN[3,1] = dN1_deta; dN[3,2] = 0.0;      dN[3,3] = dN2_deta; dN[3,4] = 0.0;      dN[3,5] = dN3_deta; dN[3,6] = 0.0;      dN[3,7] = dN4_deta;

            return dN;
        }

        /// <summary>
        /// Shape Functions: 6-Noded Quad
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static double[] ShapeFunctions_6NQ(double xi, double eta)
        {
            double[] N = new double[6];

            double N1 = -1.0 / 4.0 * xi * (xi - 1.0) * (eta - 1.0);
            double N2 = 1.0 / 2.0 * (xi + 1.0) * (xi - 1.0) * (eta - 1.0);
            double N3 = -1.0 / 4.0 * xi * (xi + 1.0) * (eta - 1.0);
            double N4 = 1.0 / 4.0 * xi * (xi + 1.0) * (eta + 1.0);
            double N5 = -1.0 / 2.0 * (xi + 1.0) * (xi - 1.0) * (eta + 1.0);
            double N6 = 1.0 / 4.0 * xi * (xi - 1.0) * (eta + 1.0);

            N[0] = N1; N[1] = N2; N[2] = N3; 
            N[3] = N4; N[4] = N5; N[5] = N6;

            return N;
        }
        /// <summary>
        /// Derivatives of Shape Functions: 6-Noded Quad
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static Matrix ShapeFunctionDerivatives_6NQ(double xi, double eta)
        {
            Matrix dN = new(4, 12);

            double dN1_dxi = 1.0 / 4.0 * (-2.0 * eta * xi + eta + 2.0 * xi - 1.0);
            double dN1_deta = -1.0 / 4.0 * (xi - 1.0) * xi;

            double dN2_dxi = xi * (eta - 1.0);
            double dN2_deta = 1.0 / 2.0 * (xi - 1.0) * (xi + 1.0);

            double dN3_dxi = -1.0 / 4.0 * (eta - 1.0) * (2.0 * xi + 1.0);
            double dN3_deta = -1.0 / 4.0 * (xi) * (xi + 1.0);

            double dN4_dxi = 1.0 / 4.0 * (eta + 1.0) * (2.0 * xi + 1.0);
            double dN4_deta = 1.0 / 4.0 * xi * (xi + 1.0);

            double dN5_dxi = -xi * (eta + 1.0);
            double dN5_deta = -1.0 / 2.0 * (xi - 1.0) * (xi + 1.0);

            double dN6_dxi = 1.0 / 4.0 * (eta + 1.0) * (2.0 * xi - 1.0);
            double dN6_deta = 1.0 / 4.0 * xi * (xi - 1.0);

            dN[0, 0] = dN1_dxi;  dN[0, 1] = 0.0; dN[0, 2] = dN2_dxi;  dN[0, 3] = 0.0; dN[0, 4] = dN3_dxi;  dN[0, 5] = 0.0; dN[0, 6] = dN4_dxi;  dN[0, 7] = 0.0; dN[0, 8] = dN5_dxi;  dN[0, 9] = 0.0; dN[0, 10] = dN6_dxi;  dN[0, 11] = 0.0;
            dN[1, 0] = dN1_deta; dN[1, 1] = 0.0; dN[1, 2] = dN2_deta; dN[1, 3] = 0.0; dN[1, 4] = dN3_deta; dN[1, 5] = 0.0; dN[1, 6] = dN4_deta; dN[1, 7] = 0.0; dN[1, 8] = dN5_deta; dN[1, 9] = 0.0; dN[1, 10] = dN6_deta; dN[1, 11] = 0.0;
            dN[2, 0] = 0.0; dN[2, 1] = dN1_dxi; dN[2, 2] = 0.0;  dN[2, 3] = dN2_dxi;  dN[2, 4] = 0.0; dN[2, 5] = dN3_dxi;  dN[2, 6] = 0.0; dN[2, 7] = dN4_dxi;  dN[2, 8] = 0.0; dN[2, 9] = dN5_dxi;  dN[2, 10] = 0.0; dN[2, 11] = dN6_dxi;
            dN[3, 0] = 0.0; dN[3, 1] = dN1_deta; dN[3, 2] = 0.0; dN[3, 3] = dN2_deta; dN[3, 4] = 0.0; dN[3, 5] = dN3_deta; dN[3, 6] = 0.0; dN[3, 7] = dN4_deta; dN[3, 8] = 0.0; dN[3, 9] = dN5_deta; dN[3, 10] = 0.0; dN[3, 11] = dN6_deta;

            return dN;
        }

        /// <summary>
        /// Shape Functions: 8-Noded Quad
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static double[] ShapeFunctions_8NQ(double xi, double eta)
        {
            double[] N = new double[8];

            double N1 = 1.0 / 4.0 * (-xi - eta - 1.0) * (xi - 1.0) * (eta - 1.0);
            double N2 = 1.0 / 2.0 * (xi + 1.0) * (xi - 1.0) * (eta - 1.0);
            double N3 = -1.0 / 4.0 * (xi - eta - 1.0) * (xi + 1.0) * (eta - 1.0);
            double N4 = -1.0 / 2.0 * (eta + 1.0) * (eta - 1.0) * (xi + 1.0);
            double N5 = -1.0 / 4.0 * (1.0 - xi - eta) * (xi + 1.0) * (eta + 1.0);
            double N6 = -1.0 / 2.0 * (xi + 1.0) * (xi - 1.0) * (eta + 1.0);
            double N7 = 1.0 / 4.0 * (xi + 1.0 - eta) * (xi - 1.0) * (eta + 1.0);
            double N8 = 1.0 / 2.0 * (eta + 1.0) * (eta - 1.0) * (xi - 1.0);

            N[0] = N1; N[1] = N2; N[2] = N3; N[3] = N4;
            N[4] = N5; N[5] = N6; N[6] = N7; N[7] = N8;

            return N;
        }
        /// <summary>
        /// Derivatives of Shape Functions: 4-Noded Quad
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        private static Matrix ShapeFunctionDerivatives_8NQ(double xi, double eta)
        {
            Matrix dN = new(4, 16);

            double dN1_dxi = -1.0 / 4.0 * (eta - 1.0) * (eta + 2.0 * xi);
            double dN1_deta = -1.0 / 4.0 * (xi - 1.0) * (2.0 * eta + xi);

            double dN2_dxi = (eta - 1.0) * xi;
            double dN2_deta = 1.0 / 2.0 * (xi - 1.0) * (xi + 1.0);

            double dN3_dxi = 1.0 / 4.0 * (eta - 1.0) * (eta - 2.0 * xi);
            double dN3_deta = 1.0 / 4.0 * (xi + 1.0) * (2.0 * eta - xi);

            double dN4_dxi = -1.0 / 2.0 * (eta - 1.0) * (eta + 1.0);
            double dN4_deta = -1.0 * eta * (xi + 1.0);

            double dN5_dxi = 1.0 / 4.0 * (eta + 1.0) * (eta + 2.0 * xi);
            double dN5_deta = 1.0 / 4.0 * (xi + 1.0) * (2.0 * eta + xi);

            double dN6_dxi = -1.0 * (eta + 1.0) * xi;
            double dN6_deta = -1.0 / 2.0 * (xi - 1.0) * (xi + 1.0);

            double dN7_dxi = -1.0 / 4.0 * (eta + 1.0) * (eta - 2.0 * xi);
            double dN7_deta = -1.0 / 4.0 * (xi - 1.0) * (2.0 * eta - xi);

            double dN8_dxi = 1.0 / 2.0 * (eta - 1.0) * (eta + 1.0);
            double dN8_deta = eta * (xi - 1.0);

            dN[0, 0] = dN1_dxi;  dN[0, 1] = 0.0; dN[0, 2] = dN2_dxi; dN[0, 3] = 0.0;  dN[0, 4] = dN3_dxi;  dN[0, 5] = 0.0; dN[0, 6] = dN4_dxi;  dN[0, 7] = 0.0; dN[0, 8] = dN5_dxi;  dN[0, 9] = 0.0; dN[0, 10] = dN6_dxi;  dN[0, 11] = 0.0; dN[0, 12] = dN7_dxi;  dN[0, 13] = 0.0; dN[0, 14] = dN8_dxi;  dN[0, 15] = 0.0;
            dN[1, 0] = dN1_deta; dN[1, 1] = 0.0; dN[1, 2] = dN2_deta; dN[1, 3] = 0.0; dN[1, 4] = dN3_deta; dN[1, 5] = 0.0; dN[1, 6] = dN4_deta; dN[1, 7] = 0.0; dN[1, 8] = dN5_deta; dN[1, 9] = 0.0; dN[1, 10] = dN6_deta; dN[1, 11] = 0.0; dN[1, 12] = dN7_deta; dN[1, 13] = 0.0; dN[1, 14] = dN8_deta; dN[1, 15] = 0.0;
            dN[2, 0] = 0.0; dN[2, 1] = dN1_dxi; dN[2, 2] = 0.0;  dN[2, 3] = dN2_dxi; dN[2, 4] = 0.0;  dN[2, 5] = dN3_dxi;  dN[2, 6] = 0.0; dN[2, 7] = dN4_dxi;  dN[2, 8] = 0.0; dN[2, 9] = dN5_dxi;  dN[2, 10] = 0.0; dN[2, 11] = dN6_dxi;  dN[2, 12] = 0.0; dN[2, 13] = dN7_dxi;  dN[2, 14] = 0.0; dN[2, 15] = dN8_dxi;
            dN[3, 0] = 0.0; dN[3, 1] = dN1_deta; dN[3, 2] = 0.0; dN[3, 3] = dN2_deta; dN[3, 4] = 0.0; dN[3, 5] = dN3_deta; dN[3, 6] = 0.0; dN[3, 7] = dN4_deta; dN[3, 8] = 0.0; dN[3, 9] = dN5_deta; dN[3, 10] = 0.0; dN[3, 11] = dN6_deta; dN[3, 12] = 0.0; dN[3, 13] = dN7_deta; dN[3, 14] = 0.0; dN[3, 15] = dN8_deta;

            return dN;
        }
    }
}
