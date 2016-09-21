﻿using System; 
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Utilities;
using System.Runtime.Serialization.Formatters.Binary;
using System.IO;
using System.Xml;

namespace UtilitiesUnitTest
{ 
   [TestClass] 
   public class UtilitiesTests
    {
        [TestMethod]
        public void IntegratorTest()
        {
            StateSpaceEOMS dynamics = new StateSpaceEOMS();
            Matrix<double> tspan = new Matrix<double>(new double[1, 2] { { 0, 20 } });
            Matrix<double> y0 = new Matrix<double>(new double[2, 1] { { 1 }, { 0 } });
            Matrix<double> result = Integrator.RK45(dynamics, tspan, y0);

            System.IO.File.WriteAllText("integratorOut.txt", result.ToString());

            Console.ReadLine();
        }

        [TestMethod]
        public void MatrixVertCatTest()
        {
            // Create Test Matricies.
            Matrix<double> A = new Matrix<double>(2, 3, 1);
            Matrix<double> B = new Matrix<double>(2, 3, 2);

            // Create Correct Answer.
            Matrix<double> C = new Matrix<double>(new double[,] { { 1, 1, 1 }, { 1, 1, 1 }, { 2, 2, 2 }, { 2, 2, 2 } });

            // Test Method.
            Matrix<double> D = Matrix<double>.Vertcat(A, B);

            // Verify Result.
            Assert.AreEqual(C, D);
        }

        [TestMethod]
        public void MatrixHorizCatTest()
        {
            // Create Test Matrices.
            Matrix<double> A = new Matrix<double>(new double[,] { { 1, 1, 1 }, { 1, 1, 1 } });
            Matrix<double> B = new Matrix<double>(new double[,] { { 2, 2, 2 }, { 2, 2, 2 } });

            // Create correct answer.
            Matrix<double> C = new Matrix<double>(new double[,] { { 1, 1, 1, 2, 2, 2 }, { 1, 1, 1, 2, 2, 2 } });

            // Test Method.
            Matrix<double> D = Matrix<double>.Horzcat(A, B);

            // Verify Result  .       
            Assert.AreEqual(C, D);
        }

        [TestMethod]
        public void MatrixCumProdTest()
        {
            // Create Test Matrix.
            Matrix<double> A = new Matrix<double>(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });

            // Create Correct Answer.
            Matrix<double> C = new Matrix<double>(new double[,] { { 1, 2, 3 }, { 4, 10, 18 }, { 28, 80, 162 } });

            //Test Method.
            Matrix<double> B = Matrix<double>.Cumprod(A);

            //Verify Result.
            Assert.AreEqual(C, B);
        }

        [TestMethod]
        public void InsertRowTest()
        {
            Matrix<double> A = new Matrix<double>(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
            A.SetRow(2, new Matrix<double>(new double[,] { { 1, 2, 3 } }));
            Matrix<double> B = new Matrix<double>(new double[,] { { 1, 2, 3 }, { 1, 2, 3 }, { 7, 8, 9 } });

            Assert.AreEqual(B, A);
        }
        [TestMethod]
        public void MatrixDeepCopyTest()
        {
            Matrix<double> A = new Matrix<double>(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });

            Matrix<double> B = DeepCopy.Copy(A);

            Assert.AreEqual(B, A);
        }

        [TestMethod]
        public void DeepCopyTest()
        {
            TestCopyClass c1 = new TestCopyClass
            {
                Num = 1,
                Name = "Morgan",
                Helper = new TestCopyHelperClass
                {
                    helperString = "Help",
                    helperInt = 2
                }
            };
            //Test copy
            TestCopyClass c2 = DeepCopy.Copy<TestCopyClass>(c1);
            Assert.AreEqual(c1.Helper.helperInt, c2.Helper.helperInt);
            Assert.AreEqual(c1.Num, c2.Num);
            //Test Modify to property
            c2.Num = 5;
            Assert.AreNotEqual(c1.Num, c2.Num);
            Assert.AreEqual(5, c2.Num);
            Assert.AreEqual(1, c1.Num);
            //Test Modification to Helper class
            c2.Helper.helperInt = 12;
            Assert.AreNotEqual(c1.Helper.helperInt, c2.Helper.helperInt);
            Assert.AreEqual(12, c2.Helper.helperInt);
            c2.Helper.helperString = "poop";
            Assert.AreNotEqual(c1.Helper.helperString, c2.Helper.helperString);
            Assert.AreEqual("poop", c2.Helper.helperString);
            //Test clone constructor
            TestCopyClass c3 = new TestCopyClass(c1);
            Assert.AreEqual(c1.Num, c3.Num);
            Assert.AreEqual(c1.Helper.helperString, c3.Helper.helperString);
            c3.Num = 7;
            c3.Helper.helperInt = 7;
            Assert.AreNotEqual(c1.Helper.helperInt, c3.Helper.helperInt);
            Assert.AreEqual(7, c3.Helper.helperInt);
            Assert.AreNotEqual(c1.Num, c3.Num);
        }
    }
    [Serializable]
    public class TestCopyClass
    {
        public int Num;
        public string Name;
        public TestCopyHelperClass Helper;
        public TestCopyClass() { }
        public TestCopyClass(TestCopyClass other)
        {
            TestCopyClass copy = DeepCopy.Copy<TestCopyClass>(other);
            Num = copy.Num;
            Name = copy.Name;
            Helper = copy.Helper;
        }
    }
    [Serializable]
    public class TestCopyHelperClass
    {
        public string helperString;
        public int helperInt;
    }

} 
  
