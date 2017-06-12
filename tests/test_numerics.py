"""Methods for testing the subroutines in the numerics module."""
import unittest as ut

class TestFactorial(ut.TestCase):
    """ Tests of the factorial subroutine."""

    def test_1(self):
        from phenum.numerics import factorial
        num = 1
        out = 1
        self.assertEqual(factorial(num),out)

    def test_2(self):
        from phenum.numerics import factorial
        num = 10
        out = 3628800
        self.assertEqual(factorial(num),out)


    def test_3(self):
        from phenum.numerics import factorial
        num = 9
        out = 362880
        self.assertEqual(factorial(num),out)


    def test_4(self):
        from phenum.numerics import factorial
        num = 42
        out = 1405006117752879898543142606244511569936384000000000
        self.assertEqual(factorial(num),out)


    def test_5(self):
        from phenum.numerics import factorial
        num = 15
        out = 1307674368000
        self.assertEqual(factorial(num),out)


    def test_6(self):
        from phenum.numerics import factorial
        num = 27
        out = 10888869450418352160768000000
        self.assertEqual(factorial(num),out)


    def test_7(self):
        from phenum.numerics import factorial
        num = 33
        out = 8683317618811886495518194401280000000
        self.assertEqual(factorial(num),out)


    def test_8(self):
        from phenum.numerics import factorial
        num = 7
        out = 5040
        self.assertEqual(factorial(num),out)


    def test_9(self):
        from phenum.numerics import factorial
        num = 2
        out = 2
        self.assertEqual(factorial(num),out)


    def test_10(self):
        from phenum.numerics import factorial
        num = 0
        out = 1
        self.assertEqual(factorial(num),out)

class TestBinomialCoefficient(ut.TestCase):
    """ Tests of the binomial_coefficient subroutine."""

    def test_1(self):
        from phenum.numerics import binomial_coefficient
        n = 47
        r = 11
        out = 17417133617
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_2(self):
        from phenum.numerics import binomial_coefficient
        n = 48
        r = 29
        out = 11541847896480
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_3(self):
        from phenum.numerics import binomial_coefficient
        n = 24
        r = 1
        out = 24
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_4(self):
        from phenum.numerics import binomial_coefficient
        n = 21
        r = 2
        out = 210
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_5(self):
        from phenum.numerics import binomial_coefficient
        n = 35
        r = 27
        out = 23535820
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_6(self):
        from phenum.numerics import binomial_coefficient
        n = 44
        r = 6
        out = 7059052
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_7(self):
        from phenum.numerics import binomial_coefficient
        n = 47
        r = 25
        out = 14833897694226
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_8(self):
        from phenum.numerics import binomial_coefficient
        n = 42
        r = 13
        out = 25518731280
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_9(self):
        from phenum.numerics import binomial_coefficient
        n = 32
        r = 22
        out = 64512240
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_10(self):
        from phenum.numerics import binomial_coefficient
        n = 1
        r = 1
        out = 1 
        self.assertEqual(binomial_coefficient(n,r),out)

    def test_11(self):
        from phenum.numerics import binomial_coefficient
        n = 0
        r = 0
        out = 1 
        self.assertEqual(binomial_coefficient(n,r),out)
        
    def test_12(self):
        from phenum.numerics import binomial_coefficient
        n = 10
        r = 0
        out = 1 
        self.assertEqual(binomial_coefficient(n,r),out)
        
    def test_13(self):
        from phenum.numerics import binomial_coefficient
        n = 1
        r = 10
        out = 0 
        self.assertEqual(binomial_coefficient(n,r),out)
