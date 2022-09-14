import unittest

from generalized_Gaussian_quadrature.test_generalized_Gaussian_quadrature import TestGeneralizedGaussianQuadrature
from trimmed_cantilever.test_trimmed_cantilever import TestTrimmedCantilever

def PyTIBRATestSuite():
    test_suite = unittest.TestSuite()
    #test_suite.addTest(unittest.makeSuite(TestGeneralizedGaussianQuadrature))
    test_suite.addTest(unittest.makeSuite(TestTrimmedCantilever))

    return test_suite


def main():

    test_suite = PyTIBRATestSuite()

    runner = unittest.TextTestRunner()
    runner.run(test_suite)

if __name__ == "__main__":
    main()