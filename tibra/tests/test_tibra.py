import unittest

from tibra.tests.ggq_cantilever_kratos.test_ggq_cantilever_kratos import TestGGQCantileverKratos
from tibra.tests.trimmed_cantilever_kratos.test_trimmed_cantilever_kratos import TestTrimmedCantileverKratos

try:
    import KratosMultiphysics as KM
    kratos_available = True
except:
    print("KratosMultiphysics is not available")
    kratos_available = False

def PyTIBRATestSuite():
    test_suite = unittest.TestSuite()
    if kratos_available:
        test_suite.addTest(unittest.makeSuite(TestGGQCantileverKratos))
        test_suite.addTest(unittest.makeSuite(TestTrimmedCantileverKratos))
    else:
        print("Warning :: Tests with KratosMultiphysics dependency are skipped.")

    return test_suite


def main():

    test_suite = PyTIBRATestSuite()

    runner = unittest.TextTestRunner()
    runner.run(test_suite)

if __name__ == "__main__":
    main()