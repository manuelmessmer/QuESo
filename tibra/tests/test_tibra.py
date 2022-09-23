from tibra.tests.ggq_tube.test_ggq_tube import TestGGQTube
from tibra.tests.ggq_rule_1d.test_ggq_rule_1d import TestGGQ1d
try:
    import KratosMultiphysics as KM
    kratos_available = True
except:
    kratos_available = False

if kratos_available:
    from tibra.tests.ggq_cantilever_kratos.test_ggq_cantilever_kratos import TestGGQCantileverKratos
    from tibra.tests.trimmed_cantilever_kratos.test_trimmed_cantilever_kratos import TestTrimmedCantileverKratos

import unittest
import sys

def PyTIBRATestSuite():
    test_suite = unittest.TestSuite()
    if kratos_available:
        test_suite.addTest(unittest.makeSuite(TestGGQCantileverKratos))
        test_suite.addTest(unittest.makeSuite(TestTrimmedCantileverKratos))
    else:
        print("Warning :: Tests with KratosMultiphysics dependencies are skipped.")

    test_suite.addTest(unittest.makeSuite(TestGGQTube))
    test_suite.addTest(unittest.makeSuite(TestGGQ1d))

    return test_suite


def main():

    test_suite = PyTIBRATestSuite()
    runner = unittest.TextTestRunner()
    result = runner.run(test_suite)

    sys.exit(not result.wasSuccessful())

if __name__ == "__main__":
    main()