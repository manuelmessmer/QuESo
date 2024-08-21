from queso.tests.ggq_tube.test_ggq_tube import TestGGQTube
from queso.tests.ggq_rule_1d.test_ggq_rule_1d import TestGGQ1d
from queso.tests.b_spline_volume.test_b_spline_volume import TestBSplineVolume
from queso.tests.parameter_container.test_parameters import TestParametersContainer
from queso.tests.boundary_conditions.test_boundary_conditions import TestBoundaryConditions

try:
    import KratosMultiphysics as KM
    kratos_available = True
except:
    kratos_available = False

if kratos_available:
    from queso.tests.ggq_cantilever_kratos.test_ggq_cantilever_kratos import TestGGQCantileverKratos
    from queso.tests.trimmed_cantilever_kratos.test_trimmed_cantilever_kratos import TestTrimmedCantileverKratos
    from queso.tests.steering_knuckle_kratos.test_strain_energy_steering_knuckle import TestStrainEnergySteeringKnuckleKratos
    from queso.tests.boundary_conditions_kratos.test_boundary_conditions_kratos import TestBoundaryConditionsKratos
import unittest
import sys

def PyQuESoTestSuite():
    test_suite = unittest.TestSuite()
    if kratos_available:
        test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGGQCantileverKratos))
        test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTrimmedCantileverKratos))
        test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStrainEnergySteeringKnuckleKratos))
        test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestBoundaryConditionsKratos))
    else:
        print("Warning :: Tests with KratosMultiphysics dependencies are skipped.")

    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGGQTube))
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGGQ1d))
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestBSplineVolume))
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestParametersContainer))
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestBoundaryConditions))

    return test_suite


def main():

    test_suite = PyQuESoTestSuite()
    runner = unittest.TextTestRunner()
    result = runner.run(test_suite)

    sys.exit(not result.wasSuccessful())

if __name__ == "__main__":
    main()