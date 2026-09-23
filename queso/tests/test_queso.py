import os
import sys
import unittest
from importlib.util import find_spec
from pathlib import Path

from b_spline_volume.test_b_spline_volume import TestBSplineVolume
from boundary_conditions.test_boundary_conditions import TestBoundaryConditions
from ggq_rule_1d.test_ggq_rule_1d import TestGGQ1d
from ggq_tube.test_ggq_tube import TestGGQTube
from multi_component.test_multi_component import (
    TestKratosAnalysis,
    TestMultiComponentModel,
    TestSettingsParser,
)
from settings_container.test_settings import TestSettingsContainer

kratos_available = find_spec("KratosMultiphysics") is not None

if kratos_available:
    from boundary_conditions_kratos.test_boundary_conditions_kratos import (
        TestBoundaryConditionsKratos,
    )
    from ggq_cantilever_kratos.test_ggq_cantilever_kratos import TestGGQCantileverKratos
    from steering_knuckle_kratos.test_strain_energy_steering_knuckle import (
        TestStrainEnergySteeringKnuckleKratos,
    )
    from trimmed_cantilever_kratos.test_trimmed_cantilever_kratos import (
        TestTrimmedCantileverKratos,
    )


def pyqueso_test_suite():
    test_suite = unittest.TestSuite()
    if kratos_available:
        test_suite.addTest(
            unittest.TestLoader().loadTestsFromTestCase(TestGGQCantileverKratos)
        )
        test_suite.addTest(
            unittest.TestLoader().loadTestsFromTestCase(TestTrimmedCantileverKratos)
        )
        test_suite.addTest(
            unittest.TestLoader().loadTestsFromTestCase(
                TestStrainEnergySteeringKnuckleKratos
            )
        )
        test_suite.addTest(
            unittest.TestLoader().loadTestsFromTestCase(TestBoundaryConditionsKratos)
        )
    else:
        print("Warning :: Tests with KratosMultiphysics dependencies are skipped.")

    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGGQTube))
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGGQ1d))
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestBSplineVolume))
    test_suite.addTest(
        unittest.TestLoader().loadTestsFromTestCase(TestBoundaryConditions)
    )
    test_suite.addTest(
        unittest.TestLoader().loadTestsFromTestCase(TestSettingsContainer)
    )
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSettingsParser))
    test_suite.addTest(
        unittest.TestLoader().loadTestsFromTestCase(TestMultiComponentModel)
    )
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestKratosAnalysis))

    return test_suite


def main():
    # Enables to open files relative to project_root
    project_root = Path(__file__).resolve().parent.parent.parent
    os.chdir(project_root)

    test_suite = pyqueso_test_suite()
    runner = unittest.TextTestRunner()
    result = runner.run(test_suite)

    sys.exit(not result.wasSuccessful())


if __name__ == "__main__":
    main()
