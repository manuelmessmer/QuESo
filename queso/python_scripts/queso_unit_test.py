from typing import List, Any
# Unittest import
import unittest

class QuESoTestCase(unittest.TestCase):
    """Extended unittest.TestCase for QuESo-specific testing.

    Provides customized assertion functions for comparing lists with numerical tolerances.
    """
    def assertListsAlmostEqual(self,
            ListA: List[float],
            ListB: List[float],
            Places: int
        ) -> None:
        """Assert that two lists are approximately equal element-wise.

        Args:
            ListA (List[float]): First list of numbers.
            ListB (List[float]): Second list of numbers.
            Places (int): Number of decimal places for comparison.

        Raises:
            AssertionError: If the lists differ in size or any corresponding elements differ beyond the given precision.
        """
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)
        for i, (a,b) in enumerate(zip(ListA, ListB)):
            error_msg_not_equal = "Test failed since elements with id " + str(i) + " of: " + str(ListA)
            error_msg_not_equal += " and " + str(ListB) + " do not match up to " + str(Places) + " places."
            self.assertAlmostEqual(a, b, Places, error_msg_not_equal)

    def assertListsEqual(self,
            ListA: List[Any],
            ListB: List[Any]
        ) -> None:
        """Assert that two lists are exactly equal element-wise.

        Args:
            ListA (List[Any]): First list.
            ListB (List[Any]): Second list.

        Raises:
            AssertionError: If the lists differ in size or any corresponding elements are not equal.
        """
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)
        for i, (a,b) in enumerate(zip(ListA, ListB)):
            error_msg_not_equal = "Test failed since elements with id " + str(i) + " of: " + str(ListA)
            error_msg_not_equal += " and " + str(ListB) + " are not equal."
            self.assertEqual(a, b, error_msg_not_equal)

    def assertListsSize(self,
            ListA: List[Any],
            ListB: List[Any]
        ) -> None:
        """Assert that two lists have the same size.

        Args:
            ListA (List[Any]): First list.
            ListB (List[Any]): Second list.

        Raises:
            AssertionError: If the lists have different sizes.
        """
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)



