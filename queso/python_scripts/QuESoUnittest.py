# Unittest import
import unittest

class QuESoTestCase(unittest.TestCase):
    ''' QuESoTestCase interface to derives from  unittest.TestCase in order
        to allow to add customized assert-Function '''

    def assertListsAlmostEqual(self, ListA, ListB, Places):
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)
        for i, (a,b) in enumerate(zip(ListA, ListB)):
            error_msg_not_equal = "Test failed since elements with id " + str(i) + " of: " + str(ListA)
            error_msg_not_equal += " and " + str(ListB) + " do not match up to " + str(Places) + " places."
            self.assertAlmostEqual(a, b, Places, error_msg_not_equal)

    def assertListsEqual(self, ListA, ListB):
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)
        for a,b in zip(ListA, ListB):
            error_msg_not_equal = "Test failed since elements with id " + str(i) + " of: " + str(ListA)
            error_msg_not_equal += " and " + str(ListB) + " are nto equal."
            self.assertEqual(a, b, error_msg_not_equal)

    def assertListsSize(self, ListA, ListB):
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)



