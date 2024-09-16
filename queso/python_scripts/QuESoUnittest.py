# Unittest import
import unittest

class QuESoTestCase(unittest.TestCase):
    ''' QuESoTestCase interface to derive from  unittest.TestCase in order
        to allow to add customized assert-Functions. '''

    def assertListsAlmostEqual(self, ListA, ListB, Places):
        ''' Check if ListA and ListB are almost equal.
        '''
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)
        for i, (a,b) in enumerate(zip(ListA, ListB)):
            error_msg_not_equal = "Test failed since elements with id " + str(i) + " of: " + str(ListA)
            error_msg_not_equal += " and " + str(ListB) + " do not match up to " + str(Places) + " places."
            self.assertAlmostEqual(a, b, Places, error_msg_not_equal)

    def assertListsEqual(self, ListA, ListB):
        ''' Check if ListA and ListB are equal.
        '''
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)
        for i, (a,b) in enumerate(zip(ListA, ListB)):
            error_msg_not_equal = "Test failed since elements with id " + str(i) + " of: " + str(ListA)
            error_msg_not_equal += " and " + str(ListB) + " are not equal."
            self.assertEqual(a, b, error_msg_not_equal)

    def assertListsSize(self, ListA, ListB):
        ''' Check if ListA and ListB have the same size.
        '''
        error_msg_size = "Test failed since the sizes of: " + str(ListA) + " and " + str(ListB) + " do not match."
        self.assertEqual(len(ListA), len(ListB), error_msg_size)



