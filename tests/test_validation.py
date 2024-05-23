# tests/test_validation.py


import unittest
from unittest.mock import patch
from importnb import Notebook

with Notebook(): 
    import code

class TestGetNonZeroInteger(unittest.TestCase):

    @patch('builtins.input', side_effect=['5'])
    def test_get_nonzero_integer(self, mock_input):
        expected_output = 3
        actual_output = code.get_nonzero_integer("Enter a non-zero integer: ")
        self.assertEqual(actual_output, expected_output)

if __name__ == '__main__':
    unittest.main()