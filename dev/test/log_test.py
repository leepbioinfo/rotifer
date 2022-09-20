#!/usr/bin/env python3
from rotifer.core.log import log
import unittest


class LogTestCase(unittest.TestCase):

    def test1(self):
        self.assertEqual('a' ,'b',msg = 'return should be a nor b')

if __name__ == '__main__':
    unittest.main()
