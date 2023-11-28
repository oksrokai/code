#!/usr/bin/env python3
from . import mupq
from .interface import M4Settings, M4

import sys


def testvectors():
    with M4() as m4:
        test = mupq.TestVectors(M4Settings(), m4)
        test.test_all(sys.argv[1:])
