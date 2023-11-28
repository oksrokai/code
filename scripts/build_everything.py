#!/usr/bin/env python3
import sys

from .interface import M4Settings
from . import mupq


def build_everything():
    mupq.BuildAll(M4Settings()).test_all(sys.argv[1:])
