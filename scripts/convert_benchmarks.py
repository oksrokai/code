#!/usr/bin/env python3
import sys
from . import mupq


def usage():
    print("Usage: python3 convert_benchmarks.py csv|md")
    sys.exit(1)


def convert_benchmarks():
    if len(sys.argv) != 2:
        usage()
    if sys.argv[1] == "csv":
        converter = mupq.CsvConverter()
    elif sys.argv[1] == "md":
        converter = mupq.MarkdownConverter()
    else:
        usage()
    converter.convert()
