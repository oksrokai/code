#!/usr/bin/env python3
from scripts.convert_benchmarks import convert_benchmarks
from scripts.test import test
from scripts.testvectors import testvectors
from scripts.benchmarks import benchmarks
from scripts.build_everything import build_everything

from rich.console import Console
from random import randrange

funcs = {
    0: test,
    1: testvectors,
    2: build_everything,
    3: benchmarks,
    4: convert_benchmarks,
}

if __name__ == "__main__":
    console = Console()
    console.rule("[bold bright_green]Available Functions")
    for i in funcs:
        console.print("[{}] {}".format(i, funcs[i].__name__), style="color({})".format(randrange(17, 230)))
    try:
        index = int(console.input("Choose function (default is 0): "))
    except:
        index = 0
    funcs[index]()

    