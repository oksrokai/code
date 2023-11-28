#!/usr/bin/env python3
import subprocess

from serial import Serial
from serial.tools import list_ports

from rich.console import Console
from random import randrange

from . import mupq


class M4Settings(mupq.PlatformSettings):
    #: Specify folders to include
    scheme_folders = [  # mupq.PlatformSettings.scheme_folders + [
        ('lattice-m4', 'schemes/crypto_kem', ''),
    ]

    #: List of dicts, in each dict specify (Scheme class) attributes of the
    #: scheme with values, if all attributes match the scheme is skipped.
    skip_list = (
        # {'scheme': 'falcon-1024-tree', 'implementation': 'opt-leaktime'}
    )


class M4(mupq.Platform):

    def __enter__(self):
        console = Console()
        ports = list_ports.comports()
        availablePorts = []
        index = 0
        console.rule("[bold bright_green]Available VCP or Serial Ports")
        for port in ports:
            if 'VCP' in port.description or 'Serial' in port.description:
                console.print("[{}] {}".format(index, port.description), style="color({})".format(randrange(17, 230)))
                availablePorts.append(port)
                index += 1
        if not availablePorts:
            console.print("No available ports!")
            exit(-1)
        else:
            try:
                index = int(console.input("Choose port (default is 0): "))
                self._dev = Serial(availablePorts[index].device, 115200, timeout=10)
            except:
                self._dev = Serial(availablePorts[0].device, 115200, timeout=10)
        return super().__enter__()

    def __exit__(self,*args, **kwargs):
        self._dev.close()
        return super().__exit__(*args, **kwargs)

    def device(self):
        return self._dev

    def flash(self, binary_path):
        super().flash(binary_path)
        subprocess.check_call(
            ["st-flash", "--reset", "write", binary_path, "0x8000000"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
