#!/usr/bin/env python3
import serial.tools.list_ports
import sys

# macos_dev = "/dev/tty.usbserial"
# linux_dev = "/dev/ttyUSB0"

ports = serial.tools.list_ports.comports()
for port in ports:
    print(port.description)
    if 'VCP' in port.description or 'Serial' in port.description:
        dev = serial.Serial(port.device, 115200, timeout=1)
        print("Listen on port {} from {}".format(port.device, port.description))
        break

print("> Returned data: ", file=sys.stderr)

while True:
    x = dev.read()
    sys.stdout.buffer.write(x)
    sys.stdout.flush()
