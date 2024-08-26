import serial

import serial.tools.list_ports

ports = serial.tools.list_ports.comports()
for port in ports:
    print(port.device)


# Set up the serial connection (adjust the port and baudrate as needed)
ser = serial.Serial('/dev/cu.usbmodem153480801', 115200, timeout=1)  # Replace '/dev/ttyUSB0' with your port
def read_signal():
     while True:
        try:
            line = ser.readline().decode('utf-8').strip()  # Read a line from the serial port
            if line:
                signal_value = int(line)  # Convert the line to an integer
                print(f'Received signal[0]: {signal_value}')
        except ValueError:
            pass  # Ignore lines that can't be converted to an integer

if __name__ == "__main__":
    read_signal()

