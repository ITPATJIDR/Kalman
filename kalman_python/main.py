import numpy as np
import time
import math

# Constants
K = 4
SAMPLING_FREQUENCY = 40000
SAMPLING_BUFFER_SIZE = 1024
MAXIMUM_SHIFT = 80

signal_frequencies = [500, 1000, 2000, 1500]
signal = np.zeros((3, SAMPLING_BUFFER_SIZE))
index_sample = 0


def generate_mask(mask_length, period_length):
    mask = np.sin(np.arange(mask_length) * 2 * np.pi / period_length)
    return mask


mask_length = [SAMPLING_BUFFER_SIZE - MAXIMUM_SHIFT for _ in range(K)]
mask_signal = [generate_mask(
    mask_length[i], SAMPLING_FREQUENCY / signal_frequencies[i]) for i in range(K)]


def read_inputs():
    global index_sample, signal

    index_sample += 1
    if index_sample >= SAMPLING_BUFFER_SIZE - 1:
        index_sample = 0

    signal[0][index_sample] = np.random.randint(0, 1024)
    signal[1][index_sample] = np.random.randint(0, 1024)
    signal[2][index_sample] = np.random.randint(0, 1024)
    print(signal)


sampling_interval = 1.0 / SAMPLING_FREQUENCY


def start_sampling():
    while True:
        read_inputs()
        time.sleep(sampling_interval)

if __name__ == "__main__":
    start_sampling()
