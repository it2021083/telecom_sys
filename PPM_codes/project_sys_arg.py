#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 22:34:40 2024

@author: pp
"""

#Complete project

import numpy as np
import matplotlib.pyplot as plt
import sys



# Maps the graycode words to the k symbols
def gray_to_k_mapping(M):
    gray_symbols = graycode(int(np.log2(M))) # Creates all the possible words of the 
                                             # Gray_code with log2(M)-bit words 
    k_mapping = {symbol: k for k, symbol in enumerate(gray_symbols)} # Creates a dictionary 
                                                                     # where the keys are the 
                                                                     # symbols and the values are their 
                                                                     # corresponding GrayCode words.
    return k_mapping




# Coverts the binary sequence to a sequence of symbols
def binary_to_ppm(binary_sequence, M):
    gray_code_sequence = binary_to_gray(binary_sequence, M) # Converts the binary sequence to a gracode sequence
    k_mapping = gray_to_k_mapping(M) # Creates a mapping of graycode-symbols
    bits_per_symbol = int(np.log2(M))  # Number of bits per symbol
    ppm_symbols = []

    # Iterate through the Gray code sequence and create PPM symbols
    for i in range(0, len(gray_code_sequence), bits_per_symbol):
        symbol_bits = gray_code_sequence[i:i + bits_per_symbol].zfill(bits_per_symbol) # Adds 0s at the end of the graycode to ensure it consists of whole log2(M)-bit words
        ppm_symbol = k_mapping[symbol_bits] # Mapping of the specific graycode word to a symbol rusing the k_mapping map
        ppm_symbols.append(ppm_symbol) # Adds the new symbol at the end of the ppm_symbols array

    return ppm_symbols # Returns the sequence of symbols as an array



# Creates the waveform(Amplitude in relevance with time) of a symbol
def symbol_ppm_wave(t, Tp, k, Ts):
    time_resolution = 1e-12
    time_points_per_symbol = int(Ts / time_resolution) # Defines the time points per symbol
    total_time_points = time_points_per_symbol # Total time of the wave
    t_per_symbol = np.linspace(0, Ts, total_time_points) # Creates a time vector for the symbol
    A = np.zeros_like(t_per_symbol) # Creates an array with 0s that represents the Amplitude of the symbol
    
    # Set the amplitude to 1 to a specific position inside the Ts(symbol duration) that is defined by the k*Tp <= t < (k+1)*Tp
    pulse_condition = np.logical_and(k * Tp <= t_per_symbol, t_per_symbol < (k + 1) * Tp)
    A[pulse_condition] = 1
    
    return A # Returns the amplitude array that is basicaly the waveform of the symbol



# Creates the waveform of a sequence of symbols
def generate_ppm_waveform(M, Ts, symbol_sequence):
    Tp = pulse_duration(Ts, M) # Sets the pulse duration Tp
    time_resolution = 1e-12
    time_points_per_symbol = int(Ts / time_resolution) # Defines the time points per symbol
    total_time_points = time_points_per_symbol * len(symbol_sequence) # Total time of the wave
    t = np.linspace(0, Ts * len(symbol_sequence), total_time_points) # Creates a time vector for the whole symbol sequence
    
    waveform = np.zeros_like(t) # Creates an array with 0s that represents the wave
    
    # Iterates the symbol sequence and creates the waveform of every symbol and append them to the waveform array
    for symbol_order, k in enumerate(symbol_sequence):
        start_index = int(symbol_order * time_points_per_symbol) # Define the start time of the symbole
        end_index = int((symbol_order + 1) * time_points_per_symbol) # Define the end time of the symbole(the start of the next)
        
        symbol_waveform = symbol_ppm_wave(t[start_index:end_index], Tp, k, Ts) # Create the waveform(Amplitude) of the symbol
        waveform[start_index:end_index] = symbol_waveform #Add the waveform of the symbol for the defined time(use of indecies)
    
    return t, waveform# waveform Return the whole duration of the wave(t) and the wave(Amplitudes)




# Set Ts(symbol duration with known Rb)
def symbol_duration(M):
    Rb = 10**9
    Ts = np.log2(M) / Rb  # Convert to seconds
    return Ts

# Set Tp(pulse duration)
def pulse_duration(Ts, M):
    Tp = Ts / M
    return Tp

# Convert text to a binary sequence
def ascii_to_binary(text):
    binary_result = ""
    for char in text:
        binary_result += format(ord(char), '08b')
    return binary_result

# Create GrayCode of M-bit words
def graycode(M):
    if M == 1:
        g = ['0', '1']
    elif M > 1:
        gs = graycode(M-1)
        gsr = gs[::-1]
        gs0 = ['0' + x for x in gs]
        gs1 = ['1' + x for x in gsr]
        g = gs0 + gs1
    return g

# Convert a binary sequence to graycode sequence in regard of M(level of PPM or else log2(M)-bit words)
def binary_to_gray(binary_piece, M):
    log2M = int(np.log2(M)) # Set log2M that is how many bits per word 
    gray_code_output = "" # Create an empty grayCode string

    # iterate the binary sequence 
    for i in range(0, len(binary_piece), log2M):
        binary_chunk = binary_piece[i:i + log2M].ljust(log2M, '0')
        gray_piece = bin(int(binary_chunk, 2) ^ (int(binary_chunk, 2) >> 1))[2:].zfill(log2M)
        gray_code_output += gray_piece

    return gray_code_output

# Values of M to consider
M_values = [2, 4, 8, 16]

# Check if the English name is provided as a command-line argument
if len(sys.argv) < 2:
    print("Usage: python script_name.py <]>")
    sys.exit(1)

text_input = sys.argv[1]

# Loop through each value of M
for M in M_values:
    

    log2M = int(np.log2(M))
    binary_input = ascii_to_binary(text_input)
    gray_code_output = binary_to_gray(binary_input, M)
    k_mapping = gray_to_k_mapping(M)

    # Parameters for plotting
    Ts = symbol_duration(M)
    Tp = pulse_duration(Ts, M)

    ppm_symbols = binary_to_ppm(binary_input, M)


    print("\n\n")
    # Print additional information on console
    print(f'Gray Code words for M={M}: {graycode(int(np.log2(M)))}')
    print('Gray to k Mapping:')
    for symbol, k_value in k_mapping.items():
        print(f'{symbol} -> {k_value}')
    print('')

    print(f'log: {int(np.log2(M))}')
    print(f'Binary Input: {binary_input}')
    print(f'Gray Code Output: {gray_code_output}')
    print(f'PPM Symbols: {ppm_symbols}')
    print('')
    print(f'Ts: {Ts}')
    print(f'Tp: {Tp}')

    # Plot the waveform for the M-PPM transmission of "P"
    t, waveform = generate_ppm_waveform(M, Ts, ppm_symbols)

    # Plot the complete waveform
    plt.figure(figsize=(50, 10))
    plt.step(t, waveform, where='post')
    plt.title(f'Complete PPM Waveform for M={M}')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')

    # Set ticks every Tp seconds
    time_resolution = Tp
    ticks = np.arange(0, len(ppm_symbols) * Ts, time_resolution)
    plt.xticks(ticks)
    
    # Save the plot as a PNG file with M in the filename
    plt.savefig(f'ppm_waveform_M{M}.png')

    plt.show(block=False)

    # Print time durations where each symbol is 1 and 0
    for symbol_order, k in enumerate(ppm_symbols):
        start_time = symbol_order * Ts
        end_time = (symbol_order + 1) * Ts
        print(f"Symbol {k}: from {start_time} s to {end_time} s")

    print("\n")
    
    
    
    
    
plt.show()
