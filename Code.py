# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 20:00:31 2019

@author: Theodore Michels
"""

import numpy as np
import math
from Bio import SeqIO
from PIL import Image

# FASTA file to load
inputSequence = 'seq.fasta'

# Define the aspect ratio of the canvas
aspectX = 16
aspectY = 9

# Associate each base with a color
colors = {
    'A': (255, 0, 0),
    'C': (0, 255, 0),
    'G': (0, 0, 255),
    'T': (255, 255, 0)
}

totalBases = 0

print("Calculating total number of bases...")
for seq_record in SeqIO.parse(inputSequence, 'fasta'):
    # Keep count of the total number of bases
    seqLength = len(seq_record)
    totalBases += seqLength

print("Total number of bases: " + str(totalBases))
print("Rendering at " + str(aspectX) + ":" + str(aspectY))

# Calculate the canvas width and height
width = math.ceil(math.sqrt((aspectX * totalBases) / aspectY))
height = math.ceil(totalBases / width)

print("Output resolution: " + str(width) + "x" + str(height))

# Create empty 3D array in which to store color data
data = np.zeros((height, width, 3), dtype=np.uint8)

count = 0
x = 0
y = 0
lastY = y
rtol = False

# Parse the input file...
for seq_record in SeqIO.parse(inputSequence, 'fasta'):
    print("Drawing sequence " + str(seq_record.id) + " with length of " + str(seqLength) + " bases...")

    # Loop through each character in the sequence
    for char in seq_record.seq:
        current = char.upper()  # Convert to uppercase to match keys

        data[y, x] = colors.get(current, (0, 0, 0))

        count += 1
        if (rtol):
            x = width - (count % width) - 1
        else:
            x = count % width

        y = int(count / width)

        if (y != lastY):
            rtol = not rtol

        lastY = y

    print("Finished drawing sequence " + str(seq_record.id))

# Save and display the image
img = Image.fromarray(data, 'RGB')
img.save('canvas.bmp')
img.show()
