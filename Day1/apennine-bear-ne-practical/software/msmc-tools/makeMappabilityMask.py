#!/usr/bin/env python

import gzip
import sys

class MaskGenerator:
  def __init__(self, filename, chr):
    self.lastCalledPos = -1
    self.lastStartPos = -1
    sys.stderr.write("making mask {}\n".format(filename))
    self.file = gzip.open(filename, "w")
    self.chr = chr
  
  # assume 1-based coordinate, output in bed format
  def addCalledPosition(self, pos):
    if self.lastCalledPos == -1:
      self.lastCalledPos = pos
      self.lastStartPos = pos
    elif pos == self.lastCalledPos + 1:
      self.lastCalledPos = pos
    else:
      self.file.write("{}\t{}\t{}\n".format(self.chr, self.lastStartPos - 1, self.lastCalledPos))
      self.lastStartPos = pos
      self.lastCalledPos = pos

  def close(self):
      if self.lastCalledPos != -1:
        self.file.write("{}\t{}\t{}\n".format(self.chr, self.lastStartPos - 1, self.lastCalledPos))
      self.file.close()

if len(sys.argv) < 3:
    print("Error: missing required argument.")
    print("Usage: makeMappabilityMask.py <input_file> <output_prefix>")
    sys.exit(1)

input_file = sys.argv[1]
output_prefix = sys.argv[2]

chromCounter = 0
with open(input_file, "r") as f:
  for line in f:
    if line.startswith('>'):
      chromCounter += 1
      if chromCounter != 1:
        mask.close()
      chr = line.split()[0][1:]
      mask = MaskGenerator("{0}.{1}.mask.bed.gz".format(output_prefix, chr), chr)
      pos = 0
      continue
    for c in line.strip():
      pos += 1
      if pos % 1000000 == 0:
        sys.stderr.write("processing pos:{}\n".format(pos))
      if c == "3":
        mask.addCalledPosition(pos)

mask.close()
