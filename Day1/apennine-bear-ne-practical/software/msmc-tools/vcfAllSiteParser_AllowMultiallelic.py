#!/usr/bin/env python

# The vcfAllSiteParser.py script ignores multiallelic sites, which are therefore considered uncalled
# irrespective of whether the genotype is missing or not. In contrast, the
# vcfAllSiteParser_AllowMultiallelic.py script regards multiallelic sites with a non-missing
# genotype as called. Such sites will be included in the output mask file. In addition, multiallelic
# sites where the genotype contains at least one alternative allele will be printed to the output
# VCF file. Note that, since the script no longer applies any constraints to the ALT field of the
# input VCF file, care should be taken to remove unwanted variants such as indels beforehand.

# Please address questions and comments about the vcfAllSiteParser_AllowMultiallelic.py script to
# @IrisLiesbethRuesinkBueno.

import sys
import gzip
import re

class MaskGenerator:
  def __init__(self, filename, chr):
    self.lastCalledPos = -1
    self.lastStartPos = -1
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
    print("too few arguments:")
    print("Usage: ./vcfCaller.py <chrom> <mask_out> > <vcf_out>")
    print("Reading VCF with all called sites, including hom-ref from stdin")
    sys.exit(1)

chr = sys.argv[1]
mask_filename = sys.argv[2]

mask = MaskGenerator(mask_filename, chr)

lastPos = 0
line_cnt = 0
for line in sys.stdin:
  if line[0] == '#':
    print line,
    continue
  fields = line.strip().split('\t')
  pos = int(fields[1])
  refAllele = fields[3]
  altAllele = fields[4]
  info = fields[7]
  genotypes = fields[9]
  if line_cnt % 10000 == 0:
    sys.stderr.write("parsing position {}\n".format(pos))
  line_cnt += 1
  if re.match("^[ACTGactg]$", refAllele):
    if genotypes[0] != '.' and genotypes[2] != '.':
      mask.addCalledPosition(pos)
      if genotypes[0] != '0' or genotypes[2] != '0':
        print line,
  
mask.close()
