#!/usr/bin/python

import sys
import re
import gzip
import fastq

RFileName = sys.argv[1]
BarcodeFileName = sys.argv[2]
OutFileName = sys.argv[3]

RFile = gzip.open(RFileName)
BarcodeFile = gzip.open(BarcodeFileName)
OutFile = open(OutFileName, "w")

while 1:

    RLine = RFile.readline()
    if not RLine:
        break

    assert RLine[0] == '@'
    Label = RLine.strip()[1:]

    RLine = RFile.readline()
    Seq = RLine.strip()

    RLine = RFile.readline()
    Plus = RLine.strip()
    assert Plus[0] == '+'

    RLine = RFile.readline()
    Qual = RLine.strip()

    BCLine = BarcodeFile.readline()
    assert BCLine[0] == '@'

    BCLine = BarcodeFile.readline()
    BC = BCLine.strip()
    
    BCLine = BarcodeFile.readline()
    assert Plus[0] == '+'
    BCLine = BarcodeFile.readline()

    Label = re.sub(" .*", "", Label) + "BC=" + BC

    fastq.WriteRec(OutFile, Label, Seq, Qual)

