#!/usr/bin/env python
#encoding: utf-8

"""
Function to merge existing PDF figures into a multi-page PDF.
"""

import glob
import os

import pyPdf


def merge_pdf_with_basename(basename, cleanup=False):
    
    output = pyPdf.PdfFileWriter()

    files = glob.glob(basename + '_*.pdf')
    if len(files) == 0:
        print('No pdf file matching %s criteria' % basename)
        return
        
    print('Merging ', files, ' in %s.pdf' % basename)

    for f in files:
        input1 = pyPdf.PdfFileReader(file(f, 'rb'))
        output.addPage(input1.getPage(0))

    fileout = file(basename+'.pdf', 'wb')
    output.write(fileout)
    fileout.close()
    
    if cleanup:
        print('Deleting %d files' % (len(files)))
        for f in files:
            os.remove(f)