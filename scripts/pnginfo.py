#!/usr/bin/env python
#encoding: utf-8

from PIL import Image

def main(figname):

    im = Image.open(figname)
    print im.info


if __name__=='__main__':
    import plac
    plac.call(main)