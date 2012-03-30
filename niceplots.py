#!/usr/bin/env python
# encoding: utf-8
"""
niceplots.py

Created by Vincent Noel on 2010-11-16, LMD/CNRS.

plotting utilities

"""

import matplotlib.pyplot as plt
import matplotlib.collections as collections
import matplotlib.dates as mdates
import matplotlib.font_manager as fm
import numpy as np

mundofont = '/users/noel/.fonts/MundoSansStd.otf'
mundofontbold = '/users/noel/.fonts/MundoSansStd-Med.otf'
myfont = {}
myfont['small'] = fm.FontProperties(fname=mundofont, size=12)
myfont['med'] = fm.FontProperties(fname=mundofont, size=14)
myfont['big'] = fm.FontProperties(fname=mundofontbold, size=18)

elev_cmap = 'RdBu_r'


def savefig(figname, author='VNoel, LMD/CNRS', title=None):
    '''
    Save figure
    if format is png or pdf (based on the extension), save
    the calling script full path as a tag + author
    '''

    import os, sys

    # save script called for figure generation
    source = os.path.abspath(sys.argv[0])

    if figname.endswith('pdf'):
    
        from matplotlib.backends.backend_pdf import PdfPages
        
        pdf = PdfPages(figname)
        pdf.savefig()
        d = pdf.infodict()
        
        d['Subject'] = source        
        if author is not None:
            d['Author'] = author
        
        pdf.close()
        
    elif figname.endswith('png'):
    
        from PIL import Image, PngImagePlugin

        plt.savefig(figname)
        im = Image.open(figname)

        meta = PngImagePlugin.PngInfo()
        meta.add_text('Source', source)
        if author is not None:
            meta.add_text('Author', author)

        im.save(figname, 'PNG', pnginfo=meta)
        
    else:
    
        plt.savefig(figname)
        

def legend(numpoints=2):
    leg = plt.legend(prop=myfont['med'], numpoints=numpoints)
    plt.setp(leg.get_frame(), lw=0.5)


def segcmap(basecmap='RdBu_r', nlev=5):
    import matplotlib.cm as cm
    palette = cm.get_cmap(basecmap, nlev)
    return palette


def beautify_colorbar(cb, title=None):
    ylabels = [yt.get_text() for yt in cb.ax.get_yticklabels()]
    cb.ax.set_yticklabels(ylabels, fontproperties=myfont['small'])
    if title:
        cb.set_label(title, fontproperties=myfont['small'])
    

def beautify_axis(ax):
    fig = plt.gcf()
    fig.canvas.draw()
    xlabels = [xt.get_text() for xt in ax.get_xticklabels()]
    ylabels = [yt.get_text() for yt in ax.get_yticklabels()]
    ax.set_xticklabels(xlabels, fontproperties=myfont['small'])
    ax.set_yticklabels(ylabels, fontproperties=myfont['small'])
    plt.setp(ax.xaxis.label, fontproperties=myfont['med'])
    plt.setp(ax.yaxis.label, fontproperties=myfont['med'])
    plt.setp(ax.title, fontproperties=myfont['big'])
        

def beautify_axis_font(fig=None, ax=None, cb=None, leg=None):
    '''
    sets various fonts to Mundo Sans
    '''

    if fig is None:
        fig = plt.gcf()
    fig.canvas.draw()
    if ax:
        beautify_axis(ax)
    # if leg:
    #     beautify_legend(leg)
    # doesn't work
    if cb:
        beautify_colorbar(cb)

def infobox(text, alpha=0.7, location='top left'):
    """
    plots a box containing text in a corner of the active figure.
    arguments:
        text: the text to display
        alpha: transparency of the box
        location: either 'top left' or 'top right'
    """
    ax = plt.gca()
    props = dict(boxstyle='round', facecolor='wheat', alpha=alpha, lw=0.5)
    if location == 'top left':
        ha = 'left'
        x = 0.05
        y = 0.95
    elif location == 'top right':
        ha = 'right'
        x = 0.95
        y = 0.95
    else:
        raise Exception, 'This location is not supported yet'
        
    ax.text(x, y, text, transform=ax.transAxes, va='top', ha=ha, bbox=props, fontproperties=myfont['med'])

def bar(x, h, color='#002299', show_stats=False, normed=False):
    """
    plots a nice histogram.
    x, h:  the histogram bins and histogram heights, as returned by np.histogram.
    color: color of the histogram bars
    show_stats: adds an infobox in a corner of the histogram displaying the average and standard deviation of the distribution
    normed: if True, the distribution will be normed so its integration gives 1.
    """
    w = (x[1]-x[0])*0.8
    delta = (x[1]-x[0])*0.1
    
    if normed:
        h =  1. * h / np.sum(h)
    
    plt.bar(x[:-1]+delta, h, width=w, color=color, ec='none')
    plt.xlim(np.min(x), np.max(x))
    plt.ylim(0, np.max(h) * 1.1)
    ax = plt.gca()
    ax.yaxis.grid(True, color='white')
    
    if show_stats:
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
        normdist = h / h.sum()
        avg = np.sum(normdist * x[:-1])
        ax.text(0.05, 0.95, '$\mu = %5.2f $' % avg, transform=ax.transAxes, va='top', bbox=props)

def hist(x, bins=20):
    '''
    shortcut to create a distribution and plot it.
    '''
    h, xe = np.histogram(x, bins=bins)
    bar(xe, h)

def ygrid(color='grey'):
    '''
    Activates the grid on the y-axis for the active axis.
    '''
    plt.gca().yaxis.grid(True, color=color)

def xgrid(color='grey'):
    '''
    Activates the grid on the x-axis for active axis.
    '''
    plt.gca().xaxis.grid(True, color=color)
    
def shade_x_areas(ax, x, areas, color='grey'):
    '''
    shade areas delimited by horizontal points
    '''
    
    ymin, ymax = ax.get_ylim()
    collection = collections.BrokenBarHCollection.span_where(x, ymin=ymin, ymax=ymax, where=areas, color=color, alpha=0.2)
    ax.add_collection(collection)
    
def shade_dates_areas(ax, dates, areas, color='grey'):
    '''
    shade date areas on plot x axis
    dates = list of dates
    areas = boolean flag of dates to shade
    
    '''
    x = mdates.date2num(dates)
    shade_x_areas(ax, x, areas, color=color)
    
def yaxis_season_dates(ax, days=[1,8,16,24]):
    '''
    Sets tick marks on the y-axis nicely chosen for display time series over a season.
    '''
    yax = ax.yaxis
    yax.axis_date()
    yax.set_major_locator(mdates.DayLocator(bymonthday=days))
    yax.set_major_formatter(mdates.DateFormatter('%m-%d'))
    
    
def axis_season_dates(ax, days=[1,8,16,24]):
    '''
    Sets tick marks on the x-axis nicely chosen for display time series over a season.
    '''
    xax = ax.xaxis
    xax.axis_date()
    xax.set_major_locator(mdates.DayLocator(bymonthday=days))
    xax.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.gcf().autofmt_xdate()
    
def axis_set_date_format(ax, format='%m-%d'):
    # cf http://docs.python.org/library/datetime.html#strftime-strptime-behavior
    xax = ax.xaxis
    xax.axis_date()
    xax.set_major_formatter(mdates.DateFormatter(format))
    
    