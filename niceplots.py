#!/usr/bin/env python
# encoding: utf-8
"""
niceplots.py
Created by Vincent Noel on 2010-11-16, LMD/CNRS.
"""

import matplotlib.pyplot as plt
import matplotlib.collections as collections
import matplotlib.dates as mdates
import numpy as np

def infobox(text, alpha=0.7, location='top left'):
    ax = plt.gca()
    props = dict(boxstyle='round', facecolor='wheat', alpha=alpha)
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
        
    ax.text(x, y, text, transform=ax.transAxes, va='top', ha=ha, bbox=props)

def bar(x, h, color='#002299', show_stats=False, normed=False):
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

def hist(x, bins=None):
    h, xe = np.histogram(x, bins=bins)
    bar(xe, h)

def main():
    pass

if __name__ == '__main__':
    main()

def ygrid(color='grey'):
    plt.gca().yaxis.grid(True, color=color)

def xgrid(color='grey'):
    plt.gca().xaxis.grid(True, color=color)
    
def shade_x_areas(ax, x, areas, color='grey'):
    '''
    shade areas delimited by horizontal points
    '''
    
    ymin, ymax = ax.get_ylim()
    collection = collections.BrokenBarHCollection.span_where(x, ymin=ymin, ymax=ymax, where=areas, color=color, alpha=0.2)
    ax.add_collection(collection)
    
def shade_dates_areas(ax, dates, areas, color='grey'):
    x = mdates.date2num(dates)
    shade_x_areas(ax, x, areas, color=color)
    
def yaxis_season_dates(ax, days=[1,8,16,24]):
    yax = ax.yaxis
    yax.axis_date()
    yax.set_major_locator(mdates.DayLocator(bymonthday=days))
    yax.set_major_formatter(mdates.DateFormatter('%m-%d'))
    
def axis_season_dates(ax, days=[1,8,16,24]):
    xax = ax.xaxis
    xax.axis_date()
    xax.set_major_locator(mdates.DayLocator(bymonthday=days))
    xax.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.gcf().autofmt_xdate()
    