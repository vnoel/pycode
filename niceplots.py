#!/usr/bin/env python
# encoding: utf-8

"""
niceplots.py

Created by Vincent Noel on 2010-11-16, LMD/CNRS.

plotting utilities

"""

import matplotlib.pyplot as plt
import matplotlib.collections as collections
from matplotlib.ticker import FuncFormatter
import matplotlib.dates as mdates
import matplotlib.font_manager as fm
import numpy as np
import os

homedir = os.path.expanduser('~')
mundofont = homedir + '/.fonts/MundoSansStd.otf'
mundofontbold = homedir + '/.fonts/MundoSansStd-Med.otf'
myfont = {
    'small': fm.FontProperties(fname=mundofont, size=12),
    'med': fm.FontProperties(fname=mundofont, size=14),
    'big': fm.FontProperties(fname=mundofontbold, size=18)
}
elev_cmap = 'RdBu_r'


def savefig(figname, author='VNoel, LMD/CNRS'):
    """
    Save figure
    if format is png or pdf (based on the extension), save
    the calling script full path as a tag + author
    """

    import os
    import sys

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


def legend(numpoints=2, **kwargs):
    leg = plt.legend(prop=myfont['med'], fancybox=True, labelspacing=0.2, numpoints=numpoints, **kwargs)
    plt.setp(leg.get_frame(), lw=0.5)


def segcmap(basecmap='jet', nlev=5):
    # use : pcolor(cmap=segcmap)
    import matplotlib.cm as cm

    palette = cm.get_cmap(basecmap, nlev)
    return palette


def segcmap_elev(basecmap='RdBu_r', nlev=5):
    # use : pcolor(cmap=segcmap)
    import matplotlib.cm as cm

    palette = cm.get_cmap(basecmap, nlev)
    return palette


def beautify_colorbar(cb, title=None, ticks=None):
    ylabels = [yt.get_text() for yt in cb.ax.get_yticklabels()]
    cb.ax.set_yticklabels(ylabels, fontproperties=myfont['small'])
    if title:
        cb.set_label(title, fontproperties=myfont['small'])
    if ticks is not None:
        cb.set_ticks(ticks)


def colorbar():
    cb = plt.colorbar()
    beautify_colorbar(cb)


def lon_formatter_func(x, pos):
    if x == 0:
        return '0°'
    elif x > 0:
        return '%3.0f°E' % x
    else:
        return '%3.0f°W' % -x


lon_formatter = FuncFormatter(lon_formatter_func)


def lat_formatter_func(x, pos):
    if x == 0:
        return '0°'
    elif x > 0:
        return '%3.0f°N' % x
    else:
        return '%3.0f°S' % -x


lat_formatter = FuncFormatter(lat_formatter_func)


def cb_right(title=None):
    ax = plt.axes([0.92, 0.25, 0.02, 0.5])
    cb = plt.colorbar(cax=ax)
    beautify_colorbar(cb, title=title)


def beautify_title(ax):
    plt.setp(ax.title, fontproperties=myfont['big'])
    # beautify_colorbar(cb, title='H$_2$O [ppmv]')


def beautify_axis(ax):
    fig = plt.gcf()
    fig.canvas.draw()
    xlabels = [xt.get_text() for xt in ax.get_xticklabels()]
    ylabels = [yt.get_text() for yt in ax.get_yticklabels()]
    ax.set_xticklabels(xlabels, fontproperties=myfont['small'])
    ax.set_yticklabels(ylabels, fontproperties=myfont['small'])
    plt.setp(ax.xaxis.label, fontproperties=myfont['med'])
    plt.setp(ax.yaxis.label, fontproperties=myfont['med'])
    beautify_title(ax)


def beautify_axis_font(fig=None, ax=None, cb=None):
    """
    sets various fonts to Mundo Sans
    """

    if fig is None:
        fig = plt.gcf()
    fig.canvas.draw()
    if ax:
        beautify_axis(ax)
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
        raise Exception('This location is not supported yet')

    ax.text(x, y, text, transform=ax.transAxes, va='top', ha=ha, bbox=props, fontproperties=myfont['med'])


def bar(x, h, color='#002299', show_stats=False, normed=False, fix_lim=True, log=False):
    """
    plots a nice histogram.
    x, h:  the histogram bins and histogram heights, as returned by np.histogram.
    color: color of the histogram bars
    show_stats: adds an infobox in a corner of the histogram displaying the average and standard deviation of the dist.
    normed: if True, the distribution will be normed so its integration gives 1.
    """
    w = (x[1] - x[0]) * 0.8
    delta = (x[1] - x[0]) * 0.1

    if normed:
        h = 1. * h / np.sum(h)

    plt.bar(x[:-1] + delta, h, width=w, color=color, ec='none', log=log)
    if fix_lim:
        plt.xlim(np.min(x), np.max(x))
        plt.ylim(0, np.max(h) * 1.1)
    ax = plt.gca()
    ax.yaxis.grid(True, color='white')

    if show_stats:
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
        normdist = h / h.sum()
        avg = np.sum(normdist * x[:-1])
        ax.text(0.05, 0.95, '$\mu = %5.2f $' % avg, transform=ax.transAxes, va='top', bbox=props)


def hist(x, bins=20, log=False):
    """
    shortcut to create a distribution and plot it.
    """
    h, xe = np.histogram(x, bins=bins)
    bar(xe, h, log=log)


def ygrid(color='grey'):
    """
    Activates the grid on the y-axis for the active axis.
    """
    plt.gca().yaxis.grid(True, color=color)


def xgrid(color='grey'):
    """
    Activates the grid on the x-axis for active axis.
    """
    plt.gca().xaxis.grid(True, color=color)


def shade_x_areas(ax, x, areas, color='grey'):
    """
    shade areas delimited by horizontal points
    """

    ymin, ymax = ax.get_ylim()
    collection = collections.BrokenBarHCollection.span_where(x, ymin=ymin, ymax=ymax, where=areas, color=color,
                                                             alpha=0.2)
    ax.add_collection(collection)


def shade_dates_areas(ax, dates, areas, color='grey'):
    """
    shade date areas on plot x axis
    dates = list of dates
    areas = boolean flag of dates to shade
    
    """
    x = mdates.date2num(dates)
    shade_x_areas(ax, x, areas, color=color)


def xaxis_day(ax):
    ax.xaxis.axis_date()

    ax.xaxis.set_minor_locator(mdates.HourLocator(12))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))

    maxlabels = [tick.label1 for tick in ax.xaxis.get_major_ticks()]
    plt.setp(maxlabels, rotation=30, fontsize=10, ha='right', weight='bold')

    minlabels = [tick.label1 for tick in ax.xaxis.get_minor_ticks()]
    plt.setp(minlabels, rotation=30, fontsize=10, ha='right')

    ax.xaxis.grid(True, 'major')


def xaxis_year_month(ax):
    ax.xaxis.axis_date()

    ax.xaxis.set_minor_locator(mdates.MonthLocator(list(range(2, 13))))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    maxlabels = [tick.label1 for tick in ax.xaxis.get_major_ticks()]
    plt.setp(maxlabels, rotation=30, fontsize=10, ha='right', weight='bold')

    minlabels = [tick.label1 for tick in ax.xaxis.get_minor_ticks()]
    plt.setp(minlabels, rotation=30, fontsize=10, ha='right')

    ax.xaxis.grid(True, 'minor')


def yaxis_season_dates(ax, days=None):
    """
    Sets tick marks on the y-axis nicely chosen for display time series over a season.
    """
    if not days: days = [1, 8, 16, 24]
    yax = ax.yaxis
    yax.axis_date()
    yax.set_major_locator(mdates.DayLocator(bymonthday=days))
    yax.set_major_formatter(mdates.DateFormatter('%m-%d'))


def axis_season_dates(ax, days=None):
    """
    Sets tick marks on the x-axis nicely chosen for display time series over a season.
    """
    if not days: days = [1, 8, 16, 24]
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
