import numpy as npy
import pylab
import matplotlib.cbook as cbook

def multihist(xvals, bins=10, normed=0, bottom=None,
    		   align='edge', orientation='vertical', width=None,
    		   log=False, type='overlap',gap=None, patch_kwargs=None, labels=None, **kwargs):

      #some integrity checks up front
      if type == 'bi' and len(xvals) != 2:
    	  raise ValueError('need exactly two data sets for "bi" multihist: %d given' % len(xvals))
      if patch_kwargs is not None and len(patch_kwargs) != len(xvals):
    	  raise ValueError('need same number of patch kwargs and data sets')

      #calculate the common bins, more or less stolen from numpy.histogram
      xvals = [npy.asarray(x).ravel() for x in xvals]
      if not npy.iterable(bins):
    	  mn = float(min([x.min() for x in xvals]))
    	  mx = float(max([x.max() for x in xvals]))
    	  if mn == mx:
    		  mn -= 0.5
    		  mx += 0.5
    	  bins = npy.linspace(mn, mx, bins, endpoint=False)

      #make the histograms using the common bins
      xn = []
      for x in xvals:
    	  n, bins2 = npy.histogram(x, bins, range=None, normed=normed)
    	  xn.append(n)

      #build the patches parameters depending on type argument
      if width is None: width = 0.9*(bins[1]-bins[0])
      delta = 0
      offset = 0
      paint_width = width
      stay_on_top = True
      if type == 'beside':
    	  if npy.iterable(width):
    		  raise ValueError('no sequence of widths allowed for "beside" multihist')
    	  width /= len(xn)
    	  delta = width
    	  if align == 'edge':
    		  offset = 0
    	  elif align == 'center':
    		  offset = ((len(xn) / -2.0 + 0.5) * width)
    	  else:
    		  raise ValueError('invalid alignment: %s' % align)
    	  if gap is None:
    		  gap = 0
    	  paint_width = width - gap
      elif type == 'bi':
    	  stay_on_top = False
      elif type != 'overlap':
    	  raise ValueError('invalid multihist type: %s' % type)

      #build the patches
      patch_list = []
      on_top = True
      for n in xn:
    	  obins = [b + offset for b in bins]
    	  if on_top:
    		  rn = n
    	  else:
    		  rn = [-v for v in n]
    	  if orientation == 'horizontal':
    		  patches = pylab.barh(obins, rn, height=paint_width, left=bottom,
    								  align=align, log=log)
    	  elif orientation == 'vertical':
			  patches = pylab.bar(obins[:-1], rn, width=paint_width, bottom=bottom,
    							  align=align, log=log, linewidth=0)
    	  else:
    		  raise ValueError('invalid orientation: %s' % orientation)
    	  patch_list.append(cbook.silent_list('Patch', patches))
    	  offset += delta
    	  on_top = on_top and stay_on_top

      for i in range(len(patch_list)):
    	  if patch_kwargs == None:
    		  kwa = kwargs
    	  else:
    		  kwa = patch_kwargs[i]
    		  if labels:
				  lab = labels[i]
    	  for p in patch_list[i]:
    		  p.update(kwa)
		  if labels:
		      patch_list[i][0].update({"label":lab})

      return xn, bins, patch_list
