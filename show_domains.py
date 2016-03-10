import sys
import numpy as np
import matplotlib.colors as pltcol
from matplotlib import cm
from matplotlib.colors import colorConverter, Normalize, SymLogNorm, cnames
from matplotlib.pyplot import *
import random

from util import clip_and_blur, read_domains, clip, topify, remove_empty_domains, interpolated
from util import for_each_diagonal, smooth, matrix_from_list, heatmap, heatmap_notlog
from find_domains import z_score_matrix 
from util import domains_affinity, inter_domain_contacts
from util import read_cc, fix_ticks, flip_to_diagonal

from clustering import k_medoids

from IPython.core.debugger import Tracer
from scipy.optimize import curve_fit
import scipy.cluster

BORDER_COLOR = 'WHITE'
COLORS = {
    '0': 'WHITE',
    'A': 'BLACK',
    'B': 'RED', 
    'C': 'GREEN',
    'D': 'PURPLE',
    'E': 'BROWN',
    'F': 'TURQUOISE'
    }
COLOR_LIST = random.sample(cnames, 100)
#print COLOR_LIST
DOM_THR = 0

def colored_with_affinity(arr, domains, affinity=False, chr=None, cc=None):
  arr = np.ma.fix_invalid(arr)
  # PuBu
  colored = cm.jet(pltcol.LogNorm()(arr))
  normalizer = Normalize(vmin=-1, vmax=3)
  non_empty_domains = remove_empty_domains(arr, topify(domains))
  if affinity:
    doms_affinity = domains_affinity(arr, non_empty_domains)
    heatmap_notlog(np.clip(doms_affinity, 0, 10))
    doms_affinity_log = np.log(doms_affinity)
    affinity_thr = np.nanmean(doms_affinity_log) + np.nanstd(doms_affinity_log)
    affinity_thr = np.exp(affinity_thr)
    affinity_thr = np.nanmean(doms_affinity) + 1 * np.nanstd(doms_affinity)
    affinity_normalizer = Normalize(vmin=affinity_thr, vmax=10)

  for i, dom in enumerate(domains):
    begin, end = dom.get_begin(), dom.get_end()
    if end - begin < DOM_THR:
      continue
    color = dom.color or BORDER_COLOR
    try:
      color_num = int(dom.color)
      if color_num >=0 :
        color = COLOR_LIST[color_num]
      else:
        color = COLORS['0']
    except:
      pass
    if cc:
      color = COLORS[cc.get((chr, dom.color), '0')]
    color = colorConverter.to_rgba(color)
    for i in range(begin, end + 1):
      if i != end:
        colored[i, end] = color
      if i != begin:
        colored[begin, i] = color
    #colored[begin, begin] = colorConverter.to_rgba('black')
    #colored[end, end] = colorConverter.to_rgba('black')
  if affinity:
    for i, dom in enumerate(non_empty_domains):
      begin, end = dom.get_begin(), dom.get_end()
      for j, dom2 in enumerate(non_empty_domains):
        if j <= i:
          continue
        begin2, end2 = dom2.get_begin(), dom2.get_end()
        if doms_affinity[i, j] < affinity_thr:
          continue
        color = cm.gnuplot(affinity_normalizer(doms_affinity[i, j]))
        for k in range(begin, end + 1):
          colored[k, begin2] = color
          colored[k, end2] = color
        for k in range(begin2, end2 + 1):
          colored[begin, k] = color
          colored[end, k] = color
  return colored

if __name__ == '__main__':
  if len(sys.argv) < 3:
    print '[-cc cc_file1 [-cc2 cc_file2]] arr domains (arr2 domains2)'
  cc1 = None
  cc2 = None
  if sys.argv[1] == '-cc':
    cc1 = read_cc(open(sys.argv[2], 'r'))
    sys.argv = sys.argv[2:]
    if sys.argv[1] == '-cc2':
      cc2 = read_cc(open(sys.argv[2], 'r'))
      sys.argv = sys.argv[2:]
  blur = True
  comparing = len(sys.argv) > 3
  colorize = False
  show_affinity = True
  arr = np.load(sys.argv[1])
  chr1 = sys.argv[1].split('/')[-1].split('.')[0].split('_')[0]
  if blur:
    arr = clip_and_blur(arr)
  else:
    arr = clip(arr)
  if comparing:
    arr = np.triu(arr)
  doms = read_domains(open(sys.argv[2], 'r'))
  result = colored_with_affinity(arr, doms, chr=chr1, cc=cc1)
  if comparing:
    arr2 = np.load(sys.argv[3])
    chr2 = sys.argv[3].split('/')[-1].split('.')[0].split('_')[0]
    doms2 = read_domains(open(sys.argv[4], 'r'))
    if blur:
      arr2 = clip_and_blur(arr2)
    else:
      arr2 = clip(arr2)
    arr2 = np.triu(arr2, 1)
    result2 = colored_with_affinity(arr2, doms2, chr=chr2, cc=cc2)
    result += np.transpose(result2, (1, 0, 2))
  else:
    result = flip_to_diagonal(result)
  figure(figsize=(8,8))
  imshow(result, origin='lower', interpolation='nearest')
  fix_ticks()
  title(sys.argv[-1])
  show()
  #savefig(sys.argv[-1]+".pdf",dpi=1200)
