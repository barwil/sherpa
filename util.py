import os
import logging
import numpy as np
from pylab import *
from matplotlib.colors import LogNorm, Normalize
from matplotlib.pyplot import *
import scipy.ndimage
from scipy.interpolate import interp2d

from config import *
from classes import Domain, Bin

from IPython.core.debugger import Tracer

logging.basicConfig(level=logging.DEBUG)

DOM_THR = 2

def gc(x):
  return x in ['G', 'C', 'g', 'c']

def valid(x):
  return x in ['A', 'G', 'C', 'T', 'a', 'g', 'c', 't']

def gc_content(seq, resolution=100):
  # Returns an int from [0, resolution]
  all = 0
  gc_count = 0
  for letter in seq:
    if valid(letter):
      all += 1
    if gc(letter):
      gc_count += 1
  if all == 0:
    return 0
  else:
    return resolution * gc_count / all

valid_chromosomes = ['2R', '2L', '3R', '3L', '4', 'X']

def valid_chromosome(x):
  return x in valid_chromosomes

def debug(*args, **kwargs):
  logging.debug(*args, **kwargs)

def info(*args, **kwargs):
  logging.info(*args, **kwargs)

def error(*args, **kwargs):
  logging.error(*args, **kwargs)

def for_each_diagonal(arr, f):
  size = arr.shape[0]
  result = np.empty(size)
  for i in range(size):
    result[i] = f(np.ma.diagonal(arr, i))
  return result

def map_over_diagonals(arr, f):
  size = arr.shape[0]
  result = np.zeros_like(arr)
  for i in range(size):
    result += np.diagflat(f(np.diagonal(arr, i)), i)
  return result

def matrix_from_list(l):
  # Creates a_{i, j} = l[j - i]
  return np.triu(np.insert(
      np.tile(np.array(l), (len(l), 1)), len(l), 0, axis=1).reshape(
          (len(l) + 1, len(l)))[:-1, :])

def blur(matrix, samples):
  result = np.zeros(matrix.shape).astype(np.float32)
  for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
      num_samples = samples(i, j)
      w = num_samples / 2
      sum_norm = 0.0
      for k in range(w * 2 + 1):
        for l in range(w * 2 + 1):
          norm = 1.0 / (abs(k - w) + abs(l - w) + 1)
          try:
            result[i, j] += matrix[i - w + k, j - w + l] * norm
            sum_norm += norm
          except:
            pass
      result[i, j] /= sum_norm
  return result

def fix_ticks(x_only=False):
  def f(x, pos):
    dist = 1.0 * x * BIN_SIZE
    if dist > 1000000:
      return '%.2fMbp' % (dist / 1000000)
    else:
      return '%dkbp' % (dist / 1000)
  axes().xaxis.set_major_formatter(FuncFormatter(f))
  if not x_only:
    axes().yaxis.set_major_formatter(FuncFormatter(f))

def flip_to_diagonal(a, n=200):
  size = a.shape[0]
  b = np.zeros_like(a)
  b = b[:n]
  for i in range(n):
    temp = None
    if a.ndim == 2:
      temp = np.diag(a, i)
      temp.resize((size - i,))
    elif a.ndim == 3:
      # Problem-specific
      color = [None] * 4
      for k in range(4):
        temp = np.diag(a[:,:,k], i*2)
        temp.resize((size - i))
        color[k] = temp
      temp = np.zeros((color[0].shape[0], 4))
      for k in range(4):
        temp[:,k] = color[k]
    b[i, i:] = temp
  return b

def read_cc(f):
  result = dict()
  for line in f:
    s = line.split(' ')
    result[(s[0], s[1])] = s[2][0]
  return result

def add_to_chr_array(result, to_add):
  for name, array in to_add.iteritems():
    result[name] = np.add(array, result[name]) \
        if name in result else np.copy(array)

def add_to_chr_dict(result, to_add):
  for name, dic in to_add.iteritems():
    if name not in result:
      result[name] = {}
    add_to_dict(result[name], dic)

def add_to_dict(d1, d2):
  for k, v in d2.iteritems():
    d1[k] = d1.get(k, 0) + v

def flatten_chr_array(d):
  result = None
  for chr, d1 in d.iteritems():
    if result is None:
      result = d1
    else:
      result = np.add(result, d1)
  return result

def flatten_chr_dict(d):
  result = {}
  for chr, d1 in d.iteritems():
    add_to_dict(result, d1)
  return result

def align(what, how):
  return ((what - 1) / how + 1) * how

def heatmap(arr, block=True):
  imshow(arr, cmap=cm.jet, interpolation='nearest', norm=LogNorm(), origin='lower')
  colorbar()
  show(block=block)

def heatmap_notlog(arr, block=True):
  imshow(arr, cmap=cm.jet, interpolation='nearest', norm=Normalize(), origin='lower')
  colorbar()
  show(block=block)

def sort_domains(domains):
  return sorted(domains, key=lambda dom: (dom.get_begin(), dom.get_end()))

def topify(domains):
  sorted_domains = sort_domains(domains)
  result = [sorted_domains[0]]
  for original_domain in sorted_domains:
    (original_domain_begin, original_domain_end) = original_domain.get_begin(), original_domain.get_end()
    if original_domain_begin == result[-1].get_begin():
      result[-1] = original_domain
    elif original_domain_begin > result[-1].get_end():
      result.append(original_domain)
  return result

def remove_empty_domains(arr, domains):
  result = []
  for domain in domains:
    begin, end = domain.get_begin(), domain.get_end()
    if end - begin < DOM_THR:
      continue
    #Tracer()()
    if ((type(arr.mask) == np.bool_ and arr.mask) or (type(arr.mask) != np.bool_ and np.all(arr.mask[begin:end+1]))) or \
        np.all(np.isnan(arr[begin:end+1])):
      continue
    result.append(domain)
  return result

def clip_and_blur(arr, stddevs=5, blur=1):
  arr = np.ma.masked_invalid(arr)
  mean = np.mean(arr)
  stddev = np.var(arr) ** 0.5
  np.clip(arr, 0, mean + stddevs * stddev, out=arr)
  arr = np.ma.filled(arr, 0)
  scipy.ndimage.gaussian_filter(arr, blur, output=arr)
  np.clip(arr, mean * 0.01, mean + stddevs * stddev, out=arr)
  return arr

def interpolated(arr):
  arr = clip(arr).copy()
  heatmap(arr, block=False)
  figure()
  nans = np.isnan(np.diag(arr))
  not_nan_indices = np.arange(arr.shape[0])[~nans]
  interp = interp2d(not_nan_indices, not_nan_indices,
      arr[np.meshgrid(not_nan_indices, not_nan_indices)])
  interpolated = interp(np.arange(arr.shape[0]), np.arange(arr.shape[0]))
  heatmap(interpolated, block=False)
  figure()
  arr[arr.mask] = interpolated[arr.mask]
  heatmap(arr)

def inter_domain_contacts(arr, domains):
  result = np.ma.asarray(arr.copy())
  result.mask = True
  size = arr.shape[0]
  for dom1 in domains:
    start1, temp_end1 = dom1.get_begin(), dom1.get_end()
    for dom2 in domains:
      if dom2 == dom1:
        continue
      start2, temp_end2 = dom2.get_begin(), dom2.get_end()
      end1 = temp_end1 + 1
      end2 = temp_end2 + 1
      if end1 >= size:
        end1 = size
      if end2 >= size:
        end2 = size
      m = np.mean(arr[start1:end1, start2:end2]) + np.mean(arr[start2:end2, start1:end1])
      result[start1:end1, start2:end2] = m
      result[start2:end2, start1:end1] = m
  return result

def domains_affinity(arr, domains):
  #debug('There is %d non-empty domains' % len(domains))
  by_dist = np.nan_to_num(for_each_diagonal(arr, lambda x: np.nanmean(x) if len(x) > 100 else 0))
  merged_domains = inter_domain_contacts(arr, topify(domains))
  expected = inter_domain_contacts(matrix_from_list(by_dist), topify(domains))
  expected_only_domains = domain_contacts(expected, domains)
  merged_only_domains = domain_contacts(merged_domains, domains)
  affinity_only_domains = merged_only_domains / expected_only_domains
  return affinity_only_domains

def domain_contacts(arr, domains):
  # TODO: not happy.
  arr = np.ma.fix_invalid(arr)
  weird_min = 0.01
  dom_len = len(domains)
  size = arr.shape[0]
  cum = np.cumsum(np.cumsum(np.ma.filled(arr, 0), axis=0), axis=1)
  count = np.cumsum(np.cumsum(~arr.mask, axis=0), axis=1)
  starts = np.array([dom.get_begin() for dom in domains])
  ends = np.array([dom.get_end() for dom in domains])
  starts1, starts2 = np.meshgrid(starts - 1, starts - 1)
  ends1, ends2 = np.meshgrid(ends, ends)
  np.clip(ends1, 0, size - 1, out=ends1)
  np.clip(ends2, 0, size - 1, out=ends2)

  cum = np.append(np.append(cum, [[0] * size], axis=0), [[0]] * (size + 1), axis=1)
  count = np.append(np.append(count, [[0] * size], axis=0), [[0]] * (size + 1), axis=1)
  plus1 = cum[ends1, ends2]
  plus2 = cum[starts1, starts2]
  minus1 = cum[starts1, ends2]
  minus2 = cum[ends1, starts2]

  plus3 = cum[ends2, ends1]
  plus4 = cum[starts2, starts1]
  minus3 = cum[starts2, ends1]
  minus4 = cum[ends2, starts1]

  count_plus1 = count[ends1, ends2]
  count_plus2 = count[starts1, starts2]
  count_minus1 = count[starts1, ends2]
  count_minus2 = count[ends1, starts2]
  result = np.ma.fix_invalid((plus1 + plus2 - minus1 - minus2 + plus3 + plus4 - minus3 - minus4) / \
      (count_plus1 + count_plus2 - count_minus1 - count_minus2))
  result[result < weird_min] = 0.0
  return np.ma.fix_invalid(result)
  return np.clip(result, weird_min, np.inf)


def clip(arr, stddevs=10):
  arr = np.ma.masked_invalid(arr)
  mean = np.mean(arr)
  stddev = np.var(arr) ** 0.5
  np.clip(arr, 0, mean + stddevs * stddev, out=arr)
  return arr

def print_domains(domains, chr=''):
  for i, dom in enumerate(domains):
    if dom.color is not None:
      print 'chr%s %d: %d %d %s' % (chr, i, dom.get_begin(), dom.get_end(), dom.color)
    else:
      print 'chr%s %d: %d %d' % (chr, i, dom.get_begin(), dom.get_end())

def read_domains(f):
  result = []
  sane_limit = 100000
  for line in f:
    l = line.split()
    begin, end, color = None, None, None
    if len(l) == 3:
      # Old format
      begin, end = int(l[1]), int(l[2])
    elif len(l) == 5:
      begin, end = int(l[2]), int(l[3])
      color = l[4]
    else:
      try:
        begin, end = int(l[2]), int(l[3])
        color = None
      except:
        begin, end = int(l[1]), int(l[2])
        color = l[3]
    result.append(Domain(Bin(begin), Bin(end), 0.0, color=color))
  if result and result[-1].get_end() > sane_limit:
    for dom in result:
      dom.binify()
  return result

def read_bed(f):
  result = []
  for line in f:
    l = line.split()
    result.append((int(l[1]) / BIN_SIZE, int(l[2]) / BIN_SIZE))
  return result

def isolators_from_domains(domains):
  # We're assuming domains are correct domains.
  # Returns list of list of isolators. The first one is the strongest isolators.
  sorted_domains = sort_domains(domains)
  events = []
  prev = None
  for dom in sorted_domains:
    if prev is not None and (dom.get_begin(), dom.get_end()) == (prev.get_begin(), prev.get_end()):
      continue
    (domain_begin, domain_end) = dom.get_begin(), dom.get_end()
    if domain_begin == domain_end:
      continue
    events.append((domain_begin, 0))
    events.append((domain_end, 1))
    prev = dom
  events.sort()
  isols = [[]]
  nest = 0
  for (where, what) in events:
    if what == 1:
      nest -= 1
    isols[nest].append(where)
    if what == 0:
      nest += 1
    while len(isols) < nest + 1:
      isols.append([])
  return isols

# You know nothing, Jon Snow
def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    import numpy as np    
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]
