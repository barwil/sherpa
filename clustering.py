import numpy as np
import heapq
import Bio.Cluster
from scipy.optimize import curve_fit
import scipy.cluster
import scipy.stats
import random
from matplotlib.colors import colorConverter, Normalize, SymLogNorm

from util import smooth, debug, topify, sort_domains
from util import clip, clip_and_blur, read_domains, topify, print_domains, remove_empty_domains
from util import domains_affinity
from config import *
from classes import Bin, Domain
import matplotlib.pyplot as plt
import sys
from copy import copy

from IPython.core.debugger import Tracer

COLORS = ['WHITE', 'RED', 'YELLOW', 'GREEN', 'BLUE', 'BLACK', 'ORANGE', 'PURPLE',
          'BROWN', 'PINK', 'MAGENTA', 'CYAN', 'CHOCOLATE', 'TURQUOISE', 'SPRINGGREEN']

def akeike(model, xs, ys, params):
  result_ys = model(xs, *params)
  errs = ys - result_ys
  distrib = scipy.stats.norm(loc=np.mean(errs), scale=np.std(errs))
  return len(params) - np.sum(np.log(1 - distrib.cdf(np.abs(errs))))

def cont_max_dist_clustering(dist_matrix, corr_thr=0.7, top_only=True, arr=None, err_thr=5):
  # Returns None if nothing is to be merged.
  bins = [Bin(i) for i in range(dist_matrix.shape[0])]
  heap = [Domain(bins[i], bins[i + 1], val) 
      for i, val in enumerate(dist_matrix.diagonal(1)) if not np.isnan(val) and not np.ma.is_masked(val)]
  heapq.heapify(heap)
  selected_lengths = []
  domains = []
  inter_domain_variance = []
  prev_elem = None
  while heap:
    elem = heapq.heappop(heap)
    if not elem.valid():
      continue
    elem.chose()
    if arr is not None:
      square = np.ma.filled(
          arr[elem.get_begin():elem.get_end()+1, elem.get_begin():elem.get_end()+1],
          np.nan)
      var = np.nanvar(square)
      inter_domain_variance.append((np.nanvar(square), np.nanmean(square)))
    selected_lengths.append(elem.val)
    #if elem.val < corr_thr:
      #break
    domains.append(elem)
    prev_elem = elem
    left_bin_number = elem.get_begin() - 1
    if left_bin_number >= 0:
      left_begin = bins[left_bin_number].get_root_domain().begin
      candidate = Domain(left_begin, elem.end, dist_matrix[left_begin.i, elem.end.i])
      if not np.ma.is_masked(candidate.val):
        heapq.heappush(heap, candidate)
    right_bin_number = elem.end.i + 1
    if right_bin_number < len(bins):
      right_end = bins[right_bin_number].get_root_domain().end
      candidate = Domain(elem.begin, right_end, dist_matrix[elem.begin.i, right_end.i])
      if not np.ma.is_masked(candidate.val):
        heapq.heappush(heap, candidate)
  #return [(bin.get_root_domain().begin.i, bin.get_root_domain().end.i) for bin in bins if bin.valid_begin()] 


  min_corr = 0.2
  lengths = np.array(selected_lengths)
  max_corr = 0.97*lengths[0]
  max_limit = np.searchsorted(-lengths, -max_corr)
  min_limit = np.searchsorted(-lengths, -min_corr)

  if not min_limit:
    min_limit = len(lengths) - 1

  break_step = 5

  def two_linear(thr):
    # return x * a1 + b1 if x < thr else x * a2 + b2
    # But the "if else" construction is not broadcastable
    def inner(x, a1, b1, a2, b2):
      ver1 = x * a1 + b1
      ver2 = x * a2 + b2
      mult1 = 1.0 * (-np.sign(x - thr) + 1) / 2
      mult2 = 1.0 * (np.sign(x - thr) + 1) / 2
      return ver1 * mult1 + ver2 + mult2
    return inner

  def one_linear(x, a, b):
    return 1.0 * a * x + b

  xs = np.array(range(len(selected_lengths)))
  lengths_limited = lengths[max_limit:min_limit]
  xs_limited = xs[max_limit:min_limit]

  tested_breaks = range(max_limit, min_limit, break_step)
  errs = np.zeros_like(tested_breaks, dtype=np.float32)
  for i, fixed_break in enumerate(tested_breaks):
    params, _ = curve_fit(two_linear(fixed_break), xs_limited, lengths_limited)
    estimated = two_linear(fixed_break)(xs_limited, *params)
    err = np.mean((estimated - lengths_limited) ** 2)
    errs[i] = err

  result_limit = tested_breaks[np.argmin(errs)]
  # For visualization purposes
  params, _ = curve_fit(two_linear(result_limit), xs_limited, lengths_limited)
  params_one, _ = curve_fit(one_linear, xs_limited, lengths_limited)
  err_one = np.mean((one_linear(xs_limited, *params_one) - lengths_limited) ** 2)
  debug('Two-linear is %f times better (%f / %f)' % (err_one / np.min(errs), err_one, np.min(errs)))
  if err_one / np.min(errs) < err_thr:
    debug('Merging below error threshold')
    return None
  debug('Two-linear has AIC: %f, one-linear has AIC: %f' % (
    akeike(two_linear(result_limit), xs_limited, lengths_limited, params),
    akeike(one_linear, xs_limited, lengths_limited, params_one)))

  smoothed_lengths = smooth(lengths, 20)
  if SHOW_PLOTS:
    plt.plot(xs, smoothed_lengths)
    plt.axvline(result_limit, c='red')
    plt.axvline(max_limit, c='red')
    plt.axvline(min_limit, c='red')
    plt.plot(xs_limited, two_linear(result_limit)(xs_limited, *params), c='purple')
    plt.plot(xs_limited, one_linear(xs_limited, *params_one), c='purple')
    ax2 = plt.twinx()
    ax2.plot(tested_breaks, errs)

    plt.show()
    plt.figure()
    gradient1 = smooth(np.gradient(smoothed_lengths), 20)
    plt.plot(range(len(selected_lengths)), gradient1)
    plt.axvline(result_limit, c='red')
    if arr is not None:
      variances = [x[0] for x in inter_domain_variance]
      means = [x[1] for x in inter_domain_variance]
      plt.show()
      plt.figure()
      plt.plot(range(len(selected_lengths)), variances)
      #plt.twinx()
      #plt.plot(range(len(selected_lengths)), means, c='green')
      #plt.twinx()
      #plt.plot(range(len(selected_lengths)), [var/mean for (var, mean) in inter_domain_variance], c='red')
      plt.axvline(result_limit, c='red')
    plt.show()

  domains = domains[:result_limit]
  #result = [(domain.begin.i, domain.end.i) for domain in domains]
  result = copy(domains)
  result.extend([Domain(bins[i], bins[i], 0.0) for i in range(dist_matrix.shape[0])])
  result = sort_domains(result)
  if top_only:
    result = topify(result)
  return result

def k_medoids(dist_matrix, k=5, big_iters=300):
  size = dist_matrix.shape[0]
  min_dist = np.inf
  result = []
  for y in range(big_iters):
    medoids = np.array(random.sample(range(size), k))
    changing = True
    iters = 0
    iters_limit = 200
    curr_sum = 0
    while changing and iters < iters_limit:
      changing = False
      assigned = np.argmin(dist_matrix[medoids], axis=0)
      assigned[medoids] = range(k)
      for i in range(k):
        #print medoids
        #print assigned
        old_medoid = medoids[i]
        indices = np.array(range(size))[assigned == i]
        #Tracer()()
        #print indices
        medoids[i] = indices[np.argmin(np.sum(dist_matrix[np.meshgrid(indices, indices)], axis=0))]
        #Tracer()()
        if medoids[i] != old_medoid:
          changing = True
      #print medoids
      #print assigned
      curr_sum = np.sum(dist_matrix[np.array(range(size)), medoids[assigned]])
      iters += 1
    print curr_sum
    if curr_sum < min_dist:
      result = assigned
      min_dist = curr_sum
      #print 'Winning curr_sum'
      #print result
  return result
  #return [np.array(range(size))[assigned == i] for i in range(k)]

def colorize_doms(arr, domains):
  max_clusters = 100
  non_empty_cluster_thr = 10
  arr = np.ma.fix_invalid(arr)
  non_empty_domains = remove_empty_domains(arr, topify(domains))
  doms_affinity = domains_affinity(arr, non_empty_domains)
  dist_matrix = - np.clip(doms_affinity, 0, 10) + 10
  doms_size = dist_matrix.shape[0]
  clusters = None
  hist = None
  Z = scipy.cluster.hierarchy.complete(dist_matrix)
  for t in range(3, max_clusters):
    clusters = scipy.cluster.hierarchy.fcluster(Z, t=t, criterion='maxclust') - 1
    hist = np.histogram(clusters, bins=range(t + 1))[0]
    if np.max(hist) < (doms_size / 3):
      break
  debug('Found %d clusters' % t) 
  debug('Division of the clusters: %s' % repr(clusters))
  debug('Clusters\' sizes: %s' % repr(hist))
  non_empty_clusters = np.where(hist >= non_empty_cluster_thr)[0]
  non_empty_clusters = np.argsort(hist)[::-1]
  #print non_empty_clusters
  debug('Non empty clusters: %s' % repr(non_empty_clusters))
  scipy.cluster.hierarchy.dendrogram(Z)
  plt.show()
  #show(block=False)
  #figure()
  # From defective clusters have color WHITE, not clustered have GREY.
  for dom in domains:
    #dom.color = 'GREY'
    dom.color = -1
  for i, dom in enumerate(non_empty_domains):
    #if clusters[i] in non_empty_clusters:
      # colorConverter.to_rgba
      #dom.color = (
      #    COLORS[
      #      np.where(non_empty_clusters == clusters[i])[0][0] + 1
      #    ])
    dom.color = np.where(non_empty_clusters == clusters[i])[0][0]
    #else:
      #dom.color = COLORS[0]
      #print 'wrong'
  #Tracer()()

if __name__ == '__main__':
  if len(sys.argv) < 3:
    print 'arr doms'
    sys.exit(1)
  blur = True
  arr = np.load(sys.argv[1])
  if blur:
    arr = clip_and_blur(arr)
  else:
    arr = clip(arr)
  doms = topify(read_domains(open(sys.argv[2], 'r')))
  colorize_doms(arr, doms)
  print_domains(doms)
