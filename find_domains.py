import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import scipy.ndimage
import pdb
from copy import copy
from IPython.core.debugger import Tracer

import scipy.stats

from util import heatmap, clip_and_blur, debug, clip, heatmap_notlog
from util import print_domains, for_each_diagonal, map_over_diagonals, smooth
from util import matrix_from_list, domain_contacts
from clustering import cont_max_dist_clustering
from classes import Domain, Bin
from config import *

def sum_inside(cumulative, start, end):
  sum = cumulative[end, end]
  sum -= cumulative[start - 1, end] if start > 0 else 0
  sum -= cumulative[end, start - 1] if start > 0 else 0
  sum += cumulative[start - 1, start - 1] if start > 0 else 0
  return sum

#<not-debuggable>
#Good luck.
#Doing data manipulation in python sucks badly. So everything needs to be 
#delegated to numpy. Which is usually possible, but very hacky.

def square_sums(cumulative):
  # A lot of magic happens. Basically it's sum_inside implemented on an array.
  size = cumulative.shape[0]
  plus1 = np.tile(np.diagonal(cumulative), (size, 1))
  plus2 = np.insert(np.tile(
    np.diagonal(cumulative).reshape((size, 1)), size), 0, 0, axis=0)[:-1, :]
  minus1 = np.insert(cumulative, 0, 0, axis=0)[:-1, :]
  minus2 = np.insert(np.transpose(cumulative), 0, 0, axis=0)[:-1, :]
  return np.triu(plus1 + plus2 - minus1 - minus2)

#</not-debuggable>

def z(q_val, mean, stddev):
  return (q_val - mean) / stddev

def q(sum, dom_size, gamma=2):
  return sum / (np.power(dom_size + 1, gamma))

def z_score_matrix(arr, gamma=2.0):
  size = arr.shape[0]
  cum = np.cumsum(np.cumsum(arr, axis=0), axis=1)
  square_sum_matrix = square_sums(cum)
  distance_matrix = matrix_from_list(range(size))
  q_matrix = q(square_sum_matrix, distance_matrix, gamma)
  means = for_each_diagonal(q_matrix, np.nanmean)

  stddevs = for_each_diagonal(q_matrix, lambda x: np.power(np.nanvar(x), 0.5))
  means_matrix = matrix_from_list(means)
  stddev_matrix = matrix_from_list(stddevs)

  z_matrix = np.nan_to_num(z(q_matrix, means_matrix, stddev_matrix))

  #z_matrix = np.nan_to_num(q_matrix - means_matrix)

  return z_matrix

def dynamic(quality_matrix):
  size = quality_matrix.shape[0]
  optimal_score = np.empty(size)
  optimal_score.fill(-np.inf)
  optimal_score[0] = 0
  previous_end = np.empty(size)
  previous_end.fill(-1)
  domain_defining = np.empty(size)
  np.set_printoptions(threshold=np.nan)
  for i in range(size):
    cand_nodomain = np.nanargmax(optimal_score)
    with_domain = optimal_score + quality_matrix[:, i]
    cand_domain = np.nanargmax(with_domain)
    if optimal_score[cand_nodomain] > with_domain[cand_domain]:
      domain_defining[i] = 0
      previous_end[i] = cand_nodomain
      optimal_score[i] = optimal_score[cand_nodomain]
    else:
      domain_defining[i] = 1
      previous_end[i] = cand_domain
      optimal_score[i] = with_domain[cand_domain]
  current_end = size - 2 
  result = []
  while current_end > 0:
    if domain_defining[current_end] == 1:
      result.append(Domain(Bin(previous_end[current_end]), Bin(current_end), 0))
    current_end = previous_end[current_end]
  return result[::-1]

def greedy_hierarchy(original_arr, gamma=2, delta=0.5):
  all_result = []
  to_do = [(0, original_arr.shape[0])]
  while to_do:
    start, end = to_do.pop()
    size = end - start
    if size < 15:
      continue
    arr = original_arr[start:end, start:end]
    z_matrix = z_score_matrix(arr, gamma) - 1
    z_matrix = Normalize()(np.clip(z_matrix, 0, np.inf)) - 0.001

    #imshow(z_matrix, cmap=cm.jet, norm=Normalize(), interpolation='nearest', origin='lower')
    #colorbar()
    #show(block=False)
    #figure()

    z_cum = np.cumsum(np.cumsum(np.clip(z_matrix, 0, np.inf), axis=0), axis=1)
    z_squares = square_sums(z_cum) #/ ((distance_matrix + 1) ** 2)
    
    to_the_end = np.tile(z_squares[:, -1], (size, 1))
    from_the_beggining = np.tile(z_squares[0, :].reshape(size, 1), (1, size))
    possibilities = np.triu(to_the_end + from_the_beggining + z_squares)

    scaled_delta = delta / (size ** 2)

    normalized_pos = Normalize()(possibilities)
    quality_matrix = np.triu(z_matrix + delta * normalized_pos)

    result = np.unravel_index(np.nanargmax(quality_matrix), quality_matrix.shape)

    #imshow(normalized_pos, cmap=cm.jet, norm=Normalize(), interpolation='nearest', origin='lower')
    #colorbar()
    #show(block=False)
    #figure()

    debug('On %d:%d, result is %d, %d' % (
        start, start + size, result[0] + start, result[1] + start))

    #imshow(quality_matrix, cmap=cm.jet, norm=Normalize(), interpolation='nearest', origin='lower')
    #colorbar()
    #show()

    if z_matrix[result] < 0:
      debug(z_matrix[result])
      continue

    if result[0] == 0 and result[1] == size -1:
      # Best solution is the whole input
      continue

    all_result.append((result[0] + start, result[1] + start))

    to_do.append((start, result[0] + start - 1))
    to_do.append((result[0] + start, result[1] + start))
    to_do.append((result[1] + start + 1, end))
  
  return all_result

def correlation(arr, dom_limit=100):
  size = arr.shape[0]
  lower_ones = np.tril(np.ones((size, size)), - 1)
  pears = np.ma.array(np.zeros((size, size)), fill_value=1)
  pval = np.zeros((size, size))
  arr = np.ma.filled(arr, 0)
  for i in range(size):
    for dom_size in range(min(dom_limit, size - i)):
      one = arr[i] + arr[:, i]
      two = arr[i + dom_size] + arr[:, i + dom_size]
      if np.all(one == 0) or np.all(two == 0):
        pears[i, i + dom_size] = np.ma.masked
        pval[i, i + dom_size] = np.ma.masked
      else:
        pear, p = scipy.stats.pearsonr(
            np.delete(one, range(i, i + dom_size + 1)), np.delete(two, range(i, i + dom_size + 1)))
        pears[i, i + dom_size] = pear
        pval[i, i + dom_size] = p
  return lower_ones + pears, lower_ones + pval

def cumulative_min_correlation(arr, dom_limit=100):
  corr, pvals = correlation(arr, dom_limit)
  xs = np.array([], dtype=np.float32)
  ys = np.array([], dtype=np.float32)
  for i in range(1, dom_limit):
    xs = np.append(xs, np.clip(np.diag(arr, i) + 1, -np.inf, np.inf))
    #xs = np.append(xs, np.log(np.clip(np.diag(arr, i) / np.nanmean(np.diag(arr, i)), -np.inf, np.inf)))
    ys = np.append(ys, np.diag(corr, i))
  if SHOW_PLOTS:
    hexbin(xs, ys, xscale='log', norm=LogNorm())
    colorbar()
    show()
    #heatmap_notlog(corr, block=False)
    #figure()
  filled_acc = np.minimum.accumulate(np.flipud(
    np.minimum.accumulate(np.flipud(np.ma.filled(corr)), axis=0)), axis=1)
  return np.ma.array(filled_acc, mask=corr.mask)

def absolutize_domains(second_order, first_order):
  # Second order is a list of domains relative to first_order, which need to be absolute
  return [Domain(Bin(first_order[second_order_dom.get_begin()].get_begin()),
                 Bin(first_order[second_order_dom.get_end()].get_end()), 0.0)
    for second_order_dom in second_order]

def greedy_correlation_merging(arr, corr_thr, max_dom_size, rounds):
  err_thr=5
  size = arr.shape[0]
  temp_arr = np.ma.copy(arr)
  prev_domains = [Domain(Bin(x), Bin(x), 0.0) for x in range(size)]
  all_domains = []
  for i in range(rounds):
    debug('Starting round %d of hierarchical clusterization' % i)
    if SHOW_PLOTS:
      heatmap(temp_arr)
    cum_corr = cumulative_min_correlation(temp_arr, max_dom_size)
    debug('Done cumulative_min_correlation')
    #heatmap_notlog(cum_corr)
    temp_domains = cont_max_dist_clustering(cum_corr, corr_thr=corr_thr,err_thr=err_thr)
    if not temp_domains:
      break
    debug('Done cont_max_dist_clustering')
    absolute_domains = absolutize_domains(temp_domains, prev_domains) 
    all_domains.extend(absolute_domains)
    temp_arr = domain_contacts(arr, absolute_domains)
    debug('Done domain_contact')
    prev_domains = absolute_domains
  return all_domains

if __name__ == '__main__':
  if len(sys.argv) < 2:
    print 'array gamma'
    print 'If gamma not specified, then 2 is assumed'
  if len(sys.argv) >= 3:
    gamma = float(sys.argv[2])
  else:
    gamma = 2
  arr = np.load(sys.argv[1])
  size = arr.shape[0]
  arr = clip(arr)

  corr_thr = 0.5
  max_dom_size = 500
  rounds = 8
  print_domains(greedy_correlation_merging(arr, corr_thr, max_dom_size, rounds))
