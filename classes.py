from config import *

class Bin(object):
  def __init__(self, i):
    self.i = i
    self._domain = Domain(self, self, 1.0)

  def valid_begin(self):
    return self.get_root_domain().begin == self
  def valid_end(self):
    return self.get_root_domain().end == self
  def set_domain(self, domain):
    self._domain.set_parent(domain)
    self._domain = domain
  def get_root_domain(self):
    return self._domain.get_root()
  def __repr__(self):
    return 'BIN-%d' % self.i

class Domain(object):
  def __init__(self, begin, end, val, color=None):
    # Both inclusive
    self.begin = begin
    self.end = end 
    self.val = val
    self.parent = None
    self.color = color
    self.zscore = 0.0
    self.size = end.i + 1 - begin.i

  def get_begin(self):
    return self.begin.i

  def get_end(self):
    return self.end.i

  def __repr__(self):
    return 'Domain (%s, %s): %f' % (repr(self.begin), repr(self.end), self.val)

  # Reverse because the heap normally returns the smallest element
  def __lt__(self, other):
    if other is None: return False
    return self.val > other.val
  def __le__(self, other):
    if other is None: return False
    return self.val >= other.val
  def __eq__(self, other):
    if other is None: return False
    return self.val == other.val and self.begin == other.begin and self.end == other.end
  def __ne__(self, other):
    if other is None: return False
    return self.val != other.val or self.begin != other.begin or self.end != other.end
  def __gt__(self, other):
    if other is None: return False
    return self.val < other.val
  def __ge__(self, other):
    if other is None: return False
    return self.val <= other.val

  def chose(self):
    self.begin.set_domain(self)
    self.end.set_domain(self)

  def get_root(self):
    if self.parent is None:
      return self
    return self.parent.get_root()

  def set_parent(self, parent):
    if self == parent:
      print self
      sys.exit(0)
    self.parent = parent

  def valid(self):
    return self.begin.valid_begin() and self.end.valid_end() and self.begin.get_root_domain() != self.end.get_root_domain()

  def top(self):
    return self.parent is None

  def binify(self):
    self.begin.i /= BIN_SIZE
    self.end.i /= BIN_SIZE
