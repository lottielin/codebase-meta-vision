#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <iterator>
#include <numeric>
#include <ostream>
#include <iomanip>
#include <vector>

// Simple 2-D Point representation
struct Point {
  int y;
  int x;
  Point(): y(0), x(0) {}
  Point(int y, int x): y(y), x(x) {}
  inline bool operator<(const Point& other) const
    { return y < other.y || (y == other.y && x < other.x); }
  inline bool operator==(const Point& other) const
    { return y == other.y && x == other.x; }
  inline bool operator!=(const Point& other) const
    { return y != other.y || x != other.x; }
  inline Point operator+(const Point& other) const
    { return Point(y + other.y, x + other.x); }
  inline Point operator-(const Point& other) const
    { return Point(y - other.y, x - other.x); }
  inline Point operator>>(int scalar) const
    { assert(scalar >= 0); return Point(y >> scalar, x >> scalar); }
  inline Point operator<<(int scalar) const
    { assert(y >= 0 && x >= 0 && scalar >= 0);
      return Point(y << scalar, x << scalar); }
  inline Point operator*(int scalar) const
    { return Point(y * scalar, x * scalar); }
  inline Point operator*(const Point& other) const
    { return Point(y * other.y, x * other.x); }
  inline Point operator/(const Point& other) const
    { return Point(y / other.y, x / other.x); }
  // Find Bi-linear interpolation points for this point
  std::pair<Point, Point> bilinear() const;
};

std::ostream& operator<<(std::ostream& out, const Point& point);

// Simple 2-D Rectangle representation and raster-scan iterator
class Tessellation;

class Rect {
public:
  Rect(int offset_y, int offset_x, int height, int width):
    min_y_(offset_y), min_x_(offset_x), max_y_(offset_y + height - 1),
    max_x_(offset_x + width - 1) {}
  Rect(const Point& top_left, int height, int width):
    min_y_(top_left.y), min_x_(top_left.x), max_y_(top_left.y + height - 1),
    max_x_(top_left.x + width - 1) {}
  Rect(const Point& top_left, const Point& bottom_right):
    min_y_(top_left.y), min_x_(top_left.x), max_y_(bottom_right.y),
    max_x_(bottom_right.x) {}
  Rect(): min_y_(0), min_x_(0), max_y_(0), max_x_(0) {}

  using tile_iterator = Tessellation;
  friend tile_iterator;
  friend std::ostream& operator<<(std::ostream&, const Rect&);

  class iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Point;
    using difference_type = long;
    using pointer = const Point*;
    using reference = Point;

    iterator(const Rect& rect, int offset_y, int offset_x):
      rect_(rect), point_(offset_y, offset_x) {}

    inline iterator& operator++() {
      if (rect_.max_x_ == point_.x++) {
        point_.x = rect_.min_x_;
        ++point_.y;
      }
      return *this;
    }
    inline iterator operator++(int)
      { iterator retval = *this; ++(*this); return retval; }
    inline bool operator==(const iterator& other) const
      { return point_ == other.point_; }
    inline bool operator!=(const iterator& other) const
      { return point_ != other.point_; }
    inline reference operator*() const { return Point(point_); }

    const Rect& rect_;
    Point point_;
  };

  inline iterator begin() const { return iterator(*this, min_y_, min_x_); }
  inline iterator end() const { return iterator(*this, max_y_ + 1, min_x_); }

  class int_iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Point;
    using difference_type = long;
    using pointer = const Point*;
    using reference = Point;

    int_iterator(const Rect& rect, int offset_y, int offset_x):
      rect_(rect), iter_(offset_x + rect.width() * offset_y) {}

    inline int_iterator& operator++() {
      ++iter_;
      return *this;
    }
    inline int_iterator operator++(int)
      { int_iterator retval = *this; ++(*this); return retval; }
    inline bool operator==(const int_iterator& other) const
      { return iter_ == other.iter_; }
    inline bool operator!=(const int_iterator& other) const
      { return iter_ != other.iter_; }
    inline reference operator*() const
      { return Point(rect_.min_x_ + iter_ / rect_.width(),
        rect_.min_y_ + iter_ % rect_.width()); }

  private:
    const Rect& rect_;
    int iter_;
  };

  inline int_iterator int_begin() const
    { return int_iterator(*this, min_y_, min_x_); }
  inline int_iterator int_end() const
    { return int_iterator(*this, max_y_ + 1, min_x_); }

  template <bool VERTICAL=false>
  class grid_iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Point;
    using difference_type = long;
    using pointer = const Point*;
    using reference = Point;

    grid_iterator(const Rect& rect, int offset_y, int offset_x,
      int stride_y, int stride_x):
      rect_(rect), point_(offset_y, offset_x),
      stride_y_(stride_y), stride_x_(stride_x) {}

    inline grid_iterator& operator++() {
      if (VERTICAL) {
        point_.y += stride_y_;
        if (rect_.max_y_ < point_.y) {
          point_.y = rect_.min_y_;
          point_.x += stride_x_;
        }
      } else {
        point_.x += stride_x_;
        if (rect_.max_x_ < point_.x) {
          point_.x = rect_.min_x_;
          point_.y += stride_y_;
        }
      }
      return *this;
    }
    inline grid_iterator operator++(int)
      { grid_iterator retval = *this; ++(*this); return retval; }
    inline bool operator==(const grid_iterator& other) const
      { return point_ == other.point_; }
    inline bool operator!=(const grid_iterator& other) const
      { return !(*this == other); }
    inline reference operator*() const { return Point(point_); }

  private:
    const Rect& rect_;
    Point point_;
    const int stride_y_;
    const int stride_x_;
  };

  inline grid_iterator<false> grid_begin(int stride_y, int stride_x) const
    { return grid_iterator<false>(*this, min_y_, min_x_, stride_y, stride_x); }
  inline grid_iterator<false> grid_end(int stride_y, int stride_x) const
    { return grid_iterator<false>(*this,
      min_y_ + stride_y * (1 + (max_y_ - min_y_) / stride_y),
      min_x_, stride_y, stride_x); }

  inline grid_iterator<true> vgrid_begin(int stride_y, int stride_x) const
    { return grid_iterator<true>(*this, min_y_, min_x_, stride_y, stride_x); }
  inline grid_iterator<true> vgrid_end(int stride_y, int stride_x) const
    { return grid_iterator<true>(*this,
      min_y_, min_x_ + stride_x * (1 + (max_x_ - min_x_) / stride_x),
      stride_y, stride_x); }

  tile_iterator tile_begin(size_t height, size_t width,
    bool vertical_first=false) const;

  inline Point offset() const { return Point(min_y_, min_x_); }
  inline Point transverse() const { return Point(max_y_, max_x_); }
  inline int offset_y() const { return min_y_; }
  inline int offset_x() const { return min_x_; }
  inline int width() const { return max_x_ - min_x_ + 1; }
  inline int height() const { return max_y_ - min_y_ + 1; }
  inline int area() const { return height() * width(); }
  inline void translate(const Point& offset) {
    min_x_ += offset.x; max_x_ += offset.x;
    min_y_ += offset.y; max_y_ += offset.y;
  }
  inline bool operator==(const Rect& other) const
  { return
      min_y_ == other.min_y_ &&
      min_x_ == other.min_x_ &&
      max_y_ == other.max_y_ &&
      max_x_ == other.max_x_; }
  inline bool operator!=(const Rect& other) const
  { return !operator==(other); }
  inline Rect operator+(const Point& point) const
  { return Rect(this->offset() + point, this->height(), this->width()); }
  inline Rect operator-(const Point& point) const
  { return Rect(this->offset() - point, this->height(), this->width()); }
  inline Rect operator*(const Point& point) const
  { return Rect(this->offset() * point,
      this->height() * point.y, this->width() * point.x); }
  inline Rect operator/(const Point& point) const
  { return Rect(this->offset() / point,
      this->height() / point.y, this->width() / point.x); }
  inline Rect intersection(const Rect& other) const {
    Rect rect;
    rect.min_y_ = std::max(min_y_, other.min_y_);
    rect.min_x_ = std::max(min_x_, other.min_x_);
    rect.max_y_ = std::max(rect.min_y_ - 1, std::min(max_y_, other.max_y_));
    rect.max_x_ = std::max(rect.min_x_ - 1, std::min(max_x_, other.max_x_));
    return rect;
  }

 private:
    int min_y_;
    int min_x_;
    int max_y_;
    int max_x_;
};

std::ostream& operator<<(std::ostream& out, const Rect& rect);

template <typename Sequence>
void DivideLength(size_t total, const Sequence& lengths,
  Sequence& division)
{
  if (lengths.size() == 0) {
    division.push_back(total);
  } else {
    auto iter = lengths.begin();
    if (total == 0 && *iter == 0) {
      division.push_back(0);
    }
    while (total > 0) {
      size_t length = std::min(total, *iter);
      total -= length;
      division.push_back(length);
      if (lengths.end() != iter + 1) {
        ++iter;
      }
    }
  }
}

// Tessellation: nestable Rect iterator class
// Given a bounding Rect, a set of tile heights and widths, and a scan order,
// this iterator will generate a sequence of tiles (Rect) that tessellate the
// bounding Rect. These iterators can be nested, i.e. the bounding Rect can
// itself be a Tessellation (i.e. a Sequence of bounding Rect).
class Tessellation
{
public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = Rect;
  using difference_type = int;
  using pointer = Rect*;
  using reference = Rect&;

  explicit Tessellation(const Rect& bounds,
    std::initializer_list<size_t> heights={1},
    std::initializer_list<size_t> widths={1},
    bool vertical_first=false):
    spec_h_(heights), spec_w_(widths),
    vertical_first_(vertical_first),
    parent_(nullptr),
    bounds_(bounds),
    tile_y_(0),
    tile_x_(0)
  {
    DivideLength(bounds_.height(), spec_h_, tile_h_);
    DivideLength(bounds_.width(), spec_w_, tile_w_);
  }

  template <typename Heights>
  explicit Tessellation(const Rect& bounds,
    const Heights& heights,
    std::initializer_list<size_t> widths,
    bool vertical_first=false):
    spec_h_(heights.cbegin(), heights.cend()),
    spec_w_(widths),
    vertical_first_(vertical_first),
    parent_(nullptr),
    bounds_(bounds),
    tile_y_(0),
    tile_x_(0)
  {
    DivideLength(bounds_.height(), spec_h_, tile_h_);
    DivideLength(bounds_.width(), spec_w_, tile_w_);
  }

  template <typename Heights, typename Widths>
  explicit Tessellation(const Rect& bounds,
    const Heights& heights,
    const Widths& widths,
    bool vertical_first=false):
    spec_h_(heights.cbegin(), heights.cend()),
    spec_w_(widths.cbegin(), widths.cend()),
    vertical_first_(vertical_first),
    parent_(nullptr),
    bounds_(bounds),
    tile_y_(0),
    tile_x_(0)
  {
    DivideLength(bounds_.height(), spec_h_, tile_h_);
    DivideLength(bounds_.width(), spec_w_, tile_w_);
  }

  explicit Tessellation(Tessellation* parent,
    std::initializer_list<size_t> heights={1},
    std::initializer_list<size_t> widths={1},
    bool vertical_first=false):
    spec_h_(heights), spec_w_(widths),
    vertical_first_(vertical_first),
    parent_(parent),
    bounds_(**parent),
    tile_y_(0),
    tile_x_(0)
  {
    DivideLength(bounds_.height(), spec_h_, tile_h_);
    DivideLength(bounds_.width(), spec_w_, tile_w_);
  }

  template <typename Heights, typename Widths>
  explicit Tessellation(Tessellation* parent,
    const Heights& heights,
    const Widths& widths,
    bool vertical_first=false):
    spec_h_(heights.cbegin(), heights.cend()),
    spec_w_(widths.cbegin(), widths.cend()),
    vertical_first_(vertical_first),
    parent_(parent),
    bounds_(**parent),
    tile_y_(0),
    tile_x_(0)
  {
    DivideLength(bounds_.height(), spec_h_, tile_h_);
    DivideLength(bounds_.width(), spec_w_, tile_w_);
  }

  Tessellation child(size_t height, size_t width, bool vertical_first=false)
  {
    return Tessellation(this, {height}, {width}, vertical_first);
  }

  void make_last()
  {
    if (vertical_first_) {
      tile_x_ = tile_w_.size();
      tile_y_ = 0;
    } else {
      tile_y_ = tile_h_.size();
      tile_x_ = 0;
    }
  }

  inline Tessellation& operator++() {
    if (vertical_first_) {
      ++tile_y_;
      if (tile_y_ == tile_h_.size()) {
        tile_y_ = 0;
        ++tile_x_;
      }
    } else {
      ++tile_x_;
      if (tile_x_ == tile_w_.size()) {
        tile_x_ = 0;
        ++tile_y_;
      }
    }
    if (done() && parent_ != nullptr) {
      if (!(++*parent_).done()) {
        bounds_ = **parent_;
        tile_y_ = 0;
        tile_x_ = 0;
        tile_h_.clear();
        tile_w_.clear();
        DivideLength(bounds_.height(), spec_h_, tile_h_);
        DivideLength(bounds_.width(), spec_w_, tile_w_);
      }
    }
    return *this;
  }

  inline Tessellation operator++(int)
    { Tessellation retval = *this; ++(*this); return retval; }

  inline bool operator==(const Tessellation& other) const
  { return done() == other.done(); }

  inline bool operator!=(const Tessellation& other) const
    { return !(*this == other); }

  inline reference operator*() {
    current_ = Rect(
      std::accumulate(tile_h_.begin(), tile_h_.begin() + tile_y_, 0u) +
        bounds_.offset_y(),
      std::accumulate(tile_w_.begin(), tile_w_.begin() + tile_x_, 0u) +
        bounds_.offset_x(),
      tile_h_.at(tile_y_),
      tile_w_.at(tile_x_));
    return current_;
  }

  bool done() const {
    return tile_y_ == tile_h_.size() || tile_x_ == tile_w_.size();
  }

  Tessellation end() const {
    if (parent_ != nullptr) {
      return parent_->end();
    } else {
      auto copy(*this);
      copy.make_last();
      return copy;
    }
  }

private:
  // User specifications for dividing a regions
  const std::vector<size_t> spec_h_;
  const std::vector<size_t> spec_w_;
  // User specified scan order
  const bool vertical_first_;
  // User specified parent iterator (or nullptr)
  Tessellation* parent_;
  // Current rectangle bounds - may be updated from parent
  Rect bounds_;
  // Current tile height and widths
  std::vector<size_t> tile_h_;
  std::vector<size_t> tile_w_;
  // Current tile location
  size_t tile_y_;
  size_t tile_x_;
  // Current tile rectangle
  Rect current_;
};

// Matrix wrapper for raw 2-D pixel-memory buffer
template <typename ValueT, int ORIGIN_Y=0, int ORIGIN_X=0>
class Mat {
public:
  Mat() {}
  Mat(ValueT* buffer, unsigned int height, unsigned int width,
    unsigned int row_width=0, unsigned int interleave=1):
    buffer_(buffer), height_(height), width_(width),
    row_width_(row_width == 0 ? width : row_width),
    interleave_(interleave) {}

  inline ValueT& at(int y, int x) const { assert(buffer_); return buffer_[clamp_x(x - ORIGIN_X)
    + row_width_ * clamp_y(y - ORIGIN_Y)]; }

  inline ValueT& at(const Point& point) const {
    return buffer_[clamp_x(point.x - ORIGIN_X)
    + row_width_ * clamp_y(point.y - ORIGIN_Y)]; }

  inline Rect rect() const { return Rect(ORIGIN_Y, ORIGIN_X, height_, width_); }

  inline void crop(unsigned height, unsigned width) {
    height_ = height;
    width_ = width;
  }

  inline ValueT* data() const { return buffer_; }
  inline void set_data(ValueT* data) { buffer_ = data; }
  inline unsigned int row_width() const { return row_width_; }

  template<int SCALE=0, int PHASES=3>
  class iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ValueT;
    using difference_type = long;
    using pointer = ValueT*;
    using reference = ValueT&;

    iterator(const Mat<ValueT, ORIGIN_Y, ORIGIN_X>& mat,
        int offset_y, int offset_x, int height, int width):
      ptr_(&mat.at(offset_y, offset_x)), done_(height == 0 || width == 0),
      min_(&mat.at(offset_y, offset_x)),
      max_(&mat.at(offset_y + ((height - 1) << SCALE),
        offset_x + ((width - 1) << SCALE)))
    {
      int start = 0;

      if (width > 1) {
        start = count_x_.span(mat.width_, offset_x - ORIGIN_X, width);
      }
      if (height > 1) {
        phase_count& count(width <=1 ? count_x_ : count_y_);
        count.span(mat.height_, offset_y - ORIGIN_Y, height, mat.row_width_,
          start);
      }
    }

    iterator(): ptr_(nullptr), done_(true) {}

    inline iterator& operator++() {
      // Increment
      ptr_ += count_x_.delta_;
      if (count_x_.increment()) {
        ptr_ += count_y_.delta_;
        done_ = count_y_.increment();
      }
      assert (done_ || ptr_ >= min_);
      assert (done_ || ptr_ <= max_);
      return *this;
    }

    inline iterator& operator+=(unsigned plus) {
      while (0 < plus--) {
        ++(*this);
      }
      return *this;
    }

    inline iterator operator++(int)
      { iterator retval = *this; ++(*this); return retval; }
    inline bool operator==(const iterator& other) const
      { return done_ == other.done_; }
    inline bool operator!=(const iterator& other) const
      { return done_ != other.done_; }
    inline ValueT& operator*() const { assert(!done_); return *ptr_; }

  private:
    class phase_count {
    public:
      int span(int width, int offset, int samples, int unit=1, int restart=0) {
        //print("span 2: Width={width}, Offset={offset}, Samples={samples},
        //Scale={scale}".format(width=width, offset=offset, samples=samples,
        //scale=scale))
        int start = 0;
        // If starting before frame data
        if (offset < 0) {
          offset = -offset;
          //print("Offset Count {count}".format(count=offset >> scale))
          int count = std::min(samples, offset >> SCALE);
          append(count, restart);
          samples -= count;
          int remainder = offset & ((1 << SCALE) - 1);
          if (remainder != 0 && samples != 0) {
            offset = (1 << SCALE) - remainder;
            //print("O->W Delta {delta}".format(delta=offset))
            append(1, offset * unit + restart);
            --samples;
          } else {
            offset = 0;
          }
          start = offset;
        }
        // Offset now positive, samples updated
        if (samples == 0) {
          return -start;
        }
        // If last sample within width
        width -= offset;
        if (((samples - 1) << SCALE) < width) {
          start += samples << SCALE;
          //print("Width Count {count}".format(count=samples))
          append(samples, (unit << SCALE) + restart);
        } else {
          // Else last sample outside width
          --width;
          width = std::max(width, 0);
          start += width;
          int count = width >> SCALE;
          samples -= count;
          //print("Width Count {count}".format(count=count))
          append(count, (unit << SCALE) + restart);
          width -= count << SCALE;
          if (0 < width) {
            //print("W->O Delta {delta}".format(delta=width))
            append(1, width * unit + restart);
            --samples;
          }
          //print("Over Count {count}".format(count=samples))
          append(samples, restart);
        }
        //print("Flyback {flyback}".format(flyback=flyback))
        return -start;
      }

      inline void append(int total, int delta=0) {
        if (total == 0) return;
        assert (phases_ < PHASES);
        totals_[phases_] = total;
        deltas_[phases_] = delta;
        if (phases_ == 0) {
          count_ = total;
          delta_ = delta;
        }
        ++phases_;
      }

      inline bool increment() {
        bool done (false);
        if (0 == --count_) {
          done = phases_ <= ++phase_;
          if (done) phase_ = 0;
          assert (phase_ < PHASES);
          count_ = totals_[phase_];
          delta_ = deltas_[phase_];
        }
        return done;
      }

      int delta_;

    private:
      std::array<int, PHASES> totals_;
      std::array<int, PHASES> deltas_;
      int count_ {1};
      int phases_ {0};
      int phase_ {0};
    };

    ValueT* ptr_;
    phase_count count_y_;
    phase_count count_x_;
    bool done_;
    // For DEBUG ONLY
    ValueT* min_;
    ValueT* max_;
  };

  template<int SCALE=0, int PHASES=3>
  inline iterator<SCALE, PHASES> begin(int offset_y, int offset_x,
    int height, int width) const
    { return iterator<SCALE, PHASES>(*this, offset_y,
      offset_x, height, width); }
  template<int SCALE=0, int PHASES=3>
  inline iterator<SCALE, PHASES> begin(const Rect& rect) const
    { return iterator<SCALE, PHASES>(*this, rect.offset_y(),
      rect.offset_x(), rect.height(), rect.width()); }
  template<int SCALE=0, int PHASES=3>
  inline iterator<SCALE, PHASES> begin() const
    { return begin(rect()); }
  template<int SCALE=0, int PHASES=3>
  inline iterator<SCALE, PHASES> end() const
    { return iterator<SCALE, PHASES>(); }
  friend std::ostream& operator<<(std::ostream& str,
    const Mat<ValueT, ORIGIN_Y, ORIGIN_X>& mat) {
      str << std::hex << std::setfill('0');
      for(int row(0); row < mat.height_; ++row) {
        for(int col(0); col < mat.width_; ++col) {
          str << " " << std::setw(2*sizeof(ValueT)) <<
            static_cast<int>(mat.buffer_[col + row * mat.row_width_]);
        }
        str << std::endl;
      }
      return str << std::dec;
    }

private:
  ValueT* buffer_ { nullptr };
  unsigned int height_ { 0 };
  unsigned int width_ { 0 };
  unsigned int row_width_ { 0 };
  unsigned int interleave_ { 1 };

  inline int clamp_x(int x) const
    { return x < 0 ? x % interleave_ :
             x >= (int)width_ ? width_ - interleave_ + x % interleave_ : x; }
  inline int clamp_y(int y) const
    { return y < 0 ? 0 : y >= (int)height_ ? height_-1 : y; }
};
