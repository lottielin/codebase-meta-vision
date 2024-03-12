#include <iostream>
#include <fstream>
#include <sstream>
#include "vid_shd.h"

using namespace std;

void vid_shd::init_hist() {
  src_hist_ [0].resize(HIST_SIZE,  0);
  ref_hist_ [0].resize(HIST_SIZE,  0);
  diff_hist_[0].resize(DHIST_SIZE, 0);
  src_hist_ [1].resize(HIST_SIZE,  0);
  ref_hist_ [1].resize(HIST_SIZE,  0);
  diff_hist_[1].resize(DHIST_SIZE, 0);
}

const vector<uint32_t>& vid_shd::get_hist(
  const HIST_TYPE type, const bool chroma) {
  switch(type) {
    case HIST_SRC:
      return (chroma)?src_hist_[1]:src_hist_[0];
    case HIST_REF:
      return (chroma)?ref_hist_[1]:ref_hist_[0];
    case HIST_DIFF:
      return (chroma)?diff_hist_[1]:diff_hist_[0];
    default:
      return (chroma)?src_hist_[1]:src_hist_[0];
  }
}

// set buffer and convert to mat
bool vid_shd::set_buffers(unsigned int* cur, unsigned int* ref) {

  unsigned int chroma_mult = (is_chroma_)?2:1;
  buffer_set_[0] = buffer_set_[1] = false;
  if (cur){
    src_ = Mat<const unsigned int>(cur, height_, chroma_mult*width_);
    buffer_set_[0] = true;
  }
  if (ref){
    ref_ = Mat<const unsigned int>(ref, height_, chroma_mult*width_);
    buffer_set_[1] = true;
  }

  return (buffer_set_[0] & buffer_set_[1]);
}

bool vid_shd::calc_core_luma() {

  unsigned int pixel0, pixel1;
  unsigned int pixel0_8, pixel1_8;
  const unsigned int num_pixels = height_*width_;
  int pixel;

  // initialize hist and sad
  fill(src_hist_[0].begin(), src_hist_[0].end(), 0);
  fill(ref_hist_[0].begin(), ref_hist_[0].end(), 0);
  fill(diff_hist_[0].begin(), diff_hist_[0].end(), 0);
  sadval_[0] = 0;

  Rect src_region(0, 0, height_, width_);
  auto src_region_iter = src_region.begin();
  Rect ref_region(0, 0, height_, width_);
  auto ref_region_iter = ref_region.begin();

  for (unsigned int ind = 0; ind < num_pixels; ++ind) {
    pixel0 = src_.at(*src_region_iter);
    pixel0_8 = (bit_depth_ == 10)?((pixel0+2)>>2):pixel0;
    pixel0_8 = (pixel0_8 > 255)?255:pixel0_8;
    src_hist_[0][pixel0_8]++;
    ++src_region_iter;

    pixel1 = ref_.at(*ref_region_iter);
    pixel1_8 = (bit_depth_ == 10)?((pixel1+2)>>2):pixel1;
    pixel1_8 = (pixel1_8 > 255)?255:pixel1_8;
    ref_hist_[0][pixel1_8]++;
    ++ref_region_iter;

    pixel = pixel0_8 - pixel1_8;
    diff_hist_[0][pixel+DHIST_SIZE/2]++;

    pixel = (pixel < 0)?-pixel:pixel;
    sadval_[0] += pixel;
  }

  return true;
}

bool vid_shd::calc_core_chroma()
{
  unsigned int pixel0, pixel1;
  unsigned int pixel0_8, pixel1_8;
  const unsigned int num_pixels = height_*width_;
  int pixel;

  // initialize histograms and sadval
  for (int idx=0; idx < 2; ++idx) {
    fill(src_hist_[idx].begin(), src_hist_[idx].end(), 0);
    fill(ref_hist_[idx].begin(), ref_hist_[idx].end(), 0);
    fill(diff_hist_[idx].begin(), diff_hist_[idx].end(), 0);
    sadval_[idx] = 0;
  }

  Rect src_region(0, 0, height_, 2*width_);
  auto src_region_iter = src_region.begin();
  Rect ref_region(0, 0, height_, 2*width_);
  auto ref_region_iter = ref_region.begin();

  for (unsigned int ind = 0; ind < (2*num_pixels); ind+=2) {
    // chroma_b
    pixel0 = src_.at(*src_region_iter);
    pixel0_8 = (bit_depth_ == 10)?((pixel0+2)>>2):pixel0;
    pixel0_8 = (pixel0_8 > 255)?255:pixel0_8;
    src_hist_[0][pixel0_8]++;
    ++src_region_iter;

    pixel1 = ref_.at(*ref_region_iter);
    pixel1_8 = (bit_depth_ == 10)?((pixel1+2)>>2):pixel1;
    pixel1_8 = (pixel1_8 > 255)?255:pixel1_8;
    ref_hist_[0][pixel1_8]++;
    ++ref_region_iter;

    pixel = pixel0_8 - pixel1_8;
    diff_hist_[0][pixel+DHIST_SIZE/2]++;

    pixel = (pixel < 0)?-pixel:pixel;
    sadval_[0] += pixel;

    // chroma_r
    pixel0 = src_.at(*src_region_iter);
    pixel0_8 = (bit_depth_ == 10)?((pixel0+2)>>2):pixel0;
    pixel0_8 = (pixel0_8 > 255)?255:pixel0_8;
    src_hist_[1][pixel0_8]++;
    ++src_region_iter;

    pixel1 = ref_.at(*ref_region_iter);
    pixel1_8 = (bit_depth_ == 10)?((pixel1+2)>>2):pixel1;
    pixel1_8 = (pixel1_8 > 255)?255:pixel1_8;
    ref_hist_[1][pixel1_8]++;
    ++ref_region_iter;

    pixel = pixel0_8 - pixel1_8;
    diff_hist_[1][pixel+DHIST_SIZE/2]++;

    pixel = (pixel < 0)?-pixel:pixel;
    sadval_[1] += pixel;
  }

  return true;
}

void vid_shd::calc_out_data() {

  unsigned int out_entries = HIST_SIZE * 2 + DHIST_SIZE;
  out_entries = (is_chroma_) ? 2 * out_entries : out_entries;
  unsigned int out_64bitwords = out_entries / 2;
  unsigned int out_size = (out_64bitwords % 2 == 1) ?
        out_64bitwords + 1 : out_64bitwords;

  shd_out_.resize(out_size, 0);
  unsigned int shd_out_index = 0;
  uint64_t word1, word2;

  // src histograms
  for (unsigned int ind = 0; ind < HIST_SIZE; ind += 2) {
    word1 = src_hist_[0][ind];
    word2 = src_hist_[0][ind + 1];
    shd_out_[shd_out_index] = (word2 << 32) + word1;
    shd_out_index++;
  }

  if (is_chroma_) {
    for (unsigned int ind = 0; ind < HIST_SIZE; ind += 2) {
      word1 = src_hist_[1][ind];
      word2 = src_hist_[1][ind + 1];
      shd_out_[shd_out_index] = (word2 << 32) + word1;
      shd_out_index++;
    }
  }

  // ref histograms
  for (unsigned int ind = 0; ind < HIST_SIZE; ind += 2) {
    word1 = ref_hist_[0][ind];
    word2 = ref_hist_[0][ind + 1];
    shd_out_[shd_out_index] = (word2 << 32) + word1;
    shd_out_index++;
  }

  if (is_chroma_) {
    for (unsigned int ind = 0; ind < HIST_SIZE; ind += 2) {
      word1 = ref_hist_[1][ind];
      word2 = ref_hist_[1][ind + 1];
      shd_out_[shd_out_index] = (word2 << 32) + word1;
      shd_out_index++;
    }
  }

  // diff histograms + sad in first position
  diff_hist_[0][0] = sadval_[0];
  for (unsigned int ind = 0; ind < DHIST_SIZE; ind += 2) {
    word1 = diff_hist_[0][ind];
    word2 = diff_hist_[0][ind + 1];
    shd_out_[shd_out_index] = (word2 << 32) + word1;
    shd_out_index++;
  }

  if (is_chroma_) {
    diff_hist_[1][0] = sadval_[1];
    for (unsigned int ind = 0; ind < DHIST_SIZE; ind += 2) {
      word1 = diff_hist_[1][ind];
      word2 = diff_hist_[1][ind + 1];
      shd_out_[shd_out_index] = (word2 << 32) + word1;
      shd_out_index++;
    }
  }
}

bool vid_shd::calc_core()
{
  if (!buffer_set_[0] | !buffer_set_[1]) {
    cout << "exited due to not both buffers set!" << endl;
    return false;
  }

  if (is_chroma_) {
    calc_core_chroma();
  } else {
    calc_core_luma();
  }

  calc_out_data();
  return true;
}
