#ifndef __FBVID_SHD_H__
#define __FBVID_SHD_H__

#include <vector>
#include <array>
#include "mat.h"

#define HIST_SIZE  256
#define DHIST_SIZE 512

typedef enum _HIST_TYPE {
  HIST_SRC  = 0,
  HIST_REF  = 1,
  HIST_DIFF = 2
}HIST_TYPE;

class vid_shd {
 public:
  vid_shd(unsigned int height, unsigned int width, unsigned int bit_depth = 0,
    unsigned int chroma = 0, unsigned int chroma_format = 0,
    unsigned int start_frame = 0):
    height_(height), width_(width), bit_depth_(bit_depth),
    is_chroma_(chroma), chroma_format_(chroma_format) {
      init_hist();
      start_frame += 0;
    }
  virtual ~vid_shd() {}
  bool set_buffers(unsigned int* cur, unsigned int* ref);
  virtual bool calc_core();
  vid_shd (const vid_shd&) = delete;
  vid_shd& operator=(const vid_shd&) = delete;

  void set_params(const int height, const int width, const int is_chroma) {
    height_    = height;
    width_     = width;
    is_chroma_ = is_chroma;
  }
  const std::vector<uint32_t>& get_hist(const HIST_TYPE type, const bool chroma);

 protected:
  unsigned int height_;
  unsigned int width_;
  unsigned int bit_depth_;
  unsigned int is_chroma_;
  unsigned int chroma_format_;
  Mat<const unsigned int> src_;
  Mat<const unsigned int> ref_;
  std::array<std::vector<uint32_t>, 2> src_hist_;
  std::array<std::vector<uint32_t>, 2> ref_hist_;
  std::array<std::vector<uint32_t>, 2> diff_hist_;
  std::array<uint32_t, 2> sadval_;
  std::array<bool, 2> buffer_set_;
  std::vector<uint64_t> shd_out_;

  void init_hist();
  bool calc_core_luma();
  bool calc_core_chroma();
  void calc_out_data();
};
#endif
