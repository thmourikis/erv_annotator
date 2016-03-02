/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef REC_MAP_H
#define REC_MAP_H

#include "annotationsketch/rec_map_api.h"

GtRecMap* gt_rec_map_new(double nw_x, double nw_y, double se_x, double se_y,
                         GtFeatureNode*);
GtRecMap* gt_rec_map_ref(GtRecMap *rm);
void      gt_rec_map_set_omitted_children(GtRecMap *rm, bool status);
int       gt_rec_map_format_html_imagemap_coords(const GtRecMap*, char*,
                                                 size_t);
void      gt_rec_map_delete(GtRecMap*);

#endif
