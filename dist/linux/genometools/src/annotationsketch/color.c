/*
  Copyright (c) 2007 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/ma.h"
#include "annotationsketch/color_api.h"

GtColor* gt_color_new(double red, double green, double blue, double alpha)
{
  GtColor *color = gt_malloc(sizeof *color);
  color->red = red;
  color->green = green;
  color->blue = blue;
  color->alpha = alpha;
  return color;
}

void gt_color_set(GtColor *color, double red, double green, double blue,
                  double alpha)
{
  gt_assert(color);
  color->red = red;
  color->green = green;
  color->blue = blue;
  color->alpha = alpha;
}

bool gt_color_equals(const GtColor *c1, const GtColor *c2)
{
  gt_assert(c1 && c2);
  return ((c1->red == c2->red) && (c1->green == c2->green) &&
          (c1->blue == c2->blue) && (c1->alpha == c2->alpha));
}

void gt_color_delete(GtColor *color)
{
  if (!color) return;
  gt_free(color);
}
