/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef DYNALLOC_H
#define DYNALLOC_H

#include <stdlib.h>

#ifndef SIZE_MAX
#define SIZE_MAX (~(size_t)0)
#endif

/*
  Do not use this function directly! It is only used internally to implement
  dynamic arrays, dynamic strings, and the like.

  It dynamically re-/allocates memory. Thereby, usually more memory
  is allocated than what was asked for (to avoid frequent realloc calls).

  <ptr> the previously allocated memory
  <allocated> allocated memory size before and after the call
  <size> requested memory size
*/

void* gt_dynalloc(void *ptr, size_t *allocated, size_t size);

#endif
