/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef TRANSCRIPT_USED_EXONS_H
#define TRANSCRIPT_USED_EXONS_H

#include "core/dlist.h"

typedef struct GtTranscriptUsedExons GtTranscriptUsedExons;

GtTranscriptUsedExons* gt_transcript_used_exons_new(void);
GtDlist*          gt_transcript_used_exons_get_all(GtTranscriptUsedExons*);
GtDlist*          gt_transcript_used_exons_get_single(GtTranscriptUsedExons*);
GtDlist*          gt_transcript_used_exons_get_initial(GtTranscriptUsedExons*);
GtDlist*          gt_transcript_used_exons_get_internal(GtTranscriptUsedExons*);
GtDlist*          gt_transcript_used_exons_get_terminal(GtTranscriptUsedExons*);
void              gt_transcript_used_exons_delete(GtTranscriptUsedExons*);

#endif
