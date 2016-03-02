/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef DICE_HMM_H
#define DICE_HMM_H

#include "core/alphabet.h"
#include "extended/hmm.h"

typedef enum {
  DICE_FAIR,
  DICE_LOADED,
  DICE_NUM_OF_STATES
} Dice_states;

typedef enum {
  ONE,
  TWO,
  THREE,
  FOUR,
  FIVE,
  SIX,
  DICE_NUM_OF_SYMBOLS
} Dice_emissions;

GtHMM*      gt_dice_hmm_loaded(void);
GtHMM*      gt_dice_hmm_fair(void);
GtAlphabet* gt_dice_hmm_alphabet(void);

#endif
