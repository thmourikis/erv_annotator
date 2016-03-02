#
# Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

# testing the Ruby bindings for AnnotationSketch (similar to the sketch tool)

require 'gtruby'

if ARGV.size != 2 then
  STDERR.puts "Usage: #{$0} PNG_file GFF3_file"
  STDERR.puts "Create PNG representation of GFF3 annotation file."
  exit(1)
end

pngfile  = ARGV[0]
gff3file = ARGV[1]

in_stream = GT::GFF3InStream.new(gff3file)
feature_index = GT::FeatureIndexMemory.new()
feature_stream = GT::FeatureStream.new(in_stream, feature_index)
gn = feature_stream.next_tree()
# fill feature index
while (gn) do
  gn = feature_stream.next_tree()
end

seqid = feature_index.get_first_seqid()
range = feature_index.get_range_for_seqid(seqid)

style = GT::Style.new()
diagram = GT::Diagram.from_index(feature_index, seqid, range, style)
layout = GT::Layout.new(diagram, 700, style)
ii = GT::ImageInfo.new()
canvas = GT::CanvasCairoFile.new(style, 700, layout.get_height, ii)
layout.sketch(canvas)
canvas.to_file(pngfile)
