--[[
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
]]

-- testing the Lua bindings for FeatureIndex and FeatureStream classes

function usage()
io.stderr:write(string.format("Usage: %s testdata_dir\n", arg[0]))
  io.stderr:write("Test the FeatureIndex and FeatureStream bindings.\n")
  os.exit(1)
end


if #arg == 1 then
  testdata = arg[1]
else
  usage()
end

-- set up the feature stream
genome_stream = gt.gff3_in_stream_new_sorted(testdata.."/gff3_file_1_short.txt")
feature_index = gt.feature_index_memory_new()
genome_stream = gt.feature_stream_new(genome_stream, feature_index)
collectgarbage()

feature = genome_stream:next_tree()
while (feature) do
  feature = genome_stream:next_tree()
end

features = feature_index:get_features_for_seqid("ctg123")
assert(features)
gff3_visitor = gt.gff3_visitor_new()

for i,feature in ipairs(features) do
  feature:accept(gff3_visitor)
end

range = gt.range_new(1, 100)

-- more tests
fi    = gt.feature_index_memory_new()
sr    = gt.region_node_new("chr1", 1, 100)
gf    = gt.feature_node_new("chr1", "gene", 1, 100 , "+")
rval, err = pcall(GenomeTools_feature_index.add_feature_node, fi, nil)
assert(not rval)
assert(string.find(err, "genome_node expected"))
rval, err = pcall(GenomeTools_feature_index.add_region_node, fi, nil)
assert(not rval)
assert(string.find(err, "genome_node expected"))
rval, err = pcall(GenomeTools_feature_index.add_region_node, fi, gf)
assert(not rval)
assert(string.find(err, "not a region node"))
fi:add_region_node(sr)
fi:add_feature_node(gf)
