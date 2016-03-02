#
# Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
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

# testing the Ruby bindings for the database-backed feature stores

require 'gtruby'

if ARGV.size != 1 then
  STDERR.puts "Usage: #{$0} GFF3_db"
  STDERR.puts "Test the FeatureIndex bindings on GFF-like db."
  exit(1)
end

gffdb = ARGV[0]

rdb = GT::RDBSqlite.new(gffdb)
adb = GT::AnnoDBGFFlike.new()
feature_index = adb.get_feature_index(rdb)

seqid = feature_index.get_first_seqid()
features = feature_index.get_features_for_seqid(seqid)
raise if not features

rng = feature_index.get_range_for_seqid(seqid)
features_r = feature_index.get_features_for_range(rng.start, rng.end, seqid)
raise unless features.length == features_r.length

gff3_visitor = GT::GFF3Visitor.new()

features.each do |feature|
  feature.accept(gff3_visitor)
end
