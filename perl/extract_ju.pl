#! /usr/bin/perl -w
# ============================================================

# Extract class_code j and u from GTF

# For class_code = j
#   - change gene_id to FBgn, FBgn is identified using gene_name
#   - keep TCONS as transcript_id
#   - remove tss_id

# For class_code = u
#   - keep XLOC and TCONS IDs
#   - remove tss_id

$gtf = $ARGV[0];

# read gtf file
open GTF, "<$gtf";

while (<GTF>) {
  chomp $_;
  @line = split /\t/, $_;
  $txinfo = $line[8];
  # remove tss_id and p_id no matter what
  $txinfo =~ s/\s+tss_id \".*?\";//;
  $txinfo =~ s/\s+p_id \".*?\";//;

  # IF class_code is u;
  if ($txinfo =~ m/class_code \"u\";/) {
    $line[8] = $txinfo;
    print join("\t", @line), "\n";
  }
  # IF class_code  is j;
  if ($txinfo =~ m/class_code \"j\";/) {
    if ($txinfo =~ m/gene_name \"(.*?)\";/) {
      $gene_id = $1;
    }
    $txinfo =~ s/gene_id \".*?\";/gene_id \"$gene_id\";/;
    $line[8] = $txinfo;
    print join("\t", @line), "\n";
  }
}

close GTF;
    
