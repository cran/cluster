#!/usr/bin/perl
foreach $f (<print.*>) {
  open(IN, "< $f");
  open(OUT, "> $f.tmp");
  select(OUT);
  while (<IN>) {
    s/\.PP/.DT/;
    s/^$/.DT/;
    print;
  }
  close(IN);
  close(OUT);
  unlink($f);
  rename($f . ".tmp", $f);
}
