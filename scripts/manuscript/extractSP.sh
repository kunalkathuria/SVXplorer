#!/bin.bash
cat $1.out.sim9 | awk 'NR % 6 == 1' > $1.dels.s
cat $1.out.sim9 | awk 'NR % 6 == 2' > $1.dels.p
cat $1.out.sim9 | awk 'NR % 6 == 3' > $1.dups.s
cat $1.out.sim9 | awk 'NR % 6 == 4' > $1.dups.p
cat $1.out.sim9 | awk 'NR % 6 == 5' > $1.invs.s
cat $1.out.sim9 | awk 'NR % 6 == 0' > $1.invs.p
