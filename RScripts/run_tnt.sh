 #!/bin/bash
N=$1
(
for j in $(seq 0 $2); do
((i=i%N)); ((i++==0)) && wait
/Users/grahamslater/tnt-mac-no-tax-limit/tnt.command mxram 1024, log tnt_log.$j.txt, run mymatrix.tnt, echo= , timeout 24:00:00, rseed0, rseed*,hold 1000,xmult= level 1, taxname=, tsave *trees_tnt.$j.tnt, save, tsave / , scores, log / , quit &
done
)