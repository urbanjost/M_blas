#!/bin/bash
(
exec 2>&1
fpm test cblat1 
#fpm test cblat1_a
fpm test cblat2    < test/cblat2.in
fpm test cblat2_a  < test/cblat2.in
fpm test cblat3    < test/cblat3.in

fpm test dblat1 
fpm test dblat1_a
fpm test dblat2    < test/dblat2.in
fpm test dblat3    < test/dblat3.in
fpm test dblat3_a  < test/dblat3.in

fpm test sblat1 
fpm test sblat2    < test/sblat2.in
fpm test sblat3    < test/sblat3.in

fpm test zblat1  
fpm test zblat2    < test/zblat2.in
fpm test zblat3    < test/zblat3.in

cat cblat2.out
cat cblat3.out
cat dblat2.out
cat dblat3.out
cat sblat2.out
cat sblat3.out
cat zblat2.out
cat zblat3.out

rm -f cblat2.out cblat3.out dblat2.out dblat3.out sblat2.out sblat3.out zblat2.out zblat3.out

)|cat -n|tee run.log

grep -i fail run.log
exit
