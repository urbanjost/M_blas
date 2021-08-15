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
)
exit
