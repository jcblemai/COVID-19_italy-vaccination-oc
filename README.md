Old version contains:
- CVO_implicit: implicit force of infection: 1st order approximation of the force of infection to mitigate
the fact that the relationship between nodes is not taken into account
- Between ocp_trial and ocp_trial-vector (diff: in vector: loop over time then over nodes, and compute foi once for all nodes. Moreover use matrix multiplication over loop (which I later found that it produced non-structural zero, so this change was undone) )

Notebooks:

- ocp_trial_vector: my starting version
- ocp_trial_vector-crtl: control out of the loop, but no structural zero (so matrix formula)
- ocp_trial_vector-crtl-struct0: same with strutural zero. Most refined version, the final program 
is based on this.
- ocp_trial-vector3-JIT: just in time trial. Not really happy yet but might be worth considering.

Segfault with gdb:
```
*** Error in `/home/chadi/miniconda3/envs/ocp-covid/bin/python': corrupted size vs. prev_size: 0x00007fff1811d4b0 ***
======= Backtrace: =========
/lib/x86_64-linux-gnu/libc.so.6(+0x777f5)[0x7ffff78677f5]
/lib/x86_64-linux-gnu/libc.so.6(+0x7e9ec)[0x7ffff786e9ec]
/lib/x86_64-linux-gnu/libc.so.6(+0x81d0c)[0x7ffff7871d0c]
/lib/x86_64-linux-gnu/libc.so.6(__libc_malloc+0x54)[0x7ffff78741d4]
/home/chadi/src/hsl/20190503/lib/libhsl.so(+0xe295f)[0x7fffe4bfd95f]
/home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/../../../././././libgomp.so.1(+0x146d5)[0x7ffff1a886d5]
/lib/x86_64-linux-gnu/libpthread.so.0(+0x76ba)[0x7ffff7bc16ba]
/lib/x86_64-linux-gnu/libc.so.6(clone+0x6d)[0x7ffff78f74dd]
======= Memory map: ========
555555554000-5555555af000 r--p 00000000 08:06 44182157                   /home/chadi/miniconda3/envs/ocp-covid/bin/python3.7
5555555af000-55555578d000 r-xp 0005b000 08:06 44182157                   /home/chadi/miniconda3/envs/ocp-covid/bin/python3.7
55555578d000-555555835000 r--p 00239000 08:06 44182157                   /home/chadi/miniconda3/envs/ocp-covid/bin/python3.7
555555835000-555555838000 r--p 002e0000 08:06 44182157                   /home/chadi/miniconda3/envs/ocp-covid/bin/python3.7
555555838000-5555558a1000 rw-p 002e3000 08:06 44182157                   /home/chadi/miniconda3/envs/ocp-covid/bin/python3.7
5555558a1000-55558136c000 rw-p 00000000 00:00 0                          [heap]
7ffef8000000-7ffef8021000 rw-p 00000000 00:00 0
7ffef8021000-7ffefc000000 ---p 00000000 00:00 0
7fff00000000-7fff0012f000 rw-p 00000000 00:00 0
7fff0012f000-7fff04000000 ---p 00000000 00:00 0
7fff08000000-7fff0812c000 rw-p 00000000 00:00 0
7fff0812c000-7fff0c000000 ---p 00000000 00:00 0
7fff10000000-7fff1012d000 rw-p 00000000 00:00 0
7fff1012d000-7fff14000000 ---p 00000000 00:00 0
7fff18000000-7fff1812c000 rw-p 00000000 00:00 0
7fff1812c000-7fff1c000000 ---p 00000000 00:00 0
7fff20000000-7fff2012d000 rw-p 00000000 00:00 0
7fff2012d000-7fff24000000 ---p 00000000 00:00 0
7fff28000000-7fff28128000 rw-p 00000000 00:00 0
7fff28128000-7fff2c000000 ---p 00000000 00:00 0
7fff2c000000-7fff2c12d000 rw-p 00000000 00:00 0
7fff2c12d000-7fff30000000 ---p 00000000 00:00 0
7fff30000000-7fff30129000 rw-p 00000000 00:00 0
7fff30129000-7fff34000000 ---p 00000000 00:00 0
7fff38000000-7fff3812d000 rw-p 00000000 00:00 0
7fff3812d000-7fff3c000000 ---p 00000000 00:00 0
7fff3e000000-7fff40000000 rw-p 00000000 00:00 0
7fff40000000-7fff4012d000 rw-p 00000000 00:00 0
7fff4012d000-7fff44000000 ---p 00000000 00:00 0
7fff44000000-7fff48000000 rw-p 00000000 00:00 0
7fff48000000-7fff4812d000 rw-p 00000000 00:00 0
7fff4812d000-7fff4c000000 ---p 00000000 00:00 0
7fff4c000000-7fff4c12e000 rw-p 00000000 00:00 0
7fff4c12e000-7fff50000000 ---p 00000000 00:00 0
7fff50000000-7fff5012f000 rw-p 00000000 00:00 0
7fff5012f000-7fff54000000 ---p 00000000 00:00 0
7fff54000000-7fff58000000 rw-p 00000000 00:00 0
7fff58000000-7fff5812d000 rw-p 00000000 00:00 0
7fff5812d000-7fff5c000000 ---p 00000000 00:00 0
7fff5c000000-7fff60000000 rw-p 00000000 00:00 0
7fff60000000-7fff60126000 rw-p 00000000 00:00 0
7fff60126000-7fff64000000 ---p 00000000 00:00 0
7fff64000000-7fff6413e000 rw-p 00000000 00:00 0
7fff6413e000-7fff68000000 ---p 00000000 00:00 0
7fff68000000-7fff6812d000 rw-p 00000000 00:00 0
7fff6812d000-7fff6c000000 ---p 00000000 00:00 0
7fff6c000000-7fff6c12c000 rw-p 00000000 00:00 0
7fff6c12c000-7fff70000000 ---p 00000000 00:00 0
7fff70000000-7fff70129000 rw-p 00000000 00:00 0
7fff70129000-7fff74000000 ---p 00000000 00:00 0
7fff74000000-7fff7412d000 rw-p 00000000 00:00 0
7fff7412d000-7fff78000000 ---p 00000000 00:00 0
7fff78000000-7fff78124000 rw-p 00000000 00:00 0
7fff78124000-7fff7c000000 ---p 00000000 00:00 0
7fff7c000000-7fff7c11e000 rw-p 00000000 00:00 0
7fff7c11e000-7fff80000000 ---p 00000000 00:00 0
7fff80000000-7fff80128000 rw-p 00000000 00:00 0
7fff80128000-7fff84000000 ---p 00000000 00:00 0
7fff8465e000-7fff8665e000 rw-p 00000000 00:00 0
7fff8665e000-7fff8665f000 ---p 00000000 00:00 0
7fff8665f000-7fff86e5f000 rw-p 00000000 00:00 0
7fff86e5f000-7fff86e60000 ---p 00000000 00:00 0
7fff86e60000-7fff87660000 rw-p 00000000 00:00 0
7fff87660000-7fff87661000 ---p 00000000 00:00 0
7fff87661000-7fff87e61000 rw-p 00000000 00:00 0
7fff87e61000-7fff87e62000 ---p 00000000 00:00 0
7fff87e62000-7fff88662000 rw-p 00000000 00:00 0
7fff88662000-7fff88663000 ---p 00000000 00:00 0
7fff88663000-7fff88e63000 rw-p 00000000 00:00 0
7fff88e63000-7fff88e64000 ---p 00000000 00:00 0
7fff88e64000-7fff89664000 rw-p 00000000 00:00 0
7fff89664000-7fff89665000 ---p 00000000 00:00 0
7fff89665000-7fff89e65000 rw-p 00000000 00:00 0
7fff89e65000-7fff89e66000 ---p 00000000 00:00 0
7fff89e66000-7fff8a666000 rw-p 00000000 00:00 0
7fff8a666000-7fff8a667000 ---p 00000000 00:00 0
7fff8a667000-7fff8ae67000 rw-p 00000000 00:00 0
7fff8ae67000-7fff8ae68000 ---p 00000000 00:00 0
7fff8ae68000-7fff8b668000 rw-p 00000000 00:00 0
7fff8b668000-7fff8b669000 ---p 00000000 00:00 0
7fff8b669000-7fff8be69000 rw-p 00000000 00:00 0
7fff8be69000-7fff8be6a000 ---p 00000000 00:00 0
7fff8be6a000-7fff8c66a000 rw-p 00000000 00:00 0
7fff8c66a000-7fff8c66b000 ---p 00000000 00:00 0
7fff8c66b000-7fff8ce6b000 rw-p 00000000 00:00 0
7fff8ce6b000-7fff8ce6c000 ---p 00000000 00:00 0
7fff8ce6c000-7fff8d66c000 rw-p 00000000 00:00 0
7fff8d66c000-7fff8d66d000 ---p 00000000 00:00 0
7fff8d66d000-7fff8de6d000 rw-p 00000000 00:00 0
7fff8de6d000-7fff8de6e000 ---p 00000000 00:00 0
7fff8de6e000-7fff8e66e000 rw-p 00000000 00:00 0
7fff8e66e000-7fff8e66f000 ---p 00000000 00:00 0
7fff8e66f000-7fff8ee6f000 rw-p 00000000 00:00 0
7fff8ee6f000-7fff8ee70000 ---p 00000000 00:00 0
7fff8ee70000-7fff8f670000 rw-p 00000000 00:00 0
7fff8f670000-7fff8f671000 ---p 00000000 00:00 0
7fff8f671000-7fff8fe71000 rw-p 00000000 00:00 0
7fff8fe71000-7fff8fe72000 ---p 00000000 00:00 0
7fff8fe72000-7fff90672000 rw-p 00000000 00:00 0
7fff90672000-7fff90673000 ---p 00000000 00:00 0
7fff90673000-7fff90e73000 rw-p 00000000 00:00 0
7fff90e73000-7fff90e74000 ---p 00000000 00:00 0
7fff90e74000-7fff91674000 rw-p 00000000 00:00 0
7fff91674000-7fff91675000 ---p 00000000 00:00 0
7fff91675000-7fff93e75000 rw-p 00000000 00:00 0
7fff93e75000-7fff93e76000 ---p 00000000 00:00 0
7fff93e76000-7fff94676000 rw-p 00000000 00:00 0
7fff94676000-7fff94677000 ---p 00000000 00:00 0
7fff94677000-7fff94e77000 rw-p 00000000 00:00 0
7fff94e77000-7fff94e78000 ---p 00000000 00:00 0
7fff94e78000-7fff95678000 rw-p 00000000 00:00 0
7fff95678000-7fff95679000 ---p 00000000 00:00 0
7fff95679000-7fff95e79000 rw-p 00000000 00:00 0
7fff95e79000-7fff95e7a000 ---p 00000000 00:00 0
7fff95e7a000-7fff9667a000 rw-p 00000000 00:00 0
7fff9667a000-7fff9667b000 ---p 00000000 00:00 0
7fff9667b000-7fff96e7b000 rw-p 00000000 00:00 0
7fff96e7b000-7fff96e7c000 ---p 00000000 00:00 0
7fff96e7c000-7fff9767c000 rw-p 00000000 00:00 0
7fff9767c000-7fff9767d000 ---p 00000000 00:00 0
7fff9767d000-7fff97e7d000 rw-p 00000000 00:00 0
7fff97e7d000-7fff97e7e000 ---p 00000000 00:00 0
7fff97e7e000-7fff9867e000 rw-p 00000000 00:00 0
7fff9867e000-7fff9867f000 ---p 00000000 00:00 0
7fff9867f000-7fff98e7f000 rw-p 00000000 00:00 0
7fff98e7f000-7fff98e80000 ---p 00000000 00:00 0
7fff98e80000-7fff99680000 rw-p 00000000 00:00 0
7fff99680000-7fff99681000 ---p 00000000 00:00 0
7fff99681000-7fff99e81000 rw-p 00000000 00:00 0
7fff99e81000-7fff99e82000 ---p 00000000 00:00 0
7fff99e82000-7fffa17da000 rw-p 00000000 00:00 0
7fffa1fd7000-7fffa8000000 rw-p 00000000 00:00 0
7fffa8000000-7fffa8021000 rw-p 00000000 00:00 0
7fffa8021000-7fffac000000 ---p 00000000 00:00 0
7fffac46b000-7fffac46c000 ---p 00000000 00:00 0
7fffac46c000-7fffacc6c000 rw-p 00000000 00:00 0
7fffacc6c000-7fffacd57000 r--p 00000000 08:06 44304412                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.10.so
7fffacd57000-7fffacd58000 r-xp 000eb000 08:06 44304412                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.10.so
7fffacd58000-7fffacd5c000 ---p 000ec000 08:06 44304412                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.10.so
7fffacd5c000-7fffae6d1000 r-xp 000f0000 08:06 44304412                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.10.so
7fffae6d1000-7fffae897000 r--p 01a65000 08:06 44304412                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.10.so
7fffae897000-7fffae8a5000 r--p 01c2a000 08:06 44304412                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.10.so
7fffae8a5000-7fffae8b2000 rw-p 01c38000 08:06 44304412                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.10.so
7fffae8b2000-7fffae943000 rw-p 00000000 00:00 0
7fffae943000-7fffae947000 r--p 00000000 08:06 44042077                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_crs.cpython-37m-x86_64-linux-gnu.so
7fffae947000-7fffae952000 r-xp 00004000 08:06 44042077                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_crs.cpython-37m-x86_64-linux-gnu.so
7fffae952000-7fffae954000 r--p 0000f000 08:06 44042077                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_crs.cpython-37m-x86_64-linux-gnu.so
7fffae954000-7fffae955000 r--p 00010000 08:06 44042077                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_crs.cpython-37m-x86_64-linux-gnu.so
7fffae955000-7fffae956000 rw-p 00011000 08:06 44042077                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_crs.cpython-37m-x86_64-linux-gnu.so
7fffae956000-7fffae996000 rw-p 00000000 00:00 0
7fffae996000-7fffae999000 r--p 00000000 08:06 44042076                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/schema.cpython-37m-x86_64-linux-gnu.so
7fffae999000-7fffae99f000 r-xp 00003000 08:06 44042076                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/schema.cpython-37m-x86_64-linux-gnu.so
7fffae99f000-7fffae9a0000 r--p 00009000 08:06 44042076                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/schema.cpython-37m-x86_64-linux-gnu.so
7fffae9a0000-7fffae9a1000 ---p 0000a000 08:06 44042076                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/schema.cpython-37m-x86_64-linux-gnu.so
7fffae9a1000-7fffae9a2000 r--p 0000a000 08:06 44042076                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/schema.cpython-37m-x86_64-linux-gnu.so
7fffae9a2000-7fffae9a3000 rw-p 0000b000 08:06 44042076                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/schema.cpython-37m-x86_64-linux-gnu.so
7fffae9a3000-7fffae9e3000 rw-p 00000000 00:00 0
7fffae9e3000-7fffae9eb000 r--p 00000000 08:06 44042082                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_env.cpython-37m-x86_64-linux-gnu.so
7fffae9eb000-7fffaea0f000 r-xp 00008000 08:06 44042082                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_env.cpython-37m-x86_64-linux-gnu.so
7fffaea0f000-7fffaea14000 r--p 0002c000 08:06 44042082                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_env.cpython-37m-x86_64-linux-gnu.so
7fffaea14000-7fffaea15000 r--p 00030000 08:06 44042082                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_env.cpython-37m-x86_64-linux-gnu.so
7fffaea15000-7fffaea19000 rw-p 00031000 08:06 44042082                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_env.cpython-37m-x86_64-linux-gnu.so
7fffaea19000-7fffaea1a000 rw-p 00000000 00:00 0
7fffaea1a000-7fffaea1f000 r--p 00000000 08:06 44042079                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_err.cpython-37m-x86_64-linux-gnu.so
7fffaea1f000-7fffaea2e000 r-xp 00005000 08:06 44042079                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_err.cpython-37m-x86_64-linux-gnu.so
7fffaea2e000-7fffaea31000 r--p 00014000 08:06 44042079                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_err.cpython-37m-x86_64-linux-gnu.so
7fffaea31000-7fffaea32000 ---p 00017000 08:06 44042079                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_err.cpython-37m-x86_64-linux-gnu.so
7fffaea32000-7fffaea33000 r--p 00017000 08:06 44042079                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_err.cpython-37m-x86_64-linux-gnu.so
7fffaea33000-7fffaea35000 rw-p 00018000 08:06 44042079                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_err.cpython-37m-x86_64-linux-gnu.so
7fffaea35000-7fffaea3b000 r--p 00000000 08:06 44042081                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_geometry.cpython-37m-x86_64-linux-gnu.so
7fffaea3b000-7fffaea52000 r-xp 00006000 08:06 44042081                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_geometry.cpython-37m-x86_64-linux-gnu.so
7fffaea52000-7fffaea57000 r--p 0001d000 08:06 44042081                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_geometry.cpython-37m-x86_64-linux-gnu.so
7fffaea57000-7fffaea58000 r--p 00021000 08:06 44042081                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_geometry.cpython-37m-x86_64-linux-gnu.so
7fffaea58000-7fffaea5b000 rw-p 00022000 08:06 44042081                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_geometry.cpython-37m-x86_64-linux-gnu.so
7fffaea5b000-7fffaea72000 r-xp 00000000 08:01 3149218                    /lib/x86_64-linux-gnu/libresolv-2.23.so
7fffaea72000-7fffaec72000 ---p 00017000 08:01 3149218                    /lib/x86_64-linux-gnu/libresolv-2.23.so
7fffaec72000-7fffaec73000 r--p 00017000 08:01 3149218                    /lib/x86_64-linux-gnu/libresolv-2.23.so
7fffaec73000-7fffaec74000 rw-p 00018000 08:01 3149218                    /lib/x86_64-linux-gnu/libresolv-2.23.so
7fffaec74000-7fffaec76000 rw-p 00000000 00:00 0
7fffaec82000-7fffaec86000 r--p 00000000 08:06 44042078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_shim.cpython-37m-x86_64-linux-gnu.so
7fffaec86000-7fffaec93000 r-xp 00004000 08:06 44042078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_shim.cpython-37m-x86_64-linux-gnu.so
7fffaec93000-7fffaec95000 r--p 00011000 08:06 44042078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_shim.cpython-37m-x86_64-linux-gnu.so
7fffaec95000-7fffaec96000 r--p 00012000 08:06 44042078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_shim.cpython-37m-x86_64-linux-gnu.so
7fffaec96000-7fffaec97000 rw-p 00013000 08:06 44042078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/_shim.cpython-37m-x86_64-linux-gnu.so
7fffaec97000-7fffaec9b000 r--p 00000000 08:06 44043804                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5support.so.0.1
7fffaec9b000-7fffaeca0000 r-xp 00004000 08:06 44043804                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5support.so.0.1
7fffaeca0000-7fffaeca2000 r--p 00009000 08:06 44043804                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5support.so.0.1
7fffaeca2000-7fffaeca3000 ---p 0000b000 08:06 44043804                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5support.so.0.1
7fffaeca3000-7fffaeca4000 r--p 0000b000 08:06 44043804                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5support.so.0.1
7fffaeca4000-7fffaeca5000 rw-p 0000c000 08:06 44043804                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5support.so.0.1
7fffaeca5000-7fffaeca7000 r--p 00000000 08:06 44043946                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcom_err.so.3.0
7fffaeca7000-7fffaeca8000 r-xp 00002000 08:06 44043946                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcom_err.so.3.0
7fffaeca8000-7fffaeca9000 r--p 00003000 08:06 44043946                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcom_err.so.3.0
7fffaeca9000-7fffaecaa000 r--p 00003000 08:06 44043946                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcom_err.so.3.0
7fffaecaa000-7fffaecab000 rw-p 00004000 08:06 44043946                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcom_err.so.3.0
7fffaecab000-7fffaecb1000 r--p 00000000 08:06 44043805                   /home/chadi/miniconda3/envs/ocp-covid/lib/libk5crypto.so.3.1
7fffaecb1000-7fffaecc0000 r-xp 00006000 08:06 44043805                   /home/chadi/miniconda3/envs/ocp-covid/lib/libk5crypto.so.3.1
7fffaecc0000-7fffaecc5000 r--p 00015000 08:06 44043805                   /home/chadi/miniconda3/envs/ocp-covid/lib/libk5crypto.so.3.1
7fffaecc5000-7fffaecc6000 ---p 0001a000 08:06 44043805                   /home/chadi/miniconda3/envs/ocp-covid/lib/libk5crypto.so.3.1
7fffaecc6000-7fffaecc8000 r--p 0001a000 08:06 44043805                   /home/chadi/miniconda3/envs/ocp-covid/lib/libk5crypto.so.3.1
7fffaecc8000-7fffaecc9000 rw-p 0001c000 08:06 44043805                   /home/chadi/miniconda3/envs/ocp-covid/lib/libk5crypto.so.3.1
7fffaecc9000-7fffaecca000 rw-p 00000000 00:00 0
7fffaecca000-7fffaecef000 r--p 00000000 08:06 44182054                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5.so.3.3
7fffaecef000-7fffaed48000 r-xp 00025000 08:06 44182054                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5.so.3.3
7fffaed48000-7fffaed90000 r--p 0007e000 08:06 44182054                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5.so.3.3
7fffaed90000-7fffaed91000 ---p 000c6000 08:06 44182054                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5.so.3.3
7fffaed91000-7fffaed9f000 r--p 000c6000 08:06 44182054                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5.so.3.3
7fffaed9f000-7fffaeda1000 rw-p 000d4000 08:06 44182054                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkrb5.so.3.3
7fffaeda1000-7fffaedab000 r--p 00000000 08:06 43647843                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssh2.so.1.0.1
7fffaedab000-7fffaedd4000 r-xp 0000a000 08:06 43647843                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssh2.so.1.0.1
7fffaedd4000-7fffaede1000 r--p 00033000 08:06 43647843                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssh2.so.1.0.1
7fffaede1000-7fffaede3000 r--p 0003f000 08:06 43647843                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssh2.so.1.0.1
7fffaede3000-7fffaede4000 rw-p 00041000 08:06 43647843                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssh2.so.1.0.1
7fffaede4000-7fffaedf1000 r--p 00000000 08:06 44182050                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgssapi_krb5.so.2.2
7fffaedf1000-7fffaee23000 r-xp 0000d000 08:06 44182050                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgssapi_krb5.so.2.2
7fffaee23000-7fffaee2f000 r--p 0003f000 08:06 44182050                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgssapi_krb5.so.2.2
7fffaee2f000-7fffaee30000 ---p 0004b000 08:06 44182050                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgssapi_krb5.so.2.2
7fffaee30000-7fffaee32000 r--p 0004b000 08:06 44182050                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgssapi_krb5.so.2.2
7fffaee32000-7fffaee33000 rw-p 0004d000 08:06 44182050                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgssapi_krb5.so.2.2
7fffaee33000-7fffaee52000 r--p 00000000 08:06 43389024                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssl.so.1.1
7fffaee52000-7fffaee9c000 r-xp 0001f000 08:06 43389024                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssl.so.1.1
7fffaee9c000-7fffaeeb5000 r--p 00069000 08:06 43389024                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssl.so.1.1
7fffaeeb5000-7fffaeeb6000 ---p 00082000 08:06 43389024                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssl.so.1.1
7fffaeeb6000-7fffaeebf000 r--p 00082000 08:06 43389024                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssl.so.1.1
7fffaeebf000-7fffaeec3000 rw-p 0008b000 08:06 43389024                   /home/chadi/miniconda3/envs/ocp-covid/lib/libssl.so.1.1
7fffaeec3000-7fffaeee4000 r-xp 00000000 08:06 44441810                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlconvenience.so.1.3.0
7fffaeee4000-7fffaf0e4000 ---p 00021000 08:06 44441810                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlconvenience.so.1.3.0
7fffaf0e4000-7fffaf0e5000 r--p 00021000 08:06 44441810                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlconvenience.so.1.3.0
7fffaf0e5000-7fffaf0e6000 rw-p 00022000 08:06 44441810                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlconvenience.so.1.3.0
7fffaf0e6000-7fffaf103000 r-xp 00000000 08:06 44441807                   /home/chadi/miniconda3/envs/ocp-covid/lib/liburiparser.so
7fffaf103000-7fffaf303000 ---p 0001d000 08:06 44441807                   /home/chadi/miniconda3/envs/ocp-covid/lib/liburiparser.so
7fffaf303000-7fffaf304000 r--p 0001d000 08:06 44441807                   /home/chadi/miniconda3/envs/ocp-covid/lib/liburiparser.so
7fffaf304000-7fffaf305000 rw-p 0001e000 08:06 44441807                   /home/chadi/miniconda3/envs/ocp-covid/lib/liburiparser.so
7fffaf305000-7fffaf312000 r-xp 00000000 08:06 44441806                   /home/chadi/miniconda3/envs/ocp-covid/lib/libminizip.so
7fffaf312000-7fffaf511000 ---p 0000d000 08:06 44441806                   /home/chadi/miniconda3/envs/ocp-covid/lib/libminizip.so
7fffaf511000-7fffaf512000 r--p 0000c000 08:06 44441806                   /home/chadi/miniconda3/envs/ocp-covid/lib/libminizip.so
7fffaf512000-7fffaf513000 rw-p 0000d000 08:06 44441806                   /home/chadi/miniconda3/envs/ocp-covid/lib/libminizip.so
7fffaf513000-7fffaf595000 r--p 00000000 08:06 44043103                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos-3.8.0.so
7fffaf595000-7fffaf681000 r-xp 00082000 08:06 44043103                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos-3.8.0.so
7fffaf681000-7fffaf6bf000 r--p 0016e000 08:06 44043103                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos-3.8.0.so
7fffaf6bf000-7fffaf6ca000 r--p 001ab000 08:06 44043103                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos-3.8.0.so
7fffaf6ca000-7fffaf6cb000 rw-p 001b6000 08:06 44043103                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos-3.8.0.so
7fffaf6cb000-7fffaf70c000 r-xp 00000000 08:06 44182081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfontconfig.so.1.11.1
7fffaf70c000-7fffaf90c000 ---p 00041000 08:06 44182081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfontconfig.so.1.11.1
7fffaf90c000-7fffaf90e000 r--p 00041000 08:06 44182081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfontconfig.so.1.11.1
7fffaf90e000-7fffaf90f000 rw-p 00043000 08:06 44182081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfontconfig.so.1.11.1
7fffaf90f000-7fffaf932000 r-xp 00000000 08:06 44434753                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblz4.so.1.8.1
7fffaf932000-7fffafb31000 ---p 00023000 08:06 44434753                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblz4.so.1.8.1
7fffafb31000-7fffafb32000 r--p 00022000 08:06 44434753                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblz4.so.1.8.1
7fffafb32000-7fffafb33000 rw-p 00023000 08:06 44434753                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblz4.so.1.8.1
7fffafb33000-7fffafb35000 r--p 00000000 08:06 44046496                   /home/chadi/miniconda3/envs/ocp-covid/lib/libbz2.so.1.0.8
7fffafb35000-7fffafb43000 r-xp 00002000 08:06 44046496                   /home/chadi/miniconda3/envs/ocp-covid/lib/libbz2.so.1.0.8
7fffafb43000-7fffafb45000 r--p 00010000 08:06 44046496                   /home/chadi/miniconda3/envs/ocp-covid/lib/libbz2.so.1.0.8
7fffafb45000-7fffafb46000 r--p 00011000 08:06 44046496                   /home/chadi/miniconda3/envs/ocp-covid/lib/libbz2.so.1.0.8
7fffafb46000-7fffafb47000 rw-p 00012000 08:06 44046496                   /home/chadi/miniconda3/envs/ocp-covid/lib/libbz2.so.1.0.8
7fffafb47000-7fffafb56000 r--p 00000000 08:06 44181959                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtbb.so.2
7fffafb56000-7fffafb74000 r-xp 0000f000 08:06 44181959                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtbb.so.2
7fffafb74000-7fffafb81000 r--p 0002d000 08:06 44181959                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtbb.so.2
7fffafb81000-7fffafb83000 r--p 00039000 08:06 44181959                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtbb.so.2
7fffafb83000-7fffafb85000 rw-p 0003b000 08:06 44181959                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtbb.so.2
7fffafb85000-7fffafb87000 rw-p 00000000 00:00 0
7fffafb87000-7fffafbb0000 r--p 00000000 08:06 44440722                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_cpp.so.103.0.0
7fffafbb0000-7fffafbe3000 r-xp 00029000 08:06 44440722                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_cpp.so.103.0.0
7fffafbe3000-7fffafbf8000 r--p 0005c000 08:06 44440722                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_cpp.so.103.0.0
7fffafbf8000-7fffafbfc000 r--p 00070000 08:06 44440722                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_cpp.so.103.0.0
7fffafbfc000-7fffafbfd000 rw-p 00074000 08:06 44440722                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_cpp.so.103.0.0
7fffafbfd000-7fffafc04000 r--p 00000000 08:06 44440757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_hl.so.100.1.1
7fffafc04000-7fffafc19000 r-xp 00007000 08:06 44440757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_hl.so.100.1.1
7fffafc19000-7fffafc1f000 r--p 0001c000 08:06 44440757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_hl.so.100.1.1
7fffafc1f000-7fffafc20000 ---p 00022000 08:06 44440757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_hl.so.100.1.1
7fffafc20000-7fffafc21000 r--p 00022000 08:06 44440757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_hl.so.100.1.1
7fffafc21000-7fffafc22000 rw-p 00023000 08:06 44440757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5_hl.so.100.1.1
7fffafc22000-7fffafc23000 rw-p 00000000 00:00 0
7fffafc23000-7fffafc66000 r--p 00000000 08:06 44181995                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5.so.103.0.0
7fffafc66000-7fffafef3000 r-xp 00043000 08:06 44181995                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5.so.103.0.0
7fffafef3000-7fffaffa7000 r--p 002d0000 08:06 44181995                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5.so.103.0.0
7fffaffa7000-7fffaffb4000 r--p 00383000 08:06 44181995                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5.so.103.0.0
7fffaffb4000-7fffaffb7000 rw-p 00390000 08:06 44181995                   /home/chadi/miniconda3/envs/ocp-covid/lib/libhdf5.so.103.0.0
7fffaffb7000-7fffaffb9000 rw-p 00000000 00:00 0
7fffaffb9000-7fffaffc6000 r--p 00000000 08:06 44303178                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkea.so.1.4.7
7fffaffc6000-7fffb0016000 r-xp 0000d000 08:06 44303178                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkea.so.1.4.7
7fffb0016000-7fffb0024000 r--p 0005d000 08:06 44303178                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkea.so.1.4.7
7fffb0024000-7fffb0026000 r--p 0006a000 08:06 44303178                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkea.so.1.4.7
7fffb0026000-7fffb0027000 rw-p 0006c000 08:06 44303178                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkea.so.1.4.7
7fffb0027000-7fffb0029000 rw-p 00000000 00:00 0
7fffb0029000-7fffb005a000 r--p 00000000 08:06 43649046                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxml2.so.2.9.10
7fffb005a000-7fffb013b000 r-xp 00031000 08:06 43649046                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxml2.so.2.9.10
7fffb013b000-7fffb0184000 r--p 00112000 08:06 43649046                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxml2.so.2.9.10
7fffb0184000-7fffb018d000 r--p 0015a000 08:06 43649046                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxml2.so.2.9.10
7fffb018d000-7fffb018e000 rw-p 00163000 08:06 43649046                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxml2.so.2.9.10
7fffb018e000-7fffb0190000 rw-p 00000000 00:00 0
7fffb0190000-7fffb019c000 r--p 00000000 08:06 44182135                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcurl.so.4.6.0
7fffb019c000-7fffb01f8000 r-xp 0000c000 08:06 44182135                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcurl.so.4.6.0
7fffb01f8000-7fffb0212000 r--p 00068000 08:06 44182135                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcurl.so.4.6.0
7fffb0212000-7fffb0215000 r--p 00081000 08:06 44182135                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcurl.so.4.6.0
7fffb0215000-7fffb0216000 rw-p 00084000 08:06 44182135                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcurl.so.4.6.0
7fffb0216000-7fffb0218000 r--p 00000000 08:06 44442157                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpcre.so.1.2.12
7fffb0218000-7fffb0242000 r-xp 00002000 08:06 44442157                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpcre.so.1.2.12
7fffb0242000-7fffb025b000 r--p 0002c000 08:06 44442157                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpcre.so.1.2.12
7fffb025b000-7fffb025c000 r--p 00044000 08:06 44442157                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpcre.so.1.2.12
7fffb025c000-7fffb025d000 rw-p 00045000 08:06 44442157                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpcre.so.1.2.12
7fffb025d000-7fffb0277000 r--p 00000000 08:06 44307156                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialite.so.7.1.0
7fffb0277000-7fffb05a0000 r-xp 0001a000 08:06 44307156                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialite.so.7.1.0
7fffb05a0000-7fffb07d7000 r--p 00343000 08:06 44307156                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialite.so.7.1.0
7fffb07d7000-7fffb07da000 r--p 00579000 08:06 44307156                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialite.so.7.1.0
7fffb07da000-7fffb07dc000 rw-p 0057c000 08:06 44307156                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialite.so.7.1.0
7fffb07dc000-7fffb07dd000 rw-p 00000000 00:00 0
7fffb07dd000-7fffb07e0000 r-xp 00000000 08:06 44309638                   /home/chadi/miniconda3/envs/ocp-covid/lib/libuuid.so.1.0.0
7fffb07e0000-7fffb09df000 ---p 00003000 08:06 44309638                   /home/chadi/miniconda3/envs/ocp-covid/lib/libuuid.so.1.0.0
7fffb09df000-7fffb09e0000 r--p 00002000 08:06 44309638                   /home/chadi/miniconda3/envs/ocp-covid/lib/libuuid.so.1.0.0
7fffb09e0000-7fffb09e1000 rw-p 00003000 08:06 44309638                   /home/chadi/miniconda3/envs/ocp-covid/lib/libuuid.so.1.0.0
7fffb09e1000-7fffb22e2000 r--p 00000000 08:06 43648343                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicudata.so.58.2
7fffb22e2000-7fffb22e3000 r--p 01900000 08:06 43648343                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicudata.so.58.2
7fffb22e3000-7fffb2341000 r--p 00000000 08:06 43648341                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicuuc.so.58.2
7fffb2341000-7fffb240d000 r-xp 0005e000 08:06 43648341                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicuuc.so.58.2
7fffb240d000-7fffb2482000 r--p 0012a000 08:06 43648341                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicuuc.so.58.2
7fffb2482000-7fffb2494000 r--p 0019e000 08:06 43648341                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicuuc.so.58.2
7fffb2494000-7fffb2495000 rw-p 001b0000 08:06 43648341                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicuuc.so.58.2
7fffb2495000-7fffb2496000 rw-p 00000000 00:00 0
7fffb2496000-7fffb2557000 r--p 00000000 08:06 43648342                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicui18n.so.58.2
7fffb2557000-7fffb2696000 r-xp 000c1000 08:06 43648342                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicui18n.so.58.2
7fffb2696000-7fffb2705000 r--p 00200000 08:06 43648342                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicui18n.so.58.2
7fffb2705000-7fffb2706000 ---p 0026f000 08:06 43648342                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicui18n.so.58.2
7fffb2706000-7fffb2714000 r--p 0026f000 08:06 43648342                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicui18n.so.58.2
7fffb2714000-7fffb2715000 rw-p 0027d000 08:06 43648342                   /home/chadi/miniconda3/envs/ocp-covid/lib/libicui18n.so.58.2
7fffb2715000-7fffb2755000 r-xp 00000000 08:06 44437643                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapclient.so.6.1.7
7fffb2755000-7fffb2954000 ---p 00040000 08:06 44437643                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapclient.so.6.1.7
7fffb2954000-7fffb2956000 r--p 0003f000 08:06 44437643                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapclient.so.6.1.7
7fffb2956000-7fffb2957000 rw-p 00041000 08:06 44437643                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapclient.so.6.1.7
7fffb2957000-7fffb2965000 r-xp 00000000 08:06 44437642                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapserver.so.7.6.7
7fffb2965000-7fffb2b65000 ---p 0000e000 08:06 44437642                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapserver.so.7.6.7
7fffb2b65000-7fffb2b66000 r--p 0000e000 08:06 44437642                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapserver.so.7.6.7
7fffb2b66000-7fffb2b67000 rw-p 0000f000 08:06 44437642                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdapserver.so.7.6.7
7fffb2b67000-7fffb2ce6000 r-xp 00000000 08:06 44182385                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdap.so.25.0.1
7fffb2ce6000-7fffb2ee6000 ---p 0017f000 08:06 44182385                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdap.so.25.0.1
7fffb2ee6000-7fffb2eef000 r--p 0017f000 08:06 44182385                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdap.so.25.0.1
7fffb2eef000-7fffb2ef0000 rw-p 00188000 08:06 44182385                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdap.so.25.0.1
7fffb2ef0000-7fffb2ef1000 rw-p 00000000 00:00 0
7fffb2ef1000-7fffb2f62000 r--p 00000000 08:06 44182151                   /home/chadi/miniconda3/envs/ocp-covid/lib/libproj.so.15.2.1
7fffb2f62000-7fffb30e2000 r-xp 00071000 08:06 44182151                   /home/chadi/miniconda3/envs/ocp-covid/lib/libproj.so.15.2.1
7fffb30e2000-7fffb314c000 r--p 001f1000 08:06 44182151                   /home/chadi/miniconda3/envs/ocp-covid/lib/libproj.so.15.2.1
7fffb314c000-7fffb314d000 ---p 0025b000 08:06 44182151                   /home/chadi/miniconda3/envs/ocp-covid/lib/libproj.so.15.2.1
7fffb314d000-7fffb3163000 r--p 0025b000 08:06 44182151                   /home/chadi/miniconda3/envs/ocp-covid/lib/libproj.so.15.2.1
7fffb3163000-7fffb3165000 rw-p 00271000 08:06 44182151                   /home/chadi/miniconda3/envs/ocp-covid/lib/libproj.so.15.2.1
7fffb3165000-7fffb3167000 rw-p 00000000 00:00 0
7fffb3167000-7fffb3172000 r--p 00000000 08:06 44182142                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpq.so.5.11
7fffb3172000-7fffb318f000 r-xp 0000b000 08:06 44182142                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpq.so.5.11
7fffb318f000-7fffb31ad000 r--p 00028000 08:06 44182142                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpq.so.5.11
7fffb31ad000-7fffb31ae000 ---p 00046000 08:06 44182142                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpq.so.5.11
7fffb31ae000-7fffb31b1000 r--p 00046000 08:06 44182142                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpq.so.5.11
7fffb31b1000-7fffb31b2000 rw-p 00049000 08:06 44182142                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpq.so.5.11
7fffb31b2000-7fffb31d5000 r--p 00000000 08:06 44433587                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcfitsio.so.8.3.47
7fffb31d5000-7fffb3303000 r-xp 00023000 08:06 44433587                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcfitsio.so.8.3.47
7fffb3303000-7fffb3349000 r--p 00151000 08:06 44433587                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcfitsio.so.8.3.47
7fffb3349000-7fffb334a000 ---p 00197000 08:06 44433587                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcfitsio.so.8.3.47
7fffb334a000-7fffb334d000 r--p 00197000 08:06 44433587                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcfitsio.so.8.3.47
7fffb334d000-7fffb334f000 rw-p 0019a000 08:06 44433587                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcfitsio.so.8.3.47
7fffb334f000-7fffb34df000 rw-p 00000000 00:00 0
7fffb34df000-7fffb34ef000 r--p 00000000 08:06 44304068                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeotiff.so.5.0.0
7fffb34ef000-7fffb34fb000 r-xp 00010000 08:06 44304068                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeotiff.so.5.0.0
7fffb34fb000-7fffb3509000 r--p 0001c000 08:06 44304068                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeotiff.so.5.0.0
7fffb3509000-7fffb350a000 ---p 0002a000 08:06 44304068                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeotiff.so.5.0.0
7fffb350a000-7fffb3513000 r--p 0002a000 08:06 44304068                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeotiff.so.5.0.0
7fffb3513000-7fffb3514000 rw-p 00033000 08:06 44304068                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeotiff.so.5.0.0
7fffb3514000-7fffb351c000 r-xp 00000000 08:06 44304914                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgif.so.7.0.0
7fffb351c000-7fffb371c000 ---p 00008000 08:06 44304914                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgif.so.7.0.0
7fffb371c000-7fffb371d000 r--p 00008000 08:06 44304914                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgif.so.7.0.0
7fffb371d000-7fffb371e000 rw-p 00009000 08:06 44304914                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgif.so.7.0.0
7fffb371e000-7fffb379d000 r-xp 00000000 08:06 44434560                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdf.so.0.0.0
7fffb379d000-7fffb399c000 ---p 0007f000 08:06 44434560                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdf.so.0.0.0
7fffb399c000-7fffb399f000 r--p 0007e000 08:06 44434560                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdf.so.0.0.0
7fffb399f000-7fffb39a0000 rw-p 00081000 08:06 44434560                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdf.so.0.0.0
7fffb39a0000-7fffb39c9000 rw-p 00000000 00:00 0
7fffb39c9000-7fffb39f0000 r-xp 00000000 08:06 44434559                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmfhdf.so.0.0.0
7fffb39f0000-7fffb3bf0000 ---p 00027000 08:06 44434559                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmfhdf.so.0.0.0
7fffb3bf0000-7fffb3bf1000 r--p 00027000 08:06 44434559                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmfhdf.so.0.0.0
7fffb3bf1000-7fffb3bf2000 rw-p 00028000 08:06 44434559                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmfhdf.so.0.0.0
7fffb3bf2000-7fffb3bf3000 rw-p 00000000 00:00 0
7fffb3bf3000-7fffb3c1a000 r--p 00000000 08:06 44436761                   /home/chadi/miniconda3/envs/ocp-covid/lib/libnetcdf.so.13
7fffb3c1a000-7fffb3cb0000 r-xp 00027000 08:06 44436761                   /home/chadi/miniconda3/envs/ocp-covid/lib/libnetcdf.so.13
7fffb3cb0000-7fffb3d37000 r--p 000bd000 08:06 44436761                   /home/chadi/miniconda3/envs/ocp-covid/lib/libnetcdf.so.13
7fffb3d37000-7fffb3d38000 ---p 00144000 08:06 44436761                   /home/chadi/miniconda3/envs/ocp-covid/lib/libnetcdf.so.13
7fffb3d38000-7fffb3d3c000 r--p 00144000 08:06 44436761                   /home/chadi/miniconda3/envs/ocp-covid/lib/libnetcdf.so.13
7fffb3d3c000-7fffb3d3f000 rw-p 00148000 08:06 44436761                   /home/chadi/miniconda3/envs/ocp-covid/lib/libnetcdf.so.13
7fffb3d3f000-7fffb3d49000 rw-p 00000000 00:00 0
7fffb3d49000-7fffb3da1000 r-xp 00000000 08:06 44303216                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenjp2.so.2.3.0
7fffb3da1000-7fffb3fa0000 ---p 00058000 08:06 44303216                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenjp2.so.2.3.0
7fffb3fa0000-7fffb3fa2000 r--p 00057000 08:06 44303216                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenjp2.so.2.3.0
7fffb3fa2000-7fffb3fa3000 rw-p 00059000 08:06 44303216                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenjp2.so.2.3.0
7fffb3fa3000-7fffb40bf000 r--p 00000000 08:06 43648434                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxerces-c-3.2.so
7fffb40bf000-7fffb4236000 r-xp 0011c000 08:06 43648434                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxerces-c-3.2.so
7fffb4236000-7fffb42eb000 r--p 00293000 08:06 43648434                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxerces-c-3.2.so
7fffb42eb000-7fffb4303000 r--p 00347000 08:06 43648434                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxerces-c-3.2.so
7fffb4303000-7fffb4325000 rw-p 0035f000 08:06 43648434                   /home/chadi/miniconda3/envs/ocp-covid/lib/libxerces-c-3.2.so
7fffb4325000-7fffb4329000 r--p 00000000 08:06 43522763                   /home/chadi/miniconda3/envs/ocp-covid/lib/libexpat.so.1.6.12
7fffb4329000-7fffb4346000 r-xp 00004000 08:06 43522763                   /home/chadi/miniconda3/envs/ocp-covid/lib/libexpat.so.1.6.12
7fffb4346000-7fffb4355000 r--p 00021000 08:06 43522763                   /home/chadi/miniconda3/envs/ocp-covid/lib/libexpat.so.1.6.12
7fffb4355000-7fffb4356000 ---p 00030000 08:06 43522763                   /home/chadi/miniconda3/envs/ocp-covid/lib/libexpat.so.1.6.12
7fffb4356000-7fffb4358000 r--p 00030000 08:06 43522763                   /home/chadi/miniconda3/envs/ocp-covid/lib/libexpat.so.1.6.12
7fffb4358000-7fffb4359000 rw-p 00032000 08:06 43522763                   /home/chadi/miniconda3/envs/ocp-covid/lib/libexpat.so.1.6.12
7fffb4359000-7fffb4364000 r-xp 00000000 08:06 44441808                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlregionator.so.1.3.0
7fffb4364000-7fffb4563000 ---p 0000b000 08:06 44441808                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlregionator.so.1.3.0
7fffb4563000-7fffb4564000 r--p 0000a000 08:06 44441808                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlregionator.so.1.3.0
7fffb4564000-7fffb4565000 rw-p 0000b000 08:06 44441808                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlregionator.so.1.3.0
7fffb4565000-7fffb4577000 r-xp 00000000 08:06 44441809                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlxsd.so.1.3.0
7fffb4577000-7fffb4777000 ---p 00012000 08:06 44441809                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlxsd.so.1.3.0
7fffb4777000-7fffb4778000 r--p 00012000 08:06 44441809                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlxsd.so.1.3.0
7fffb4778000-7fffb4779000 rw-p 00013000 08:06 44441809                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlxsd.so.1.3.0
7fffb4779000-7fffb47b6000 r-xp 00000000 08:06 44441812                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlengine.so.1.3.0
7fffb47b6000-7fffb49b5000 ---p 0003d000 08:06 44441812                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlengine.so.1.3.0
7fffb49b5000-7fffb49b7000 r--p 0003c000 08:06 44441812                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlengine.so.1.3.0
7fffb49b7000-7fffb49b8000 rw-p 0003e000 08:06 44441812                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlengine.so.1.3.0
7fffb49b8000-7fffb4a78000 r-xp 00000000 08:06 44441813                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmldom.so.1.3.0
7fffb4a78000-7fffb4c77000 ---p 000c0000 08:06 44441813                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmldom.so.1.3.0
7fffb4c77000-7fffb4c7e000 r--p 000bf000 08:06 44441813                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmldom.so.1.3.0
7fffb4c7e000-7fffb4c80000 rw-p 000c6000 08:06 44441813                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmldom.so.1.3.0
7fffb4c80000-7fffb4ca2000 r-xp 00000000 08:06 44441811                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlbase.so.1.3.0
7fffb4ca2000-7fffb4ea2000 ---p 00022000 08:06 44441811                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlbase.so.1.3.0
7fffb4ea2000-7fffb4ea3000 r--p 00022000 08:06 44441811                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlbase.so.1.3.0
7fffb4ea3000-7fffb4ea4000 rw-p 00023000 08:06 44441811                   /home/chadi/miniconda3/envs/ocp-covid/lib/libkmlbase.so.1.3.0
7fffb4ea4000-7fffb4eb4000 r--p 00000000 08:06 43387081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libsqlite3.so.0.8.6
7fffb4eb4000-7fffb4f90000 r-xp 00010000 08:06 43387081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libsqlite3.so.0.8.6
7fffb4f90000-7fffb4fc1000 r--p 000ec000 08:06 43387081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libsqlite3.so.0.8.6
7fffb4fc1000-7fffb4fc2000 ---p 0011d000 08:06 43387081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libsqlite3.so.0.8.6
7fffb4fc2000-7fffb4fc6000 r--p 0011d000 08:06 43387081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libsqlite3.so.0.8.6
7fffb4fc6000-7fffb4fc9000 rw-p 00121000 08:06 43387081                   /home/chadi/miniconda3/envs/ocp-covid/lib/libsqlite3.so.0.8.6
7fffb4fc9000-7fffb4fdc000 r--p 00000000 08:06 44043105                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos_c.so.1.13.1
7fffb4fdc000-7fffb4ff9000 r-xp 00013000 08:06 44043105                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos_c.so.1.13.1
7fffb4ff9000-7fffb5003000 r--p 00030000 08:06 44043105                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos_c.so.1.13.1
7fffb5003000-7fffb5005000 r--p 00039000 08:06 44043105                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos_c.so.1.13.1
7fffb5005000-7fffb5006000 rw-p 0003b000 08:06 44043105                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgeos_c.so.1.13.1
7fffb5006000-7fffb500e000 r-xp 00000000 08:06 44441489                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreexl.so.1.1.0
7fffb500e000-7fffb520e000 ---p 00008000 08:06 44441489                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreexl.so.1.1.0
7fffb520e000-7fffb520f000 r--p 00008000 08:06 44441489                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreexl.so.1.1.0
7fffb520f000-7fffb5210000 rw-p 00009000 08:06 44441489                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreexl.so.1.1.0
7fffb5210000-7fffb521e000 r-xp 00000000 08:06 43392044                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjson-c.so.4.0.0
7fffb521e000-7fffb541e000 ---p 0000e000 08:06 43392044                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjson-c.so.4.0.0
7fffb541e000-7fffb541f000 r--p 0000e000 08:06 43392044                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjson-c.so.4.0.0
7fffb541f000-7fffb5420000 rw-p 0000f000 08:06 43392044                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjson-c.so.4.0.0
7fffb5420000-7fffb5675000 r-xp 00000000 08:06 44182396                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpoppler.so.76.0.0
7fffb5675000-7fffb5875000 ---p 00255000 08:06 44182396                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpoppler.so.76.0.0
7fffb5875000-7fffb5896000 r--p 00255000 08:06 44182396                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpoppler.so.76.0.0
7fffb5896000-7fffb58bc000 rw-p 00276000 08:06 44182396                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpoppler.so.76.0.0
7fffb58bc000-7fffb59ec000 r--p 00000000 08:06 44436447                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiledb.so.1.6.3
7fffb59ec000-7fffb5eb3000 r-xp 00130000 08:06 44436447                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiledb.so.1.6.3
7fffb5eb3000-7fffb600e000 r--p 005f7000 08:06 44436447                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiledb.so.1.6.3
7fffb600e000-7fffb600f000 ---p 00752000 08:06 44436447                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiledb.so.1.6.3
7fffb600f000-7fffb6031000 r--p 00752000 08:06 44436447                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiledb.so.1.6.3
7fffb6031000-7fffb6032000 rw-p 00774000 08:06 44436447                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiledb.so.1.6.3
7fffb6032000-7fffb6049000 rw-p 00000000 00:00 0
7fffb6049000-7fffb64aa000 r--p 00000000 08:06 44182416                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgdal.so.26.0.2
7fffb64aa000-7fffb7054000 r-xp 00461000 08:06 44182416                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgdal.so.26.0.2
7fffb7054000-7fffb73a4000 r--p 0100b000 08:06 44182416                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgdal.so.26.0.2
7fffb73a4000-7fffb73a5000 ---p 0135b000 08:06 44182416                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgdal.so.26.0.2
7fffb73a5000-7fffb749c000 r--p 0135b000 08:06 44182416                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgdal.so.26.0.2
7fffb749c000-7fffb749d000 rw-p 01452000 08:06 44182416                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgdal.so.26.0.2
7fffb749d000-7fffb74e6000 rw-p 00000000 00:00 0
7fffb74e6000-7fffb74f3000 r--p 00000000 08:06 44042083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/ogrext.cpython-37m-x86_64-linux-gnu.so
7fffb74f3000-7fffb7557000 r-xp 0000d000 08:06 44042083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/ogrext.cpython-37m-x86_64-linux-gnu.so
7fffb7557000-7fffb755f000 r--p 00071000 08:06 44042083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/ogrext.cpython-37m-x86_64-linux-gnu.so
7fffb755f000-7fffb7560000 ---p 00079000 08:06 44042083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/ogrext.cpython-37m-x86_64-linux-gnu.so
7fffb7560000-7fffb7561000 r--p 00079000 08:06 44042083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/ogrext.cpython-37m-x86_64-linux-gnu.so
7fffb7561000-7fffb7568000 rw-p 0007a000 08:06 44042083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/fiona/ogrext.cpython-37m-x86_64-linux-gnu.so
7fffb7568000-7fffb756a000 rw-p 00000000 00:00 0
7fffb756a000-7fffb756e000 r--p 00000000 08:06 44044105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/testing.cpython-37m-x86_64-linux-gnu.so
7fffb756e000-7fffb757c000 r-xp 00004000 08:06 44044105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/testing.cpython-37m-x86_64-linux-gnu.so
7fffb757c000-7fffb757e000 r--p 00012000 08:06 44044105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/testing.cpython-37m-x86_64-linux-gnu.so
7fffb757e000-7fffb757f000 r--p 00013000 08:06 44044105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/testing.cpython-37m-x86_64-linux-gnu.so
7fffb757f000-7fffb7581000 rw-p 00014000 08:06 44044105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/testing.cpython-37m-x86_64-linux-gnu.so
7fffb7581000-7fffb7701000 rw-p 00000000 00:00 0
7fffb7701000-7fffb7707000 r--p 00000000 08:06 44044104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/json.cpython-37m-x86_64-linux-gnu.so
7fffb7707000-7fffb7712000 r-xp 00006000 08:06 44044104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/json.cpython-37m-x86_64-linux-gnu.so
7fffb7712000-7fffb7716000 r--p 00011000 08:06 44044104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/json.cpython-37m-x86_64-linux-gnu.so
7fffb7716000-7fffb7717000 ---p 00015000 08:06 44044104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/json.cpython-37m-x86_64-linux-gnu.so
7fffb7717000-7fffb7718000 r--p 00015000 08:06 44044104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/json.cpython-37m-x86_64-linux-gnu.so
7fffb7718000-7fffb7719000 rw-p 00016000 08:06 44044104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/json.cpython-37m-x86_64-linux-gnu.so
7fffb7719000-7fffb7759000 rw-p 00000000 00:00 0
7fffb7759000-7fffb7766000 r--p 00000000 08:06 44044132                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/parsers.cpython-37m-x86_64-linux-gnu.so
7fffb7766000-7fffb77bd000 r-xp 0000d000 08:06 44044132                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/parsers.cpython-37m-x86_64-linux-gnu.so
7fffb77bd000-7fffb77ca000 r--p 00064000 08:06 44044132                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/parsers.cpython-37m-x86_64-linux-gnu.so
7fffb77ca000-7fffb77cb000 r--p 00070000 08:06 44044132                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/parsers.cpython-37m-x86_64-linux-gnu.so
7fffb77cb000-7fffb77d1000 rw-p 00071000 08:06 44044132                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/parsers.cpython-37m-x86_64-linux-gnu.so
7fffb77d1000-7fffb7893000 rw-p 00000000 00:00 0
7fffb7893000-7fffb789b000 r--p 00000000 08:06 44044125                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reduction.cpython-37m-x86_64-linux-gnu.so
7fffb789b000-7fffb78cb000 r-xp 00008000 08:06 44044125                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reduction.cpython-37m-x86_64-linux-gnu.so
7fffb78cb000-7fffb78d4000 r--p 00038000 08:06 44044125                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reduction.cpython-37m-x86_64-linux-gnu.so
7fffb78d4000-7fffb78d5000 ---p 00041000 08:06 44044125                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reduction.cpython-37m-x86_64-linux-gnu.so
7fffb78d5000-7fffb78d6000 r--p 00041000 08:06 44044125                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reduction.cpython-37m-x86_64-linux-gnu.so
7fffb78d6000-7fffb78da000 rw-p 00042000 08:06 44044125                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reduction.cpython-37m-x86_64-linux-gnu.so
7fffb78da000-7fffb78db000 rw-p 00000000 00:00 0
7fffb78db000-7fffb78e8000 r--p 00000000 08:06 44044138                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/groupby.cpython-37m-x86_64-linux-gnu.so
7fffb78e8000-7fffb79d3000 r-xp 0000d000 08:06 44044138                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/groupby.cpython-37m-x86_64-linux-gnu.so
7fffb79d3000-7fffb79e5000 r--p 000f8000 08:06 44044138                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/groupby.cpython-37m-x86_64-linux-gnu.so
7fffb79e5000-7fffb79e6000 r--p 00109000 08:06 44044138                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/groupby.cpython-37m-x86_64-linux-gnu.so
7fffb79e6000-7fffb79ef000 rw-p 0010a000 08:06 44044138                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/groupby.cpython-37m-x86_64-linux-gnu.so
7fffb79ef000-7fffb7ab1000 rw-p 00000000 00:00 0
7fffb7ab1000-7fffb7ab7000 r--p 00000000 08:06 44044108                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/indexers.cpython-37m-x86_64-linux-gnu.so
7fffb7ab7000-7fffb7acf000 r-xp 00006000 08:06 44044108                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/indexers.cpython-37m-x86_64-linux-gnu.so
7fffb7acf000-7fffb7ad5000 r--p 0001e000 08:06 44044108                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/indexers.cpython-37m-x86_64-linux-gnu.so
7fffb7ad5000-7fffb7ad6000 r--p 00023000 08:06 44044108                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/indexers.cpython-37m-x86_64-linux-gnu.so
7fffb7ad6000-7fffb7ad9000 rw-p 00024000 08:06 44044108                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/indexers.cpython-37m-x86_64-linux-gnu.so
7fffb7ad9000-7fffb7b19000 rw-p 00000000 00:00 0
7fffb7b19000-7fffb7b23000 r--p 00000000 08:06 44044129                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/aggregations.cpython-37m-x86_64-linux-gnu.so
7fffb7b23000-7fffb7b68000 r-xp 0000a000 08:06 44044129                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/aggregations.cpython-37m-x86_64-linux-gnu.so
7fffb7b68000-7fffb7b72000 r--p 0004f000 08:06 44044129                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/aggregations.cpython-37m-x86_64-linux-gnu.so
7fffb7b72000-7fffb7b73000 ---p 00059000 08:06 44044129                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/aggregations.cpython-37m-x86_64-linux-gnu.so
7fffb7b73000-7fffb7b74000 r--p 00059000 08:06 44044129                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/aggregations.cpython-37m-x86_64-linux-gnu.so
7fffb7b74000-7fffb7b79000 rw-p 0005a000 08:06 44044129                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/window/aggregations.cpython-37m-x86_64-linux-gnu.so
7fffb7b79000-7fffb7bfa000 rw-p 00000000 00:00 0
7fffb7bfa000-7fffb7bfc000 r--p 00000000 08:06 44567078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/mmap.cpython-37m-x86_64-linux-gnu.so
7fffb7bfc000-7fffb7bff000 r-xp 00002000 08:06 44567078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/mmap.cpython-37m-x86_64-linux-gnu.so
7fffb7bff000-7fffb7c00000 r--p 00005000 08:06 44567078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/mmap.cpython-37m-x86_64-linux-gnu.so
7fffb7c00000-7fffb7c01000 ---p 00006000 08:06 44567078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/mmap.cpython-37m-x86_64-linux-gnu.so
7fffb7c01000-7fffb7c02000 r--p 00006000 08:06 44567078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/mmap.cpython-37m-x86_64-linux-gnu.so
7fffb7c02000-7fffb7c03000 rw-p 00007000 08:06 44567078                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/mmap.cpython-37m-x86_64-linux-gnu.so
7fffb7c03000-7fffb7c43000 rw-p 00000000 00:00 0
7fffb7c43000-7fffb7c4a000 r--p 00000000 08:06 44044120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/internals.cpython-37m-x86_64-linux-gnu.so
7fffb7c4a000-7fffb7c70000 r-xp 00007000 08:06 44044120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/internals.cpython-37m-x86_64-linux-gnu.so
7fffb7c70000-7fffb7c78000 r--p 0002d000 08:06 44044120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/internals.cpython-37m-x86_64-linux-gnu.so
7fffb7c78000-7fffb7c79000 ---p 00035000 08:06 44044120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/internals.cpython-37m-x86_64-linux-gnu.so
7fffb7c79000-7fffb7c7a000 r--p 00035000 08:06 44044120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/internals.cpython-37m-x86_64-linux-gnu.so
7fffb7c7a000-7fffb7c7e000 rw-p 00036000 08:06 44044120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/internals.cpython-37m-x86_64-linux-gnu.so
7fffb7c7e000-7fffb7c7f000 rw-p 00000000 00:00 0
7fffb7c7f000-7fffb7c86000 r--p 00000000 08:06 44044113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/writers.cpython-37m-x86_64-linux-gnu.so
7fffb7c86000-7fffb7ca3000 r-xp 00007000 08:06 44044113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/writers.cpython-37m-x86_64-linux-gnu.so
7fffb7ca3000-7fffb7caa000 r--p 00024000 08:06 44044113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/writers.cpython-37m-x86_64-linux-gnu.so
7fffb7caa000-7fffb7cab000 ---p 0002b000 08:06 44044113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/writers.cpython-37m-x86_64-linux-gnu.so
7fffb7cab000-7fffb7cac000 r--p 0002b000 08:06 44044113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/writers.cpython-37m-x86_64-linux-gnu.so
7fffb7cac000-7fffb7cb0000 rw-p 0002c000 08:06 44044113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/writers.cpython-37m-x86_64-linux-gnu.so
7fffb7cb0000-7fffb7cf0000 rw-p 00000000 00:00 0
7fffb7cf0000-7fffb7cf3000 r--p 00000000 08:06 44044099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/indexing.cpython-37m-x86_64-linux-gnu.so
7fffb7cf3000-7fffb7cf8000 r-xp 00003000 08:06 44044099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/indexing.cpython-37m-x86_64-linux-gnu.so
7fffb7cf8000-7fffb7cfa000 r--p 00008000 08:06 44044099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/indexing.cpython-37m-x86_64-linux-gnu.so
7fffb7cfa000-7fffb7cfb000 r--p 00009000 08:06 44044099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/indexing.cpython-37m-x86_64-linux-gnu.so
7fffb7cfb000-7fffb7cfc000 rw-p 0000a000 08:06 44044099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/indexing.cpython-37m-x86_64-linux-gnu.so
7fffb7cfc000-7fffb7d7c000 rw-p 00000000 00:00 0
7fffb7d7c000-7fffb7d84000 r--p 00000000 08:06 44044122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reshape.cpython-37m-x86_64-linux-gnu.so
7fffb7d84000-7fffb7daf000 r-xp 00008000 08:06 44044122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reshape.cpython-37m-x86_64-linux-gnu.so
7fffb7daf000-7fffb7db7000 r--p 00033000 08:06 44044122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reshape.cpython-37m-x86_64-linux-gnu.so
7fffb7db7000-7fffb7db8000 ---p 0003b000 08:06 44044122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reshape.cpython-37m-x86_64-linux-gnu.so
7fffb7db8000-7fffb7db9000 r--p 0003b000 08:06 44044122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reshape.cpython-37m-x86_64-linux-gnu.so
7fffb7db9000-7fffb7dbd000 rw-p 0003c000 08:06 44044122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/reshape.cpython-37m-x86_64-linux-gnu.so
7fffb7dbd000-7fffb7f7e000 rw-p 00000000 00:00 0
7fffb7f7e000-7fffb7f88000 r--p 00000000 08:06 44044136                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/sparse.cpython-37m-x86_64-linux-gnu.so
7fffb7f88000-7fffb803b000 r-xp 0000a000 08:06 44044136                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/sparse.cpython-37m-x86_64-linux-gnu.so
7fffb803b000-7fffb8049000 r--p 000bd000 08:06 44044136                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/sparse.cpython-37m-x86_64-linux-gnu.so
7fffb8049000-7fffb804a000 r--p 000ca000 08:06 44044136                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/sparse.cpython-37m-x86_64-linux-gnu.so
7fffb804a000-7fffb8050000 rw-p 000cb000 08:06 44044136                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/sparse.cpython-37m-x86_64-linux-gnu.so
7fffb8050000-7fffb8052000 rw-p 00000000 00:00 0
7fffb8052000-7fffb8061000 r--p 00000000 08:06 44044141                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/join.cpython-37m-x86_64-linux-gnu.so
7fffb8061000-7fffb82b1000 r-xp 0000f000 08:06 44044141                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/join.cpython-37m-x86_64-linux-gnu.so
7fffb82b1000-7fffb82cd000 r--p 0025f000 08:06 44044141                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/join.cpython-37m-x86_64-linux-gnu.so
7fffb82cd000-7fffb82ce000 ---p 0027b000 08:06 44044141                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/join.cpython-37m-x86_64-linux-gnu.so
7fffb82ce000-7fffb82cf000 r--p 0027b000 08:06 44044141                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/join.cpython-37m-x86_64-linux-gnu.so
7fffb82cf000-7fffb82d7000 rw-p 0027c000 08:06 44044141                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/join.cpython-37m-x86_64-linux-gnu.so
7fffb82d7000-7fffb831b000 rw-p 00000000 00:00 0
7fffb831b000-7fffb8325000 r--p 00000000 08:06 44044133                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/index.cpython-37m-x86_64-linux-gnu.so
7fffb8325000-7fffb8387000 r-xp 0000a000 08:06 44044133                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/index.cpython-37m-x86_64-linux-gnu.so
7fffb8387000-7fffb8395000 r--p 0006c000 08:06 44044133                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/index.cpython-37m-x86_64-linux-gnu.so
7fffb8395000-7fffb8396000 r--p 00079000 08:06 44044133                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/index.cpython-37m-x86_64-linux-gnu.so
7fffb8396000-7fffb839e000 rw-p 0007a000 08:06 44044133                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/index.cpython-37m-x86_64-linux-gnu.so
7fffb839e000-7fffb845f000 rw-p 00000000 00:00 0
7fffb845f000-7fffb8461000 r--p 00000000 08:06 44567089                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_csv.cpython-37m-x86_64-linux-gnu.so
7fffb8461000-7fffb8465000 r-xp 00002000 08:06 44567089                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_csv.cpython-37m-x86_64-linux-gnu.so
7fffb8465000-7fffb8466000 r--p 00006000 08:06 44567089                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_csv.cpython-37m-x86_64-linux-gnu.so
7fffb8466000-7fffb8467000 ---p 00007000 08:06 44567089                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_csv.cpython-37m-x86_64-linux-gnu.so
7fffb8467000-7fffb8468000 r--p 00007000 08:06 44567089                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_csv.cpython-37m-x86_64-linux-gnu.so
7fffb8468000-7fffb846a000 rw-p 00008000 08:06 44567089                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_csv.cpython-37m-x86_64-linux-gnu.so
7fffb846a000-7fffb84ea000 rw-p 00000000 00:00 0
7fffb84ea000-7fffb84f1000 r--p 00000000 08:06 44044117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops.cpython-37m-x86_64-linux-gnu.so
7fffb84f1000-7fffb8512000 r-xp 00007000 08:06 44044117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops.cpython-37m-x86_64-linux-gnu.so
7fffb8512000-7fffb8519000 r--p 00028000 08:06 44044117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops.cpython-37m-x86_64-linux-gnu.so
7fffb8519000-7fffb851a000 r--p 0002e000 08:06 44044117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops.cpython-37m-x86_64-linux-gnu.so
7fffb851a000-7fffb851d000 rw-p 0002f000 08:06 44044117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops.cpython-37m-x86_64-linux-gnu.so
7fffb851d000-7fffb859e000 rw-p 00000000 00:00 0
7fffb859e000-7fffb85a4000 r--p 00000000 08:06 44044110                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashing.cpython-37m-x86_64-linux-gnu.so
7fffb85a4000-7fffb85bd000 r-xp 00006000 08:06 44044110                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashing.cpython-37m-x86_64-linux-gnu.so
7fffb85bd000-7fffb85c3000 r--p 0001f000 08:06 44044110                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashing.cpython-37m-x86_64-linux-gnu.so
7fffb85c3000-7fffb85c4000 ---p 00025000 08:06 44044110                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashing.cpython-37m-x86_64-linux-gnu.so
7fffb85c4000-7fffb85c5000 r--p 00025000 08:06 44044110                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashing.cpython-37m-x86_64-linux-gnu.so
7fffb85c5000-7fffb85c8000 rw-p 00026000 08:06 44044110                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashing.cpython-37m-x86_64-linux-gnu.so
7fffb85c8000-7fffb8648000 rw-p 00000000 00:00 0
7fffb8648000-7fffb864d000 r--p 00000000 08:06 44044109                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslib.cpython-37m-x86_64-linux-gnu.so
7fffb864d000-7fffb8668000 r-xp 00005000 08:06 44044109                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslib.cpython-37m-x86_64-linux-gnu.so
7fffb8668000-7fffb866c000 r--p 00020000 08:06 44044109                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslib.cpython-37m-x86_64-linux-gnu.so
7fffb866c000-7fffb866d000 ---p 00024000 08:06 44044109                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslib.cpython-37m-x86_64-linux-gnu.so
7fffb866d000-7fffb866e000 r--p 00024000 08:06 44044109                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslib.cpython-37m-x86_64-linux-gnu.so
7fffb866e000-7fffb8671000 rw-p 00025000 08:06 44044109                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslib.cpython-37m-x86_64-linux-gnu.so
7fffb8671000-7fffb867f000 r--p 00000000 08:06 44044134                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/lib.cpython-37m-x86_64-linux-gnu.so
7fffb867f000-7fffb86d5000 r-xp 0000e000 08:06 44044134                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/lib.cpython-37m-x86_64-linux-gnu.so
7fffb86d5000-7fffb86e5000 r--p 00064000 08:06 44044134                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/lib.cpython-37m-x86_64-linux-gnu.so
7fffb86e5000-7fffb86e6000 r--p 00073000 08:06 44044134                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/lib.cpython-37m-x86_64-linux-gnu.so
7fffb86e6000-7fffb86f2000 rw-p 00074000 08:06 44044134                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/lib.cpython-37m-x86_64-linux-gnu.so
7fffb86f2000-7fffb86f4000 rw-p 00000000 00:00 0
7fffb86f4000-7fffb8706000 r--p 00000000 08:06 44044140                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/algos.cpython-37m-x86_64-linux-gnu.so
7fffb8706000-7fffb883d000 r-xp 00012000 08:06 44044140                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/algos.cpython-37m-x86_64-linux-gnu.so
7fffb883d000-7fffb8858000 r--p 00149000 08:06 44044140                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/algos.cpython-37m-x86_64-linux-gnu.so
7fffb8858000-7fffb8859000 ---p 00164000 08:06 44044140                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/algos.cpython-37m-x86_64-linux-gnu.so
7fffb8859000-7fffb885a000 r--p 00164000 08:06 44044140                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/algos.cpython-37m-x86_64-linux-gnu.so
7fffb885a000-7fffb8865000 rw-p 00165000 08:06 44044140                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/algos.cpython-37m-x86_64-linux-gnu.so
7fffb8865000-7fffb8868000 rw-p 00000000 00:00 0
7fffb8868000-7fffb886c000 r--p 00000000 08:06 44044101                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops_dispatch.cpython-37m-x86_64-linux-gnu.so
7fffb886c000-7fffb8872000 r-xp 00004000 08:06 44044101                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops_dispatch.cpython-37m-x86_64-linux-gnu.so
7fffb8872000-7fffb8874000 r--p 0000a000 08:06 44044101                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops_dispatch.cpython-37m-x86_64-linux-gnu.so
7fffb8874000-7fffb8875000 r--p 0000b000 08:06 44044101                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops_dispatch.cpython-37m-x86_64-linux-gnu.so
7fffb8875000-7fffb8877000 rw-p 0000c000 08:06 44044101                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/ops_dispatch.cpython-37m-x86_64-linux-gnu.so
7fffb8877000-7fffb887d000 r--p 00000000 08:06 44044118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/vectorized.cpython-37m-x86_64-linux-gnu.so
7fffb887d000-7fffb88a0000 r-xp 00006000 08:06 44044118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/vectorized.cpython-37m-x86_64-linux-gnu.so
7fffb88a0000-7fffb88a7000 r--p 00029000 08:06 44044118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/vectorized.cpython-37m-x86_64-linux-gnu.so
7fffb88a7000-7fffb88a8000 ---p 00030000 08:06 44044118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/vectorized.cpython-37m-x86_64-linux-gnu.so
7fffb88a8000-7fffb88a9000 r--p 00030000 08:06 44044118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/vectorized.cpython-37m-x86_64-linux-gnu.so
7fffb88a9000-7fffb88ad000 rw-p 00031000 08:06 44044118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/vectorized.cpython-37m-x86_64-linux-gnu.so
7fffb88ad000-7fffb88ed000 rw-p 00000000 00:00 0
7fffb88ed000-7fffb88f8000 r--p 00000000 08:06 44044128                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/period.cpython-37m-x86_64-linux-gnu.so
7fffb88f8000-7fffb8933000 r-xp 0000b000 08:06 44044128                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/period.cpython-37m-x86_64-linux-gnu.so
7fffb8933000-7fffb8945000 r--p 00046000 08:06 44044128                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/period.cpython-37m-x86_64-linux-gnu.so
7fffb8945000-7fffb8946000 r--p 00057000 08:06 44044128                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/period.cpython-37m-x86_64-linux-gnu.so
7fffb8946000-7fffb894e000 rw-p 00058000 08:06 44044128                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/period.cpython-37m-x86_64-linux-gnu.so
7fffb894e000-7fffb8950000 rw-p 00000000 00:00 0
7fffb8950000-7fffb8954000 r--p 00000000 08:06 44044102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/properties.cpython-37m-x86_64-linux-gnu.so
7fffb8954000-7fffb895c000 r-xp 00004000 08:06 44044102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/properties.cpython-37m-x86_64-linux-gnu.so
7fffb895c000-7fffb895e000 r--p 0000c000 08:06 44044102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/properties.cpython-37m-x86_64-linux-gnu.so
7fffb895e000-7fffb895f000 ---p 0000e000 08:06 44044102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/properties.cpython-37m-x86_64-linux-gnu.so
7fffb895f000-7fffb8960000 r--p 0000e000 08:06 44044102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/properties.cpython-37m-x86_64-linux-gnu.so
7fffb8960000-7fffb8961000 rw-p 0000f000 08:06 44044102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/properties.cpython-37m-x86_64-linux-gnu.so
7fffb8961000-7fffb896c000 r--p 00000000 08:06 44044126                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/strptime.cpython-37m-x86_64-linux-gnu.so
7fffb896c000-7fffb89aa000 r-xp 0000b000 08:06 44044126                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/strptime.cpython-37m-x86_64-linux-gnu.so
7fffb89aa000-7fffb89b4000 r--p 00049000 08:06 44044126                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/strptime.cpython-37m-x86_64-linux-gnu.so
7fffb89b4000-7fffb89b5000 r--p 00052000 08:06 44044126                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/strptime.cpython-37m-x86_64-linux-gnu.so
7fffb89b5000-7fffb89bc000 rw-p 00053000 08:06 44044126                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/strptime.cpython-37m-x86_64-linux-gnu.so
7fffb89bc000-7fffb89fd000 rw-p 00000000 00:00 0
7fffb89fd000-7fffb8a05000 r--p 00000000 08:06 44044121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/fields.cpython-37m-x86_64-linux-gnu.so
7fffb8a05000-7fffb8a2e000 r-xp 00008000 08:06 44044121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/fields.cpython-37m-x86_64-linux-gnu.so
7fffb8a2e000-7fffb8a35000 r--p 00031000 08:06 44044121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/fields.cpython-37m-x86_64-linux-gnu.so
7fffb8a35000-7fffb8a36000 ---p 00038000 08:06 44044121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/fields.cpython-37m-x86_64-linux-gnu.so
7fffb8a36000-7fffb8a37000 r--p 00038000 08:06 44044121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/fields.cpython-37m-x86_64-linux-gnu.so
7fffb8a37000-7fffb8a3b000 rw-p 00039000 08:06 44044121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/fields.cpython-37m-x86_64-linux-gnu.so
7fffb8a3b000-7fffb8a3c000 rw-p 00000000 00:00 0
7fffb8a3c000-7fffb8a48000 r--p 00000000 08:06 44044131                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timestamps.cpython-37m-x86_64-linux-gnu.so
7fffb8a48000-7fffb8a8b000 r-xp 0000c000 08:06 44044131                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timestamps.cpython-37m-x86_64-linux-gnu.so
7fffb8a8b000-7fffb8a99000 r--p 0004f000 08:06 44044131                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timestamps.cpython-37m-x86_64-linux-gnu.so
7fffb8a99000-7fffb8a9a000 ---p 0005d000 08:06 44044131                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timestamps.cpython-37m-x86_64-linux-gnu.so
7fffb8a9a000-7fffb8a9b000 r--p 0005d000 08:06 44044131                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timestamps.cpython-37m-x86_64-linux-gnu.so
7fffb8a9b000-7fffb8aa4000 rw-p 0005e000 08:06 44044131                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timestamps.cpython-37m-x86_64-linux-gnu.so
7fffb8aa4000-7fffb8aa6000 rw-p 00000000 00:00 0
7fffb8aa6000-7fffb8ab1000 r--p 00000000 08:06 44044130                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timedeltas.cpython-37m-x86_64-linux-gnu.so
7fffb8ab1000-7fffb8af5000 r-xp 0000b000 08:06 44044130                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timedeltas.cpython-37m-x86_64-linux-gnu.so
7fffb8af5000-7fffb8b03000 r--p 0004f000 08:06 44044130                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timedeltas.cpython-37m-x86_64-linux-gnu.so
7fffb8b03000-7fffb8b04000 r--p 0005c000 08:06 44044130                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timedeltas.cpython-37m-x86_64-linux-gnu.so
7fffb8b04000-7fffb8b0b000 rw-p 0005d000 08:06 44044130                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timedeltas.cpython-37m-x86_64-linux-gnu.so
7fffb8b0b000-7fffb8b0d000 rw-p 00000000 00:00 0
7fffb8b0d000-7fffb8b20000 r--p 00000000 08:06 44044137                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/offsets.cpython-37m-x86_64-linux-gnu.so
7fffb8b20000-7fffb8bc9000 r-xp 00013000 08:06 44044137                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/offsets.cpython-37m-x86_64-linux-gnu.so
7fffb8bc9000-7fffb8be1000 r--p 000bc000 08:06 44044137                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/offsets.cpython-37m-x86_64-linux-gnu.so
7fffb8be1000-7fffb8be2000 ---p 000d4000 08:06 44044137                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/offsets.cpython-37m-x86_64-linux-gnu.so
7fffb8be2000-7fffb8be3000 r--p 000d4000 08:06 44044137                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/offsets.cpython-37m-x86_64-linux-gnu.so
7fffb8be3000-7fffb8bf5000 rw-p 000d5000 08:06 44044137                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/offsets.cpython-37m-x86_64-linux-gnu.so
7fffb8bf5000-7fffb8bf8000 rw-p 00000000 00:00 0
7fffb8bf8000-7fffb8c03000 r--p 00000000 08:06 44044127                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/parsing.cpython-37m-x86_64-linux-gnu.so
7fffb8c03000-7fffb8c45000 r-xp 0000b000 08:06 44044127                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/parsing.cpython-37m-x86_64-linux-gnu.so
7fffb8c45000-7fffb8c51000 r--p 0004d000 08:06 44044127                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/parsing.cpython-37m-x86_64-linux-gnu.so
7fffb8c51000-7fffb8c52000 r--p 00058000 08:06 44044127                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/parsing.cpython-37m-x86_64-linux-gnu.so
7fffb8c52000-7fffb8c59000 rw-p 00059000 08:06 44044127                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/parsing.cpython-37m-x86_64-linux-gnu.so
7fffb8c59000-7fffb8c5b000 rw-p 00000000 00:00 0
7fffb8c5b000-7fffb8c5f000 r--p 00000000 08:06 44044103                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/ccalendar.cpython-37m-x86_64-linux-gnu.so
7fffb8c5f000-7fffb8c69000 r-xp 00004000 08:06 44044103                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/ccalendar.cpython-37m-x86_64-linux-gnu.so
7fffb8c69000-7fffb8c6b000 r--p 0000e000 08:06 44044103                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/ccalendar.cpython-37m-x86_64-linux-gnu.so
7fffb8c6b000-7fffb8c6c000 ---p 00010000 08:06 44044103                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/ccalendar.cpython-37m-x86_64-linux-gnu.so
7fffb8c6c000-7fffb8c6d000 r--p 00010000 08:06 44044103                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/ccalendar.cpython-37m-x86_64-linux-gnu.so
7fffb8c6d000-7fffb8c6f000 rw-p 00011000 08:06 44044103                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/ccalendar.cpython-37m-x86_64-linux-gnu.so
7fffb8c6f000-7fffb8caf000 rw-p 00000000 00:00 0
7fffb8caf000-7fffb8cb6000 r--p 00000000 08:06 44044124                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/tzconversion.cpython-37m-x86_64-linux-gnu.so
7fffb8cb6000-7fffb8ce7000 r-xp 00007000 08:06 44044124                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/tzconversion.cpython-37m-x86_64-linux-gnu.so
7fffb8ce7000-7fffb8cef000 r--p 00038000 08:06 44044124                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/tzconversion.cpython-37m-x86_64-linux-gnu.so
7fffb8cef000-7fffb8cf0000 r--p 0003f000 08:06 44044124                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/tzconversion.cpython-37m-x86_64-linux-gnu.so
7fffb8cf0000-7fffb8cf4000 rw-p 00040000 08:06 44044124                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/tzconversion.cpython-37m-x86_64-linux-gnu.so
7fffb8cf4000-7fffb8cf5000 rw-p 00000000 00:00 0
7fffb8cf5000-7fffb8cfc000 r--p 00000000 08:06 44044119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timezones.cpython-37m-x86_64-linux-gnu.so
7fffb8cfc000-7fffb8d1e000 r-xp 00007000 08:06 44044119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timezones.cpython-37m-x86_64-linux-gnu.so
7fffb8d1e000-7fffb8d25000 r--p 00029000 08:06 44044119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timezones.cpython-37m-x86_64-linux-gnu.so
7fffb8d25000-7fffb8d26000 ---p 00030000 08:06 44044119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timezones.cpython-37m-x86_64-linux-gnu.so
7fffb8d26000-7fffb8d27000 r--p 00030000 08:06 44044119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timezones.cpython-37m-x86_64-linux-gnu.so
7fffb8d27000-7fffb8d2b000 rw-p 00031000 08:06 44044119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/timezones.cpython-37m-x86_64-linux-gnu.so
7fffb8d2b000-7fffb8d33000 r--p 00000000 08:06 44044112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/nattype.cpython-37m-x86_64-linux-gnu.so
7fffb8d33000-7fffb8d4b000 r-xp 00008000 08:06 44044112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/nattype.cpython-37m-x86_64-linux-gnu.so
7fffb8d4b000-7fffb8d53000 r--p 00020000 08:06 44044112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/nattype.cpython-37m-x86_64-linux-gnu.so
7fffb8d53000-7fffb8d54000 r--p 00027000 08:06 44044112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/nattype.cpython-37m-x86_64-linux-gnu.so
7fffb8d54000-7fffb8d58000 rw-p 00028000 08:06 44044112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/nattype.cpython-37m-x86_64-linux-gnu.so
7fffb8d58000-7fffb8d59000 rw-p 00000000 00:00 0
7fffb8d59000-7fffb8d61000 r--p 00000000 08:06 44044123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/conversion.cpython-37m-x86_64-linux-gnu.so
7fffb8d61000-7fffb8d8c000 r-xp 00008000 08:06 44044123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/conversion.cpython-37m-x86_64-linux-gnu.so
7fffb8d8c000-7fffb8d95000 r--p 00033000 08:06 44044123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/conversion.cpython-37m-x86_64-linux-gnu.so
7fffb8d95000-7fffb8d96000 r--p 0003b000 08:06 44044123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/conversion.cpython-37m-x86_64-linux-gnu.so
7fffb8d96000-7fffb8d9a000 rw-p 0003c000 08:06 44044123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/conversion.cpython-37m-x86_64-linux-gnu.so
7fffb8d9a000-7fffb8d9b000 rw-p 00000000 00:00 0
7fffb8d9b000-7fffb8da1000 r--p 00000000 08:06 44044106                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/dtypes.cpython-37m-x86_64-linux-gnu.so
7fffb8da1000-7fffb8db0000 r-xp 00006000 08:06 44044106                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/dtypes.cpython-37m-x86_64-linux-gnu.so
7fffb8db0000-7fffb8db4000 r--p 00015000 08:06 44044106                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/dtypes.cpython-37m-x86_64-linux-gnu.so
7fffb8db4000-7fffb8db5000 r--p 00018000 08:06 44044106                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/dtypes.cpython-37m-x86_64-linux-gnu.so
7fffb8db5000-7fffb8db8000 rw-p 00019000 08:06 44044106                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/dtypes.cpython-37m-x86_64-linux-gnu.so
7fffb8db8000-7fffb8db9000 rw-p 00000000 00:00 0
7fffb8db9000-7fffb8dc0000 r--p 00000000 08:06 44044111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/missing.cpython-37m-x86_64-linux-gnu.so
7fffb8dc0000-7fffb8ddb000 r-xp 00007000 08:06 44044111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/missing.cpython-37m-x86_64-linux-gnu.so
7fffb8ddb000-7fffb8de1000 r--p 00022000 08:06 44044111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/missing.cpython-37m-x86_64-linux-gnu.so
7fffb8de1000-7fffb8de2000 r--p 00027000 08:06 44044111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/missing.cpython-37m-x86_64-linux-gnu.so
7fffb8de2000-7fffb8de6000 rw-p 00028000 08:06 44044111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/missing.cpython-37m-x86_64-linux-gnu.so
7fffb8de6000-7fffb8de7000 rw-p 00000000 00:00 0
7fffb8de7000-7fffb8df2000 r--p 00000000 08:06 44044135                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashtable.cpython-37m-x86_64-linux-gnu.so
7fffb8df2000-7fffb8e59000 r-xp 0000b000 08:06 44044135                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashtable.cpython-37m-x86_64-linux-gnu.so
7fffb8e59000-7fffb8e68000 r--p 00072000 08:06 44044135                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashtable.cpython-37m-x86_64-linux-gnu.so
7fffb8e68000-7fffb8e69000 r--p 00080000 08:06 44044135                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashtable.cpython-37m-x86_64-linux-gnu.so
7fffb8e69000-7fffb8e74000 rw-p 00081000 08:06 44044135                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/hashtable.cpython-37m-x86_64-linux-gnu.so
7fffb8e74000-7fffb8eb6000 rw-p 00000000 00:00 0
7fffb8eb6000-7fffb8ec6000 r--p 00000000 08:06 44044139                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/interval.cpython-37m-x86_64-linux-gnu.so
7fffb8ec6000-7fffb8fb6000 r-xp 00010000 08:06 44044139                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/interval.cpython-37m-x86_64-linux-gnu.so
7fffb8fb6000-7fffb8fd7000 r--p 00100000 08:06 44044139                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/interval.cpython-37m-x86_64-linux-gnu.so
7fffb8fd7000-7fffb8fd8000 ---p 00121000 08:06 44044139                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/interval.cpython-37m-x86_64-linux-gnu.so
7fffb8fd8000-7fffb8fd9000 r--p 00121000 08:06 44044139                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/interval.cpython-37m-x86_64-linux-gnu.so
7fffb8fd9000-7fffb8fe5000 rw-p 00122000 08:06 44044139                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/interval.cpython-37m-x86_64-linux-gnu.so
7fffb8fe5000-7fffb9067000 rw-p 00000000 00:00 0
7fffb9067000-7fffb906a000 r--p 00000000 08:06 44438182                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/backends/_backend_agg.cpython-37m-x86_64-linux-gnu.so
7fffb906a000-7fffb9092000 r-xp 00003000 08:06 44438182                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/backends/_backend_agg.cpython-37m-x86_64-linux-gnu.so
7fffb9092000-7fffb9096000 r--p 0002b000 08:06 44438182                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/backends/_backend_agg.cpython-37m-x86_64-linux-gnu.so
7fffb9096000-7fffb9097000 ---p 0002f000 08:06 44438182                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/backends/_backend_agg.cpython-37m-x86_64-linux-gnu.so
7fffb9097000-7fffb9098000 r--p 0002f000 08:06 44438182                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/backends/_backend_agg.cpython-37m-x86_64-linux-gnu.so
7fffb9098000-7fffb9099000 rw-p 00030000 08:06 44438182                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/backends/_backend_agg.cpython-37m-x86_64-linux-gnu.so
7fffb9099000-7fffb965b000 rw-p 00000000 00:00 0
7fffb965b000-7fffb9665000 r--p 00000000 08:06 44303934                   /home/chadi/miniconda3/envs/ocp-covid/lib/libzstd.so.1.3.7
7fffb9665000-7fffb96f1000 r-xp 0000a000 08:06 44303934                   /home/chadi/miniconda3/envs/ocp-covid/lib/libzstd.so.1.3.7
7fffb96f1000-7fffb96fd000 r--p 00096000 08:06 44303934                   /home/chadi/miniconda3/envs/ocp-covid/lib/libzstd.so.1.3.7
7fffb96fd000-7fffb96fe000 ---p 000a2000 08:06 44303934                   /home/chadi/miniconda3/envs/ocp-covid/lib/libzstd.so.1.3.7
7fffb96fe000-7fffb96ff000 r--p 000a2000 08:06 44303934                   /home/chadi/miniconda3/envs/ocp-covid/lib/libzstd.so.1.3.7
7fffb96ff000-7fffb9700000 rw-p 000a3000 08:06 44303934                   /home/chadi/miniconda3/envs/ocp-covid/lib/libzstd.so.1.3.7
7fffb9700000-7fffb970a000 r--p 00000000 08:06 44565252                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiff.so.5.5.0
7fffb970a000-7fffb974f000 r-xp 0000a000 08:06 44565252                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiff.so.5.5.0
7fffb974f000-7fffb977b000 r--p 0004f000 08:06 44565252                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiff.so.5.5.0
7fffb977b000-7fffb977f000 r--p 0007a000 08:06 44565252                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiff.so.5.5.0
7fffb977f000-7fffb9780000 rw-p 0007e000 08:06 44565252                   /home/chadi/miniconda3/envs/ocp-covid/lib/libtiff.so.5.5.0
7fffb9780000-7fffb97bb000 r-xp 00000000 08:06 44041878                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjpeg.so.9.2.0
7fffb97bb000-7fffb99ba000 ---p 0003b000 08:06 44041878                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjpeg.so.9.2.0
7fffb99ba000-7fffb99bb000 r--p 0003a000 08:06 44041878                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjpeg.so.9.2.0
7fffb99bb000-7fffb99bc000 rw-p 0003b000 08:06 44041878                   /home/chadi/miniconda3/envs/ocp-covid/lib/libjpeg.so.9.2.0
7fffb99bc000-7fffb99d0000 r--p 00000000 08:06 43647992                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/PIL/_imaging.cpython-37m-x86_64-linux-gnu.so
7fffb99d0000-7fffb9a25000 r-xp 00014000 08:06 43647992                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/PIL/_imaging.cpython-37m-x86_64-linux-gnu.so
7fffb9a25000-7fffb9a34000 r--p 00069000 08:06 43647992                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/PIL/_imaging.cpython-37m-x86_64-linux-gnu.so
7fffb9a34000-7fffb9a35000 ---p 00078000 08:06 43647992                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/PIL/_imaging.cpython-37m-x86_64-linux-gnu.so
7fffb9a35000-7fffb9a39000 r--p 00078000 08:06 43647992                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/PIL/_imaging.cpython-37m-x86_64-linux-gnu.so
7fffb9a39000-7fffb9a3c000 rw-p 0007c000 08:06 43647992                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/PIL/_imaging.cpython-37m-x86_64-linux-gnu.so
7fffb9a3c000-7fffb9a3d000 rw-p 00000000 00:00 0
7fffb9a3d000-7fffb9a44000 r--p 00000000 08:06 44567117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/pyexpat.cpython-37m-x86_64-linux-gnu.so
7fffb9a44000-7fffb9a6f000 r-xp 00007000 08:06 44567117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/pyexpat.cpython-37m-x86_64-linux-gnu.so
7fffb9a6f000-7fffb9a79000 r--p 00032000 08:06 44567117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/pyexpat.cpython-37m-x86_64-linux-gnu.so
7fffb9a79000-7fffb9a7c000 r--p 0003b000 08:06 44567117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/pyexpat.cpython-37m-x86_64-linux-gnu.so
7fffb9a7c000-7fffb9a7e000 rw-p 0003e000 08:06 44567117                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/pyexpat.cpython-37m-x86_64-linux-gnu.so
7fffb9a7e000-7fffb9a82000 r--p 00000000 08:06 44567104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_elementtree.cpython-37m-x86_64-linux-gnu.so
7fffb9a82000-7fffb9a8a000 r-xp 00004000 08:06 44567104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_elementtree.cpython-37m-x86_64-linux-gnu.so
7fffb9a8a000-7fffb9a8d000 r--p 0000c000 08:06 44567104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_elementtree.cpython-37m-x86_64-linux-gnu.so
7fffb9a8d000-7fffb9a8e000 ---p 0000f000 08:06 44567104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_elementtree.cpython-37m-x86_64-linux-gnu.so
7fffb9a8e000-7fffb9a8f000 r--p 0000f000 08:06 44567104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_elementtree.cpython-37m-x86_64-linux-gnu.so
7fffb9a8f000-7fffb9a91000 rw-p 00010000 08:06 44567104                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_elementtree.cpython-37m-x86_64-linux-gnu.so
7fffb9a91000-7fffb9b11000 rw-p 00000000 00:00 0
7fffb9b11000-7fffb9b14000 r--p 00000000 08:06 44567118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/unicodedata.cpython-37m-x86_64-linux-gnu.so
7fffb9b14000-7fffb9b19000 r-xp 00003000 08:06 44567118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/unicodedata.cpython-37m-x86_64-linux-gnu.so
7fffb9b19000-7fffb9bfa000 r--p 00008000 08:06 44567118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/unicodedata.cpython-37m-x86_64-linux-gnu.so
7fffb9bfa000-7fffb9bfb000 r--p 000e8000 08:06 44567118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/unicodedata.cpython-37m-x86_64-linux-gnu.so
7fffb9bfb000-7fffb9c18000 rw-p 000e9000 08:06 44567118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/unicodedata.cpython-37m-x86_64-linux-gnu.so
7fffb9c18000-7fffb9fd8000 rw-p 00000000 00:00 0
7fffb9fd8000-7fffb9fdb000 r--p 00000000 08:06 44438181                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_path.cpython-37m-x86_64-linux-gnu.so
7fffb9fdb000-7fffb9ff4000 r-xp 00003000 08:06 44438181                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_path.cpython-37m-x86_64-linux-gnu.so
7fffb9ff4000-7fffb9ff8000 r--p 0001c000 08:06 44438181                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_path.cpython-37m-x86_64-linux-gnu.so
7fffb9ff8000-7fffb9ff9000 r--p 0001f000 08:06 44438181                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_path.cpython-37m-x86_64-linux-gnu.so
7fffb9ff9000-7fffb9ffa000 rw-p 00020000 08:06 44438181                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_path.cpython-37m-x86_64-linux-gnu.so
7fffb9ffa000-7fffba0ba000 rw-p 00000000 00:00 0
7fffba0ba000-7fffba0c2000 r--p 00000000 08:06 44442279                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/kiwisolver.cpython-37m-x86_64-linux-gnu.so
7fffba0c2000-7fffba0eb000 r-xp 00008000 08:06 44442279                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/kiwisolver.cpython-37m-x86_64-linux-gnu.so
7fffba0eb000-7fffba0f0000 r--p 00031000 08:06 44442279                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/kiwisolver.cpython-37m-x86_64-linux-gnu.so
7fffba0f0000-7fffba0f1000 ---p 00036000 08:06 44442279                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/kiwisolver.cpython-37m-x86_64-linux-gnu.so
7fffba0f1000-7fffba0f2000 r--p 00036000 08:06 44442279                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/kiwisolver.cpython-37m-x86_64-linux-gnu.so
7fffba0f2000-7fffba0f3000 rw-p 00037000 08:06 44442279                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/kiwisolver.cpython-37m-x86_64-linux-gnu.so
7fffba0f3000-7fffba0f9000 r--p 00000000 08:06 44565536                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpng16.so.16.37.0
7fffba0f9000-7fffba11f000 r-xp 00006000 08:06 44565536                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpng16.so.16.37.0
7fffba11f000-7fffba12a000 r--p 0002c000 08:06 44565536                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpng16.so.16.37.0
7fffba12a000-7fffba12b000 r--p 00036000 08:06 44565536                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpng16.so.16.37.0
7fffba12b000-7fffba12c000 rw-p 00037000 08:06 44565536                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpng16.so.16.37.0
7fffba12c000-7fffba1ce000 r--p 00000000 08:06 43386894                   /home/chadi/miniconda3/envs/ocp-covid/lib/libstdc++.so.6.0.26
7fffba1ce000-7fffba24d000 r-xp 000a2000 08:06 43386894                   /home/chadi/miniconda3/envs/ocp-covid/lib/libstdc++.so.6.0.26
7fffba24d000-7fffba28e000 r--p 00121000 08:06 43386894                   /home/chadi/miniconda3/envs/ocp-covid/lib/libstdc++.so.6.0.26
7fffba28e000-7fffba299000 r--p 00161000 08:06 43386894                   /home/chadi/miniconda3/envs/ocp-covid/lib/libstdc++.so.6.0.26
7fffba299000-7fffba29d000 rw-p 0016c000 08:06 43386894                   /home/chadi/miniconda3/envs/ocp-covid/lib/libstdc++.so.6.0.26
7fffba29d000-7fffba2a0000 rw-p 00000000 00:00 0
7fffba2a0000-7fffba2ae000 r--p 00000000 08:06 44041757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreetype.so.6.17.4
7fffba2ae000-7fffba325000 r-xp 0000e000 08:06 44041757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreetype.so.6.17.4
7fffba325000-7fffba34d000 r--p 00085000 08:06 44041757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreetype.so.6.17.4
7fffba34d000-7fffba354000 r--p 000ac000 08:06 44041757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreetype.so.6.17.4
7fffba354000-7fffba355000 rw-p 000b3000 08:06 44041757                   /home/chadi/miniconda3/envs/ocp-covid/lib/libfreetype.so.6.17.4
7fffba355000-7fffba359000 r--p 00000000 08:06 44438179                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so
7fffba359000-7fffba362000 r-xp 00004000 08:06 44438179                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so
7fffba362000-7fffba367000 r--p 0000d000 08:06 44438179                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so
7fffba367000-7fffba368000 ---p 00012000 08:06 44438179                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so
7fffba368000-7fffba369000 r--p 00012000 08:06 44438179                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so
7fffba369000-7fffba36a000 rw-p 00013000 08:06 44438179                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so
7fffba36a000-7fffba36b000 rw-p 00000000 00:00 0
7fffba36b000-7fffba36d000 r--p 00000000 08:06 44567107                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_json.cpython-37m-x86_64-linux-gnu.so
7fffba36d000-7fffba37b000 r-xp 00002000 08:06 44567107                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_json.cpython-37m-x86_64-linux-gnu.so
7fffba37b000-7fffba37d000 r--p 00010000 08:06 44567107                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_json.cpython-37m-x86_64-linux-gnu.so
7fffba37d000-7fffba37e000 r--p 00011000 08:06 44567107                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_json.cpython-37m-x86_64-linux-gnu.so
7fffba37e000-7fffba37f000 rw-p 00012000 08:06 44567107                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_json.cpython-37m-x86_64-linux-gnu.so
7fffba37f000-7fffba57f000 rw-p 00000000 00:00 0
7fffba57f000-7fffba581000 r--p 00000000 08:06 44567083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/binascii.cpython-37m-x86_64-linux-gnu.so
7fffba581000-7fffba584000 r-xp 00002000 08:06 44567083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/binascii.cpython-37m-x86_64-linux-gnu.so
7fffba584000-7fffba586000 r--p 00005000 08:06 44567083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/binascii.cpython-37m-x86_64-linux-gnu.so
7fffba586000-7fffba587000 r--p 00006000 08:06 44567083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/binascii.cpython-37m-x86_64-linux-gnu.so
7fffba587000-7fffba588000 rw-p 00007000 08:06 44567083                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/binascii.cpython-37m-x86_64-linux-gnu.so
7fffba588000-7fffba648000 rw-p 00000000 00:00 0
7fffba648000-7fffba64a000 r--p 00000000 08:06 44567085                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/select.cpython-37m-x86_64-linux-gnu.so
7fffba64a000-7fffba64e000 r-xp 00002000 08:06 44567085                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/select.cpython-37m-x86_64-linux-gnu.so
7fffba64e000-7fffba64f000 r--p 00006000 08:06 44567085                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/select.cpython-37m-x86_64-linux-gnu.so
7fffba64f000-7fffba650000 r--p 00006000 08:06 44567085                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/select.cpython-37m-x86_64-linux-gnu.so
7fffba650000-7fffba652000 rw-p 00007000 08:06 44567085                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/select.cpython-37m-x86_64-linux-gnu.so
7fffba652000-7fffba654000 r--p 00000000 08:06 44567073                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_posixsubprocess.cpython-37m-x86_64-linux-gnu.so
7fffba654000-7fffba656000 r-xp 00002000 08:06 44567073                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_posixsubprocess.cpython-37m-x86_64-linux-gnu.so
7fffba656000-7fffba657000 r--p 00004000 08:06 44567073                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_posixsubprocess.cpython-37m-x86_64-linux-gnu.so
7fffba657000-7fffba658000 r--p 00004000 08:06 44567073                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_posixsubprocess.cpython-37m-x86_64-linux-gnu.so
7fffba658000-7fffba659000 rw-p 00005000 08:06 44567073                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_posixsubprocess.cpython-37m-x86_64-linux-gnu.so
7fffba659000-7fffba6d9000 rw-p 00000000 00:00 0
7fffba6d9000-7fffba6db000 r--p 00000000 08:06 44567068                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_random.cpython-37m-x86_64-linux-gnu.so
7fffba6db000-7fffba6dd000 r-xp 00002000 08:06 44567068                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_random.cpython-37m-x86_64-linux-gnu.so
7fffba6dd000-7fffba6de000 r--p 00004000 08:06 44567068                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_random.cpython-37m-x86_64-linux-gnu.so
7fffba6de000-7fffba6df000 r--p 00004000 08:06 44567068                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_random.cpython-37m-x86_64-linux-gnu.so
7fffba6df000-7fffba6e0000 rw-p 00005000 08:06 44567068                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_random.cpython-37m-x86_64-linux-gnu.so
7fffba6e0000-7fffba75b000 r--p 00000000 08:06 44181950                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcrypto.so.1.1
7fffba75b000-7fffba8ee000 r-xp 0007b000 08:06 44181950                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcrypto.so.1.1
7fffba8ee000-7fffba97b000 r--p 0020e000 08:06 44181950                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcrypto.so.1.1
7fffba97b000-7fffba9a6000 r--p 0029a000 08:06 44181950                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcrypto.so.1.1
7fffba9a6000-7fffba9a8000 rw-p 002c5000 08:06 44181950                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcrypto.so.1.1
7fffba9a8000-7fffba9ec000 rw-p 00000000 00:00 0
7fffba9ef000-7fffba9f2000 r--p 00000000 08:06 44044100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/np_datetime.cpython-37m-x86_64-linux-gnu.so
7fffba9f2000-7fffba9f8000 r-xp 00003000 08:06 44044100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/np_datetime.cpython-37m-x86_64-linux-gnu.so
7fffba9f8000-7fffba9fa000 r--p 00009000 08:06 44044100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/np_datetime.cpython-37m-x86_64-linux-gnu.so
7fffba9fa000-7fffba9fb000 ---p 0000b000 08:06 44044100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/np_datetime.cpython-37m-x86_64-linux-gnu.so
7fffba9fb000-7fffba9fc000 r--p 0000b000 08:06 44044100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/np_datetime.cpython-37m-x86_64-linux-gnu.so
7fffba9fc000-7fffba9fd000 rw-p 0000c000 08:06 44044100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/np_datetime.cpython-37m-x86_64-linux-gnu.so
7fffba9fd000-7fffbaa00000 r--p 00000000 08:06 44044098                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/base.cpython-37m-x86_64-linux-gnu.so
7fffbaa00000-7fffbaa05000 r-xp 00003000 08:06 44044098                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/base.cpython-37m-x86_64-linux-gnu.so
7fffbaa05000-7fffbaa07000 r--p 00008000 08:06 44044098                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/base.cpython-37m-x86_64-linux-gnu.so
7fffbaa07000-7fffbaa08000 r--p 00009000 08:06 44044098                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/base.cpython-37m-x86_64-linux-gnu.so
7fffbaa08000-7fffbaa09000 rw-p 0000a000 08:06 44044098                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pandas/_libs/tslibs/base.cpython-37m-x86_64-linux-gnu.so
7fffbaa09000-7fffbaa0b000 r--p 00000000 08:06 44438180                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_image.cpython-37m-x86_64-linux-gnu.so
7fffbaa0b000-7fffbaa25000 r-xp 00002000 08:06 44438180                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_image.cpython-37m-x86_64-linux-gnu.so
7fffbaa25000-7fffbaa28000 r--p 0001c000 08:06 44438180                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_image.cpython-37m-x86_64-linux-gnu.so
7fffbaa28000-7fffbaa29000 r--p 0001e000 08:06 44438180                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_image.cpython-37m-x86_64-linux-gnu.so
7fffbaa29000-7fffbaa2a000 rw-p 0001f000 08:06 44438180                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/matplotlib/_image.cpython-37m-x86_64-linux-gnu.so
7fffbaa2a000-7fffbaa2e000 r--p 00000000 08:06 44567100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/array.cpython-37m-x86_64-linux-gnu.so
7fffbaa2e000-7fffbaa35000 r-xp 00004000 08:06 44567100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/array.cpython-37m-x86_64-linux-gnu.so
7fffbaa35000-7fffbaa38000 r--p 0000b000 08:06 44567100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/array.cpython-37m-x86_64-linux-gnu.so
7fffbaa38000-7fffbaa39000 ---p 0000e000 08:06 44567100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/array.cpython-37m-x86_64-linux-gnu.so
7fffbaa39000-7fffbaa3a000 r--p 0000e000 08:06 44567100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/array.cpython-37m-x86_64-linux-gnu.so
7fffbaa3a000-7fffbaa3d000 rw-p 0000f000 08:06 44567100                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/array.cpython-37m-x86_64-linux-gnu.so
7fffbaa3d000-7fffbabbd000 rw-p 00000000 00:00 0
7fffbabbe000-7fffbabbf000 r--p 00000000 08:06 44567054                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_opcode.cpython-37m-x86_64-linux-gnu.so
7fffbabbf000-7fffbabc0000 r-xp 00001000 08:06 44567054                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_opcode.cpython-37m-x86_64-linux-gnu.so
7fffbabc0000-7fffbabc1000 r--p 00002000 08:06 44567054                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_opcode.cpython-37m-x86_64-linux-gnu.so
7fffbabc1000-7fffbabc2000 r--p 00002000 08:06 44567054                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_opcode.cpython-37m-x86_64-linux-gnu.so
7fffbabc2000-7fffbabc3000 rw-p 00003000 08:06 44567054                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_opcode.cpython-37m-x86_64-linux-gnu.so
7fffbabc3000-7fffbabc4000 r--p 00000000 08:06 44567057                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bisect.cpython-37m-x86_64-linux-gnu.so
7fffbabc4000-7fffbabc5000 r-xp 00001000 08:06 44567057                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bisect.cpython-37m-x86_64-linux-gnu.so
7fffbabc5000-7fffbabc6000 r--p 00002000 08:06 44567057                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bisect.cpython-37m-x86_64-linux-gnu.so
7fffbabc6000-7fffbabc7000 r--p 00002000 08:06 44567057                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bisect.cpython-37m-x86_64-linux-gnu.so
7fffbabc7000-7fffbabc8000 rw-p 00003000 08:06 44567057                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bisect.cpython-37m-x86_64-linux-gnu.so
7fffbabc8000-7fffbabcb000 r--p 00000000 08:06 44567114                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_sha3.cpython-37m-x86_64-linux-gnu.so
7fffbabcb000-7fffbabde000 r-xp 00003000 08:06 44567114                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_sha3.cpython-37m-x86_64-linux-gnu.so
7fffbabde000-7fffbabdf000 r--p 00016000 08:06 44567114                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_sha3.cpython-37m-x86_64-linux-gnu.so
7fffbabdf000-7fffbabe0000 ---p 00017000 08:06 44567114                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_sha3.cpython-37m-x86_64-linux-gnu.so
7fffbabe0000-7fffbabe1000 r--p 00017000 08:06 44567114                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_sha3.cpython-37m-x86_64-linux-gnu.so
7fffbabe1000-7fffbabe3000 rw-p 00018000 08:06 44567114                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_sha3.cpython-37m-x86_64-linux-gnu.so
7fffbabe3000-7fffbabe5000 r--p 00000000 08:06 44567102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_blake2.cpython-37m-x86_64-linux-gnu.so
7fffbabe5000-7fffbabed000 r-xp 00002000 08:06 44567102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_blake2.cpython-37m-x86_64-linux-gnu.so
7fffbabed000-7fffbabee000 r--p 0000a000 08:06 44567102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_blake2.cpython-37m-x86_64-linux-gnu.so
7fffbabee000-7fffbabef000 ---p 0000b000 08:06 44567102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_blake2.cpython-37m-x86_64-linux-gnu.so
7fffbabef000-7fffbabf0000 r--p 0000b000 08:06 44567102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_blake2.cpython-37m-x86_64-linux-gnu.so
7fffbabf0000-7fffbabf1000 rw-p 0000c000 08:06 44567102                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_blake2.cpython-37m-x86_64-linux-gnu.so
7fffbabf1000-7fffbabfc000 r--p 00000000 08:06 44439122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/random/mtrand.cpython-37m-x86_64-linux-gnu.so
7fffbabfc000-7fffbac81000 r-xp 0000b000 08:06 44439122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/random/mtrand.cpython-37m-x86_64-linux-gnu.so
7fffbac81000-7fffbaca4000 r--p 00090000 08:06 44439122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/random/mtrand.cpython-37m-x86_64-linux-gnu.so
7fffbaca4000-7fffbaca5000 ---p 000b3000 08:06 44439122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/random/mtrand.cpython-37m-x86_64-linux-gnu.so
7fffbaca5000-7fffbaca6000 r--p 000b3000 08:06 44439122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/random/mtrand.cpython-37m-x86_64-linux-gnu.so
7fffbaca6000-7fffbacca000 rw-p 000b4000 08:06 44439122                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/random/mtrand.cpython-37m-x86_64-linux-gnu.so
7fffbacca000-7fffbad0d000 rw-p 00000000 00:00 0
7fffbad0d000-7fffbad0e000 r--p 00000000 08:06 44439118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/fft/fftpack_lite.cpython-37m-x86_64-linux-gnu.so
7fffbad0e000-7fffbad16000 r-xp 00001000 08:06 44439118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/fft/fftpack_lite.cpython-37m-x86_64-linux-gnu.so
7fffbad16000-7fffbad17000 r--p 00009000 08:06 44439118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/fft/fftpack_lite.cpython-37m-x86_64-linux-gnu.so
7fffbad17000-7fffbad18000 r--p 00009000 08:06 44439118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/fft/fftpack_lite.cpython-37m-x86_64-linux-gnu.so
7fffbad18000-7fffbad19000 rw-p 0000a000 08:06 44439118                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/fft/fftpack_lite.cpython-37m-x86_64-linux-gnu.so
7fffbad19000-7fffbad59000 rw-p 00000000 00:00 0
7fffbad59000-7fffbad60000 r--p 00000000 08:06 44567119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_decimal.cpython-37m-x86_64-linux-gnu.so
7fffbad60000-7fffbad95000 r-xp 00007000 08:06 44567119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_decimal.cpython-37m-x86_64-linux-gnu.so
7fffbad95000-7fffbad9f000 r--p 0003c000 08:06 44567119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_decimal.cpython-37m-x86_64-linux-gnu.so
7fffbad9f000-7fffbada0000 ---p 00046000 08:06 44567119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_decimal.cpython-37m-x86_64-linux-gnu.so
7fffbada0000-7fffbada1000 r--p 00046000 08:06 44567119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_decimal.cpython-37m-x86_64-linux-gnu.so
7fffbada1000-7fffbada9000 rw-p 00047000 08:06 44567119                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_decimal.cpython-37m-x86_64-linux-gnu.so
7fffbada9000-7fffbadab000 r--p 00000000 08:06 44567062                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/grp.cpython-37m-x86_64-linux-gnu.so
7fffbadab000-7fffbadac000 r-xp 00002000 08:06 44567062                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/grp.cpython-37m-x86_64-linux-gnu.so
7fffbadac000-7fffbadad000 r--p 00003000 08:06 44567062                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/grp.cpython-37m-x86_64-linux-gnu.so
7fffbadad000-7fffbadae000 r--p 00003000 08:06 44567062                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/grp.cpython-37m-x86_64-linux-gnu.so
7fffbadae000-7fffbadaf000 rw-p 00004000 08:06 44567062                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/grp.cpython-37m-x86_64-linux-gnu.so
7fffbadaf000-7fffbadb3000 r--p 00000000 08:06 43391082                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblzma.so.5.2.5
7fffbadb3000-7fffbadca000 r-xp 00004000 08:06 43391082                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblzma.so.5.2.5
7fffbadca000-7fffbadd5000 r--p 0001b000 08:06 43391082                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblzma.so.5.2.5
7fffbadd5000-7fffbadd6000 ---p 00026000 08:06 43391082                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblzma.so.5.2.5
7fffbadd6000-7fffbadd7000 r--p 00026000 08:06 43391082                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblzma.so.5.2.5
7fffbadd7000-7fffbadd8000 rw-p 00027000 08:06 43391082                   /home/chadi/miniconda3/envs/ocp-covid/lib/liblzma.so.5.2.5
7fffbadd8000-7fffbaddb000 r--p 00000000 08:06 44567087                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_lzma.cpython-37m-x86_64-linux-gnu.so
7fffbaddb000-7fffbaddf000 r-xp 00003000 08:06 44567087                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_lzma.cpython-37m-x86_64-linux-gnu.so
7fffbaddf000-7fffbade1000 r--p 00007000 08:06 44567087                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_lzma.cpython-37m-x86_64-linux-gnu.so
7fffbade1000-7fffbade2000 r--p 00008000 08:06 44567087                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_lzma.cpython-37m-x86_64-linux-gnu.so
7fffbade2000-7fffbade4000 rw-p 00009000 08:06 44567087                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_lzma.cpython-37m-x86_64-linux-gnu.so
7fffbade4000-7fffbade7000 r--p 00000000 08:06 44567091                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bz2.cpython-37m-x86_64-linux-gnu.so
7fffbade7000-7fffbadf6000 r-xp 00003000 08:06 44567091                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bz2.cpython-37m-x86_64-linux-gnu.so
7fffbadf6000-7fffbadf8000 r--p 00012000 08:06 44567091                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bz2.cpython-37m-x86_64-linux-gnu.so
7fffbadf8000-7fffbadf9000 ---p 00014000 08:06 44567091                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bz2.cpython-37m-x86_64-linux-gnu.so
7fffbadf9000-7fffbadfa000 r--p 00014000 08:06 44567091                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bz2.cpython-37m-x86_64-linux-gnu.so
7fffbadfa000-7fffbadfc000 rw-p 00015000 08:06 44567091                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_bz2.cpython-37m-x86_64-linux-gnu.so
7fffbadfc000-7fffbae3c000 rw-p 00000000 00:00 0
7fffbae3c000-7fffbae3f000 r--p 00000000 08:06 43386567                   /home/chadi/miniconda3/envs/ocp-covid/lib/libz.so.1.2.11
7fffbae3f000-7fffbae53000 r-xp 00003000 08:06 43386567                   /home/chadi/miniconda3/envs/ocp-covid/lib/libz.so.1.2.11
7fffbae53000-7fffbae5a000 r--p 00017000 08:06 43386567                   /home/chadi/miniconda3/envs/ocp-covid/lib/libz.so.1.2.11
7fffbae5a000-7fffbae5b000 r--p 0001d000 08:06 43386567                   /home/chadi/miniconda3/envs/ocp-covid/lib/libz.so.1.2.11
7fffbae5b000-7fffbae5c000 rw-p 0001e000 08:06 43386567                   /home/chadi/miniconda3/envs/ocp-covid/lib/libz.so.1.2.11
7fffbae5c000-7fffbae5e000 r--p 00000000 08:06 44567086                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/zlib.cpython-37m-x86_64-linux-gnu.so
7fffbae5e000-7fffbae62000 r-xp 00002000 08:06 44567086                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/zlib.cpython-37m-x86_64-linux-gnu.so
7fffbae62000-7fffbae63000 r--p 00006000 08:06 44567086                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/zlib.cpython-37m-x86_64-linux-gnu.so
7fffbae63000-7fffbae64000 ---p 00007000 08:06 44567086                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/zlib.cpython-37m-x86_64-linux-gnu.so
7fffbae64000-7fffbae65000 r--p 00007000 08:06 44567086                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/zlib.cpython-37m-x86_64-linux-gnu.so
7fffbae65000-7fffbae67000 rw-p 00008000 08:06 44567086                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/zlib.cpython-37m-x86_64-linux-gnu.so
7fffbae67000-7fffbaf27000 rw-p 00000000 00:00 0
7fffbaf27000-7fffbaf2f000 r--p 00000000 08:06 44439121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/_umath_linalg.cpython-37m-x86_64-linux-gnu.so
7fffbaf2f000-7fffbaf4a000 r-xp 00008000 08:06 44439121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/_umath_linalg.cpython-37m-x86_64-linux-gnu.so
7fffbaf4a000-7fffbaf4f000 r--p 00023000 08:06 44439121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/_umath_linalg.cpython-37m-x86_64-linux-gnu.so
7fffbaf4f000-7fffbaf50000 ---p 00028000 08:06 44439121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/_umath_linalg.cpython-37m-x86_64-linux-gnu.so
7fffbaf50000-7fffbaf51000 r--p 00028000 08:06 44439121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/_umath_linalg.cpython-37m-x86_64-linux-gnu.so
7fffbaf51000-7fffbaf52000 rw-p 00029000 08:06 44439121                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/_umath_linalg.cpython-37m-x86_64-linux-gnu.so
7fffbaf52000-7fffbafd2000 rw-p 00000000 00:00 0
7fffbafd2000-7fffbafdb000 r--p 00000000 08:06 44439120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_tests.cpython-37m-x86_64-linux-gnu.so
7fffbafdb000-7fffbafed000 r-xp 00009000 08:06 44439120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_tests.cpython-37m-x86_64-linux-gnu.so
7fffbafed000-7fffbaff2000 r--p 0001b000 08:06 44439120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_tests.cpython-37m-x86_64-linux-gnu.so
7fffbaff2000-7fffbaff3000 r--p 0001f000 08:06 44439120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_tests.cpython-37m-x86_64-linux-gnu.so
7fffbaff3000-7fffbaff4000 rw-p 00020000 08:06 44439120                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_tests.cpython-37m-x86_64-linux-gnu.so
7fffbaff4000-7fffbb074000 rw-p 00000000 00:00 0
7fffbb074000-7fffbb079000 r--p 00000000 08:06 44567113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_pickle.cpython-37m-x86_64-linux-gnu.so
7fffbb079000-7fffbb08c000 r-xp 00005000 08:06 44567113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_pickle.cpython-37m-x86_64-linux-gnu.so
7fffbb08c000-7fffbb090000 r--p 00018000 08:06 44567113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_pickle.cpython-37m-x86_64-linux-gnu.so
7fffbb090000-7fffbb091000 r--p 0001b000 08:06 44567113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_pickle.cpython-37m-x86_64-linux-gnu.so
7fffbb091000-7fffbb094000 rw-p 0001c000 08:06 44567113                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_pickle.cpython-37m-x86_64-linux-gnu.so
7fffbb094000-7fffbb0d4000 rw-p 00000000 00:00 0
7fffbb0d4000-7fffbb0d7000 r--p 00000000 08:06 44567099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_struct.cpython-37m-x86_64-linux-gnu.so
7fffbb0d7000-7fffbb0dd000 r-xp 00003000 08:06 44567099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_struct.cpython-37m-x86_64-linux-gnu.so
7fffbb0dd000-7fffbb0df000 r--p 00009000 08:06 44567099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_struct.cpython-37m-x86_64-linux-gnu.so
7fffbb0df000-7fffbb0e0000 ---p 0000b000 08:06 44567099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_struct.cpython-37m-x86_64-linux-gnu.so
7fffbb0e0000-7fffbb0e1000 r--p 0000b000 08:06 44567099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_struct.cpython-37m-x86_64-linux-gnu.so
7fffbb0e1000-7fffbb0e3000 rw-p 0000c000 08:06 44567099                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_struct.cpython-37m-x86_64-linux-gnu.so
7fffbb0e3000-7fffbb0e5000 r--p 00000000 08:06 43386437                   /home/chadi/miniconda3/envs/ocp-covid/lib/libffi.so.7.1.0
7fffbb0e5000-7fffbb0eb000 r-xp 00002000 08:06 43386437                   /home/chadi/miniconda3/envs/ocp-covid/lib/libffi.so.7.1.0
7fffbb0eb000-7fffbb0ec000 r--p 00008000 08:06 43386437                   /home/chadi/miniconda3/envs/ocp-covid/lib/libffi.so.7.1.0
7fffbb0ec000-7fffbb0ed000 ---p 00009000 08:06 43386437                   /home/chadi/miniconda3/envs/ocp-covid/lib/libffi.so.7.1.0
7fffbb0ed000-7fffbb0ee000 r--p 00009000 08:06 43386437                   /home/chadi/miniconda3/envs/ocp-covid/lib/libffi.so.7.1.0
7fffbb0ee000-7fffbb0ef000 rw-p 0000a000 08:06 43386437                   /home/chadi/miniconda3/envs/ocp-covid/lib/libffi.so.7.1.0
7fffbb0ef000-7fffbb0f1000 r--p 00000000 08:06 44567080                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_hashlib.cpython-37m-x86_64-linux-gnu.so
7fffbb0f1000-7fffbb0f4000 r-xp 00002000 08:06 44567080                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_hashlib.cpython-37m-x86_64-linux-gnu.so
7fffbb0f4000-7fffbb0f6000 r--p 00005000 08:06 44567080                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_hashlib.cpython-37m-x86_64-linux-gnu.so
7fffbb0f6000-7fffbb0f7000 r--p 00006000 08:06 44567080                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_hashlib.cpython-37m-x86_64-linux-gnu.so
7fffbb0f7000-7fffbb0f8000 rw-p 00007000 08:06 44567080                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_hashlib.cpython-37m-x86_64-linux-gnu.so
7fffbb0f8000-7fffbb238000 rw-p 00000000 00:00 0
7fffbb238000-7fffbb23d000 r--p 00000000 08:06 44567115                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_datetime.cpython-37m-x86_64-linux-gnu.so
7fffbb23d000-7fffbb24d000 r-xp 00005000 08:06 44567115                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_datetime.cpython-37m-x86_64-linux-gnu.so
7fffbb24d000-7fffbb253000 r--p 00015000 08:06 44567115                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_datetime.cpython-37m-x86_64-linux-gnu.so
7fffbb253000-7fffbb254000 r--p 0001a000 08:06 44567115                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_datetime.cpython-37m-x86_64-linux-gnu.so
7fffbb254000-7fffbb256000 rw-p 0001b000 08:06 44567115                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_datetime.cpython-37m-x86_64-linux-gnu.so
7fffbb256000-7fffbf256000 rw-p 00000000 00:00 0
7fffbf256000-7fffbf257000 ---p 00000000 00:00 0
7fffbf257000-7fffc1a57000 rw-p 00000000 00:00 0
7fffc1a57000-7fffc1a58000 ---p 00000000 00:00 0
7fffc1a58000-7fffc4258000 rw-p 00000000 00:00 0
7fffc4258000-7fffc4259000 ---p 00000000 00:00 0
7fffc4259000-7fffc6a59000 rw-p 00000000 00:00 0
7fffc6a59000-7fffc6a5a000 ---p 00000000 00:00 0
7fffc6a5a000-7fffc725a000 rw-p 00000000 00:00 0
7fffc725a000-7fffc725b000 ---p 00000000 00:00 0
7fffc725b000-7fffcba5b000 rw-p 00000000 00:00 0
7fffcba5b000-7fffcba5c000 ---p 00000000 00:00 0
7fffcba5c000-7fffce25c000 rw-p 00000000 00:00 0
7fffce25c000-7fffce25d000 ---p 00000000 00:00 0
7fffce25d000-7fffd0a5d000 rw-p 00000000 00:00 0
7fffd10f2000-7fffd5a5f000 rw-p 00000000 00:00 0
7fffd60f4000-7fffd8a61000 rw-p 00000000 00:00 0
7fffd8a61000-7fffd8a62000 ---p 00000000 00:00 0
7fffd8a62000-7fffdb262000 rw-p 00000000 00:00 0
7fffdb262000-7fffdb263000 ---p 00000000 00:00 0
7fffdb263000-7fffdda63000 rw-p 00000000 00:00 0
7fffdda63000-7fffdda64000 ---p 00000000 00:00 0
7fffdda64000-7fffe2264000 rw-p 00000000 00:00 0
7fffe267a000-7fffe273a000 rw-p 00000000 00:00 0
7fffe273a000-7fffe2863000 r-xp 00000000 08:01 903                        /usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0
7fffe2863000-7fffe2a62000 ---p 00129000 08:01 903                        /usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0
7fffe2a62000-7fffe2a63000 r--p 00128000 08:01 903                        /usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0
7fffe2a63000-7fffe2a65000 rw-p 00129000 08:01 903                        /usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0
7fffe2a65000-7fffe4a65000 rw-p 00000000 00:00 0
7fffe4a9b000-7fffe4b1b000 rw-p 00000000 00:00 0
7fffe4b1b000-7fffe4c53000 r-xp 00000000 08:06 42078660                   /home/chadi/src/hsl/20190503/lib/libcoinhsl.so.0.0.0
7fffe4c53000-7fffe4e52000 ---p 00138000 08:06 42078660                   /home/chadi/src/hsl/20190503/lib/libcoinhsl.so.0.0.0
7fffe4e52000-7fffe4e53000 r--p 00137000 08:06 42078660                   /home/chadi/src/hsl/20190503/lib/libcoinhsl.so.0.0.0
7fffe4e53000-7fffe4e54000 rw-p 00138000 08:06 42078660                   /home/chadi/src/hsl/20190503/lib/libcoinhsl.so.0.0.0
7fffe4e54000-7fffe4ec3000 r-xp 00000000 08:06 44437449                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmetis.so
7fffe4ec3000-7fffe50c2000 ---p 0006f000 08:06 44437449                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmetis.so
7fffe50c2000-7fffe50c4000 r--p 0006e000 08:06 44437449                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmetis.so
7fffe50c4000-7fffe50c5000 rw-p 00070000 08:06 44437449                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmetis.so
7fffe50c5000-7fffe50de000 r--p 00000000 08:06 43522285                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdmumps_seq-5.1.2.so
7fffe50de000-7fffe5230000 r-xp 00019000 08:06 43522285                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdmumps_seq-5.1.2.so
7fffe5230000-7fffe525d000 r--p 0016b000 08:06 43522285                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdmumps_seq-5.1.2.so
7fffe525d000-7fffe5260000 r--p 00197000 08:06 43522285                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdmumps_seq-5.1.2.so
7fffe5260000-7fffe5261000 rw-p 0019a000 08:06 43522285                   /home/chadi/miniconda3/envs/ocp-covid/lib/libdmumps_seq-5.1.2.so
7fffe5261000-7fffe7266000 rw-p 00000000 00:00 0
7fffe7290000-7fffe72d0000 rw-p 00000000 00:00 0
7fffe72d0000-7fffe72ec000 r--p 00000000 08:06 44043001                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotch-6.so
7fffe72ec000-7fffe7353000 r-xp 0001c000 08:06 44043001                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotch-6.so
7fffe7353000-7fffe736c000 r--p 00083000 08:06 44043001                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotch-6.so
7fffe736c000-7fffe736f000 r--p 0009b000 08:06 44043001                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotch-6.so
7fffe736f000-7fffe7373000 rw-p 0009e000 08:06 44043001                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotch-6.so
7fffe7373000-7fffe7374000 rw-p 00000000 00:00 0
7fffe7374000-7fffe7381000 r--p 00000000 08:06 43522283                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmumps_common_seq-5.1.2.so
7fffe7381000-7fffe73c9000 r-xp 0000d000 08:06 43522283                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmumps_common_seq-5.1.2.so
7fffe73c9000-7fffe73d3000 r--p 00055000 08:06 43522283                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmumps_common_seq-5.1.2.so
7fffe73d3000-7fffe73d4000 ---p 0005f000 08:06 43522283                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmumps_common_seq-5.1.2.so
7fffe73d4000-7fffe73d5000 r--p 0005f000 08:06 43522283                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmumps_common_seq-5.1.2.so
7fffe73d5000-7fffe73d6000 rw-p 00060000 08:06 43522283                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmumps_common_seq-5.1.2.so
7fffe73d6000-7fffe73d7000 rw-p 00000000 00:00 0
7fffe73d7000-7fffe7435000 r--p 00000000 08:06 44565883                   /home/chadi/miniconda3/envs/ocp-covid/lib/libipopt.so.1.10.12
7fffe7435000-7fffe75b7000 r-xp 0005e000 08:06 44565883                   /home/chadi/miniconda3/envs/ocp-covid/lib/libipopt.so.1.10.12
7fffe75b7000-7fffe760d000 r--p 001e0000 08:06 44565883                   /home/chadi/miniconda3/envs/ocp-covid/lib/libipopt.so.1.10.12
7fffe760d000-7fffe7616000 r--p 00235000 08:06 44565883                   /home/chadi/miniconda3/envs/ocp-covid/lib/libipopt.so.1.10.12
7fffe7616000-7fffe7617000 rw-p 0023e000 08:06 44565883                   /home/chadi/miniconda3/envs/ocp-covid/lib/libipopt.so.1.10.12
7fffe7617000-7fffe7622000 r--p 00000000 08:06 44305075                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi_nlpsol_ipopt.so.3.5
7fffe7622000-7fffe7636000 r-xp 0000b000 08:06 44305075                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi_nlpsol_ipopt.so.3.5
7fffe7636000-7fffe7664000 r--p 0001f000 08:06 44305075                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi_nlpsol_ipopt.so.3.5
7fffe7664000-7fffe7665000 ---p 0004d000 08:06 44305075                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi_nlpsol_ipopt.so.3.5
7fffe7665000-7fffe7666000 r--p 0004d000 08:06 44305075                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi_nlpsol_ipopt.so.3.5
7fffe7666000-7fffe7667000 rw-p 0004e000 08:06 44305075                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi_nlpsol_ipopt.so.3.5
7fffe7667000-7fffe9a67000 rw-p 00000000 00:00 0
7fffe9a69000-7fffeca69000 rw-p 00000000 00:00 0
7fffeca6a000-7fffef26a000 rw-p 00000000 00:00 0
7fffef281000-7fffef481000 rw-p 00000000 00:00 0
7fffef481000-7fffef539000 r--p 00000000 08:06 44182409                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi.so.3.5
7fffef539000-7fffef967000 r-xp 000b8000 08:06 44182409                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi.so.3.5
7fffef967000-7fffefa4f000 r--p 004e6000 08:06 44182409                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi.so.3.5
7fffefa4f000-7fffefa69000 r--p 005cd000 08:06 44182409                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi.so.3.5
7fffefa69000-7fffefa6a000 rw-p 005e7000 08:06 44182409                   /home/chadi/miniconda3/envs/ocp-covid/lib/libcasadi.so.3.5
7fffefa6a000-7ffff1a6b000 rw-p 00000000 00:00 0
7ffff1a74000-7ffff1a7e000 r--p 00000000 08:06 43388426                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgomp.so.1.0.0
7ffff1a7e000-7ffff1a95000 r-xp 0000a000 08:06 43388426                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgomp.so.1.0.0
7ffff1a95000-7ffff1a9e000 r--p 00021000 08:06 43388426                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgomp.so.1.0.0
7ffff1a9e000-7ffff1a9f000 ---p 0002a000 08:06 43388426                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgomp.so.1.0.0
7ffff1a9f000-7ffff1aa0000 r--p 0002a000 08:06 43388426                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgomp.so.1.0.0
7ffff1aa0000-7ffff1aa1000 rw-p 0002b000 08:06 43388426                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgomp.so.1.0.0
7ffff1aa1000-7ffff1f61000 rw-p 00000000 00:00 0
7ffff1f61000-7ffff1fb3000 r--p 00000000 08:06 44305105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/_casadi.so
7ffff1fb3000-7ffff2182000 r-xp 00052000 08:06 44305105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/_casadi.so
7ffff2182000-7ffff2258000 r--p 00221000 08:06 44305105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/_casadi.so
7ffff2258000-7ffff2259000 ---p 002f7000 08:06 44305105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/_casadi.so
7ffff2259000-7ffff225e000 r--p 002f7000 08:06 44305105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/_casadi.so
7ffff225e000-7ffff226c000 rw-p 002fc000 08:06 44305105                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/_casadi.so
7ffff226c000-7ffff426c000 rw-p 00000000 00:00 0
7ffff4276000-7ffff4277000 r--p 00000000 08:06 44042992                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotcherr-6.so
7ffff4277000-7ffff4278000 r-xp 00001000 08:06 44042992                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotcherr-6.so
7ffff4278000-7ffff4279000 r--p 00002000 08:06 44042992                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotcherr-6.so
7ffff4279000-7ffff427a000 r--p 00002000 08:06 44042992                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotcherr-6.so
7ffff427a000-7ffff427b000 rw-p 00003000 08:06 44042992                   /home/chadi/miniconda3/envs/ocp-covid/lib/libscotcherr-6.so
7ffff427b000-7ffff427e000 r--p 00000000 08:06 43522280                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmpiseq_seq-5.1.2.so
7ffff427e000-7ffff4282000 r-xp 00003000 08:06 43522280                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmpiseq_seq-5.1.2.so
7ffff4282000-7ffff4285000 r--p 00007000 08:06 43522280                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmpiseq_seq-5.1.2.so
7ffff4285000-7ffff4286000 r--p 00009000 08:06 43522280                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmpiseq_seq-5.1.2.so
7ffff4286000-7ffff4287000 rw-p 0000a000 08:06 43522280                   /home/chadi/miniconda3/envs/ocp-covid/lib/libmpiseq_seq-5.1.2.so
7ffff4287000-7ffff428b000 r--p 00000000 08:06 43522281                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpord_seq-5.1.2.so
7ffff428b000-7ffff429c000 r-xp 00004000 08:06 43522281                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpord_seq-5.1.2.so
7ffff429c000-7ffff42a0000 r--p 00015000 08:06 43522281                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpord_seq-5.1.2.so
7ffff42a0000-7ffff42a1000 r--p 00018000 08:06 43522281                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpord_seq-5.1.2.so
7ffff42a1000-7ffff42a2000 rw-p 00019000 08:06 43522281                   /home/chadi/miniconda3/envs/ocp-covid/lib/libpord_seq-5.1.2.so
7ffff42a2000-7ffff4422000 rw-p 00000000 00:00 0
7ffff4423000-7ffff4425000 r--p 00000000 08:06 44042998                   /home/chadi/miniconda3/envs/ocp-covid/lib/libesmumps-6.so
7ffff4425000-7ffff4427000 r-xp 00002000 08:06 44042998                   /home/chadi/miniconda3/envs/ocp-covid/lib/libesmumps-6.so
7ffff4427000-7ffff4428000 r--p 00004000 08:06 44042998                   /home/chadi/miniconda3/envs/ocp-covid/lib/libesmumps-6.so
7ffff4428000-7ffff4429000 r--p 00004000 08:06 44042998                   /home/chadi/miniconda3/envs/ocp-covid/lib/libesmumps-6.so
7ffff4429000-7ffff442a000 rw-p 00005000 08:06 44042998                   /home/chadi/miniconda3/envs/ocp-covid/lib/libesmumps-6.so
7ffff442a000-7ffff442b000 r--p 00000000 08:06 44567070                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/fcntl.cpython-37m-x86_64-linux-gnu.so
7ffff442b000-7ffff442d000 r-xp 00001000 08:06 44567070                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/fcntl.cpython-37m-x86_64-linux-gnu.so
7ffff442d000-7ffff442e000 r--p 00003000 08:06 44567070                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/fcntl.cpython-37m-x86_64-linux-gnu.so
7ffff442e000-7ffff442f000 r--p 00003000 08:06 44567070                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/fcntl.cpython-37m-x86_64-linux-gnu.so
7ffff442f000-7ffff4430000 rw-p 00004000 08:06 44567070                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/fcntl.cpython-37m-x86_64-linux-gnu.so
7ffff4430000-7ffff4433000 r--p 00000000 08:06 44567064                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/termios.cpython-37m-x86_64-linux-gnu.so
7ffff4433000-7ffff4434000 r-xp 00003000 08:06 44567064                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynl*** Error in `/home/chadi/miniconda3/envs/ocp-covid/bin/python': free(): invalid next size (normal): 0x00007fff7411d4c0 ***
oad/termios.cpython-37m-x86_64-linux-gnu.so
7ffff4434000-7ffff4435000 r--p 00004000 08:06 44567064                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/termios.cpython-37m-x86_64-linux-gnu.so
7ffff4435000-7ffff4436000 r--p 00004000 08:06 44567064                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/termios.cpython-37m-x86_64-linux-gnu.so
7ffff4436000-7ffff4438000 rw-p 00005000 08:06 44567064                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/termios.cpython-37m-x86_64-linux-gnu.so
7ffff4438000-7ffff443a000 r--p 00000000 08:06 44567065                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_multiprocessing.cpython-37m-x86_64-linux-gnu.so
7ffff443a000-7ffff443b000 r-xp 00002000 08:06 44567065                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_multiprocessing.cpython-37m-x86_64-linux-gnu.so
7ffff443b000-7ffff443c000 r--p 00003000 08:06 44567065                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_multiprocessing.cpython-37m-x86_64-linux-gnu.so
7ffff443c000-7ffff443d000 r--p 00003000 08:06 44567065                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_multiprocessing.cpython-37m-x86_64-linux-gnu.so
7ffff443d000-7ffff443e000 rw-p 00004000 08:06 44567065                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_multiprocessing.cpython-37m-x86_64-linux-gnu.so
7ffff443e000-7ffff453e000 rw-p 00000000 00:00 0
7ffff453e000-7ffff4548000 r--p 00000000 08:06 44567111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ssl.cpython-37m-x86_64-linux-gnu.so
7ffff4548000-7ffff4553000 r-xp 0000a000 08:06 44567111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ssl.cpython-37m-x86_64-linux-gnu.so
7ffff4553000-7ffff4559000 r--p 00015000 08:06 44567111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ssl.cpython-37m-x86_64-linux-gnu.so
7ffff4559000-7ffff455a000 ---p 0001b000 08:06 44567111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ssl.cpython-37m-x86_64-linux-gnu.so
7ffff455a000-7ffff455b000 r--p 0001b000 08:06 44567111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ssl.cpython-37m-x86_64-linux-gnu.so
7ffff455b000-7ffff4560000 rw-p 0001c000 08:06 44567111                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ssl.cpython-37m-x86_64-linux-gnu.so
7ffff4560000-7ffff45a0000 rw-p 00000000 00:00 0
7ffff45a0000-7ffff45a5000 r--p 00000000 08:06 44567112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_socket.cpython-37m-x86_64-linux-gnu.so
7ffff45a5000-7ffff45b3000 r-xp 00005000 08:06 44567112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_socket.cpython-37m-x86_64-linux-gnu.so
7ffff45b3000-7ffff45b8000 r--p 00013000 08:06 44567112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_socket.cpython-37m-x86_64-linux-gnu.so
7ffff45b8000-7ffff45b9000 r--p 00017000 08:06 44567112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_socket.cpython-37m-x86_64-linux-gnu.so
7ffff45b9000-7ffff45be000 rw-p 00018000 08:06 44567112                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_socket.cpython-37m-x86_64-linux-gnu.so
7ffff45be000-7ffff46be000 rw-p 00000000 00:00 0
7ffff46be000-7ffff46c4000 r--p 00000000 08:06 44307339                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/shapely/speedups/_speedups.cpython-37m-x86_64-linux-gnu.so
7ffff46c4000-7ffff46dd000 r-xp 00006000 08:06 44307339                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/shapely/speedups/_speedups.cpython-37m-x86_64-linux-gnu.so
7ffff46dd000-7ffff46e1000 r--p 0001f000 08:06 44307339                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/shapely/speedups/_speedups.cpython-37m-x86_64-linux-gnu.so
7ffff46e1000-7ffff46e2000 r--p 00022000 08:06 44307339                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/shapely/speedups/_speedups.cpython-37m-x86_64-linux-gnu.so
7ffff46e2000-7ffff46e4000 rw-p 00023000 08:06 44307339                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/shapely/speedups/_speedups.cpython-37m-x86_64-linux-gnu.so
7ffff46e4000-7ffff47a5000 rw-p 00000000 00:00 0
7ffff47a5000-7ffff47ab000 r--p 00000000 08:06 44304354                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_transformer.cpython-37m-x86_64-linux-gnu.so
7ffff47ab000-7ffff47c2000 r-xp 00006000 08:06 44304354                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_transformer.cpython-37m-x86_64-linux-gnu.so
7ffff47c2000-7ffff47c6000 r--p 0001d000 08:06 44304354                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_transformer.cpython-37m-x86_64-linux-gnu.so
7ffff47c6000-7ffff47c7000 r--p 00020000 08:06 44304354                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_transformer.cpython-37m-x86_64-linux-gnu.so
7ffff47c7000-7ffff47ca000 rw-p 00021000 08:06 44304354                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_transformer.cpython-37m-x86_64-linux-gnu.so
7ffff47ca000-7ffff480a000 rw-p 00000000 00:00 0
7ffff480a000-7ffff480f000 r--p 00000000 08:06 44304353                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_proj.cpython-37m-x86_64-linux-gnu.so
7ffff480f000-7ffff481f000 r-xp 00005000 08:06 44304353                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_proj.cpython-37m-x86_64-linux-gnu.so
7ffff481f000-7ffff4822000 r--p 00015000 08:06 44304353                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_proj.cpython-37m-x86_64-linux-gnu.so
7ffff4822000-7ffff4823000 r--p 00017000 08:06 44304353                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_proj.cpython-37m-x86_64-linux-gnu.so
7ffff4823000-7ffff4825000 rw-p 00018000 08:06 44304353                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_proj.cpython-37m-x86_64-linux-gnu.so
7ffff4825000-7ffff4866000 rw-p 00000000 00:00 0
7ffff4866000-7ffff486a000 r--p 00000000 08:06 44304351                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_geod.cpython-37m-x86_64-linux-gnu.so
7ffff486a000-7ffff4876000 r-xp 00004000 08:06 44304351                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_geod.cpython-37m-x86_64-linux-gnu.so
7ffff4876000-7ffff4878000 r--p 00010000 08:06 44304351                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_geod.cpython-37m-x86_64-linux-gnu.so
7ffff4878000-7ffff4879000 ---p 00012000 08:06 44304351                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_geod.cpython-37m-x86_64-linux-gnu.so
7ffff4879000-7ffff487a000 r--p 00012000 08:06 44304351                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_geod.cpython-37m-x86_64-linux-gnu.so
7ffff487a000-7ffff487c000 rw-p 00013000 08:06 44304351                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_geod.cpython-37m-x86_64-linux-gnu.so
7ffff487c000-7ffff4888000 r--p 00000000 08:06 44304355                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_crs.cpython-37m-x86_64-linux-gnu.so
7ffff4888000-7ffff48d2000 r-xp 0000c000 08:06 44304355                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_crs.cpython-37m-x86_64-linux-gnu.so
7ffff48d2000-7ffff48de000 r--p 00056000 08:06 44304355                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_crs.cpython-37m-x86_64-linux-gnu.so
7ffff48de000-7ffff48df000 ---p 00062000 08:06 44304355                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_crs.cpython-37m-x86_64-linux-gnu.so
7ffff48df000-7ffff48e0000 r--p 00062000 08:06 44304355                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_crs.cpython-37m-x86_64-linux-gnu.so
7ffff48e0000-7ffff48ec000 rw-p 00063000 08:06 44304355                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_crs.cpython-37m-x86_64-linux-gnu.so
7ffff48ec000-7ffff492e000 rw-p 00000000 00:00 0
7ffff492e000-7ffff4932000 r--p 00000000 08:06 44304352                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_list.cpython-37m-x86_64-linux-gnu.so
7ffff4932000-7ffff4943000 r-xp 00004000 08:06 44304352                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_list.cpython-37m-x86_64-linux-gnu.so
7ffff4943000-7ffff4945000 r--p 00015000 08:06 44304352                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_list.cpython-37m-x86_64-linux-gnu.so
7ffff4945000-7ffff4946000 r--p 00016000 08:06 44304352                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_list.cpython-37m-x86_64-linux-gnu.so
7ffff4946000-7ffff4948000 rw-p 00017000 08:06 44304352                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_list.cpython-37m-x86_64-linux-gnu.so
7ffff4948000-7ffff494b000 r--p 00000000 08:06 44304350                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_datadir.cpython-37m-x86_64-linux-gnu.so
7ffff494b000-7ffff4950000 r-xp 00003000 08:06 44304350                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_datadir.cpython-37m-x86_64-linux-gnu.so
7ffff4950000-7ffff4951000 r--p 00008000 08:06 44304350                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_datadir.cpython-37m-x86_64-linux-gnu.so
7ffff4951000-7ffff4952000 ---p 00009000 08:06 44304350                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_datadir.cpython-37m-x86_64-linux-gnu.so
7ffff4952000-7ffff4953000 r--p 00009000 08:06 44304350                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_datadir.cpython-37m-x86_64-linux-gnu.so
7ffff4953000-7ffff4954000 rw-p 0000a000 08:06 44304350                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/pyproj/_datadir.cpython-37m-x86_64-linux-gnu.so
7ffff4954000-7ffff4988000 r--p 00000000 08:06 44043732                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex.so.6.1.1
7ffff4988000-7ffff4a14000 r-xp 00034000 08:06 44043732                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex.so.6.1.1
7ffff4a14000-7ffff4a34000 r--p 000c0000 08:06 44043732                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex.so.6.1.1
7ffff4a34000-7ffff4a3a000 r--p 000df000 08:06 44043732                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex.so.6.1.1
7ffff4a3a000-7ffff4a3b000 rw-p 000e5000 08:06 44043732                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex.so.6.1.1
7ffff4a3b000-7ffff4a45000 r--p 00000000 08:06 44043731                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex_c.so.6.1.1
7ffff4a45000-7ffff4a62000 r-xp 0000a000 08:06 44043731                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex_c.so.6.1.1
7ffff4a62000-7ffff4a6a000 r--p 00027000 08:06 44043731                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex_c.so.6.1.1
7ffff4a6a000-7ffff4a6b000 ---p 0002f000 08:06 44043731                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex_c.so.6.1.1
7ffff4a6b000-7ffff4a6c000 r--p 0002f000 08:06 44043731                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex_c.so.6.1.1
7ffff4a6c000-7ffff4a6d000 rw-p 00030000 08:06 44043731                   /home/chadi/miniconda3/envs/ocp-covid/lib/libspatialindex_c.so.6.1.1
7ffff4a6d000-7ffff4a70000 r--p 00000000 08:06 43387743                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgcc_s.so.1
7ffff4a70000-7ffff4a7c000 r-xp 00003000 08:06 43387743                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgcc_s.so.1
7ffff4a7c000-7ffff4a7f000 r--p 0000f000 08:06 43387743                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgcc_s.so.1
7ffff4a7f000-7ffff4a80000 r--p 00011000 08:06 43387743                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgcc_s.so.1
7ffff4a80000-7ffff4a81000 rw-p 00012000 08:06 43387743                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgcc_s.so.1
7ffff4a81000-7ffff4a84000 r--p 00000000 08:06 43388104                   /home/chadi/miniconda3/envs/ocp-covid/lib/libquadmath.so.0.0.0
7ffff4a84000-7ffff4aa2000 r-xp 00003000 08:06 43388104                   /home/chadi/miniconda3/envs/ocp-covid/lib/libquadmath.so.0.0.0
7ffff4aa2000-7ffff4ab8000 r--p 00021000 08:06 43388104                   /home/chadi/miniconda3/envs/ocp-covid/lib/libquadmath.so.0.0.0
7ffff4ab8000-7ffff4ab9000 ---p 00037000 08:06 43388104                   /home/chadi/miniconda3/envs/ocp-covid/lib/libquadmath.so.0.0.0
7ffff4ab9000-7ffff4aba000 r--p 00037000 08:06 43388104                   /home/chadi/miniconda3/envs/ocp-covid/lib/libquadmath.so.0.0.0
7ffff4aba000-7ffff4abb000 rw-p 00038000 08:06 43388104                   /home/chadi/miniconda3/envs/ocp-covid/lib/libquadmath.so.0.0.0
7ffff4abb000-7ffff4ad6000 r--p 00000000 08:06 44303965                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgfortran.so.4.0.0
7ffff4ad6000-7ffff4bc5000 r-xp 0001b000 08:06 44303965                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgfortran.so.4.0.0
7ffff4bc5000-7ffff4be6000 r--p 0010a000 08:06 44303965                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgfortran.so.4.0.0
7ffff4be6000-7ffff4be7000 ---p 0012b000 08:06 44303965                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgfortran.so.4.0.0
7ffff4be7000-7ffff4be8000 r--p 0012b000 08:06 44303965                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgfortran.so.4.0.0
7ffff4be8000-7ffff4be9000 rw-p 0012c000 08:06 44303965                   /home/chadi/miniconda3/envs/ocp-covid/lib/libgfortran.so.4.0.0
7ffff4be9000-7ffff4bea000 rw-p 00000000 00:00 0
7ffff4bea000-7ffff4ccc000 r--p 00000000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff4ccc000-7ffff4ccd000 r-xp 000e2000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff4ccd000-7ffff4cd2000 ---p 000e3000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff4cd2000-7ffff65d9000 r-xp 000e8000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff65d9000-7ffff678b000 r--p 019ef000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff678b000-7ffff678c000 ---p 01ba1000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff678c000-7ffff6798000 r--p 01ba1000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff6798000-7ffff67a5000 rw-p 01bad000 08:06 44042192                   /home/chadi/miniconda3/envs/ocp-covid/lib/libopenblasp-r0.3.3.so
7ffff67a5000-7ffff67c3000 rw-p 00000000 00:00 0
7ffff67c3000-7ffff67eb000 r--p 00000000 08:06 44439123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so
7ffff67eb000-7ffff6a37000 r-xp 00028000 08:06 44439123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so
7ffff6a37000-7ffff6ab2000 r--p 00274000 08:06 44439123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so
7ffff6ab2000-7ffff6ab3000 ---p 002ef000 08:06 44439123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so
7ffff6ab3000-7ffff6ab6000 r--p 002ef000 08:06 44439123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so
7ffff6ab6000-7ffff6ad2000 rw-p 002f2000 08:06 44439123                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so
7ffff6ad2000-7ffff6bb2000 rw-p 00000000 00:00 0
7ffff6bb2000-7ffff6ed8000 r--p 00000000 08:01 1090                       /usr/lib/locale/locale-archive
7ffff6ed8000-7ffff6fe0000 r-xp 00000000 08:01 3148652                    /lib/x86_64-linux-gnu/libm-2.23.so
7ffff6fe0000-7ffff71df000 ---p 00108000 08:01 3148652                    /lib/x86_64-linux-gnu/libm-2.23.so
7ffff71df000-7ffff71e0000 r--p 00107000 08:01 3148652                    /lib/x86_64-linux-gnu/libm-2.23.so
7ffff71e0000-7ffff71e1000 rw-p 00108000 08:01 3148652                    /lib/x86_64-linux-gnu/libm-2.23.so
7ffff71e1000-7ffff71e8000 r-xp 00000000 08:01 3149220                    /lib/x86_64-linux-gnu/librt-2.23.so
7ffff71e8000-7ffff73e7000 ---p 00007000 08:01 3149220                    /lib/x86_64-linux-gnu/librt-2.23.so
7ffff73e7000-7ffff73e8000 r--p 00006000 08:01 3149220                    /lib/x86_64-linux-gnu/librt-2.23.so
7ffff73e8000-7ffff73e9000 rw-p 00007000 08:01 3149220                    /lib/x86_64-linux-gnu/librt-2.23.so
7ffff73e9000-7ffff73eb000 r-xp 00000000 08:01 3148647                    /lib/x86_64-linux-gnu/libutil-2.23.so
7ffff73eb000-7ffff75ea000 ---p 00002000 08:01 3148647                    /lib/x86_64-linux-gnu/libutil-2.23.so
7ffff75ea000-7ffff75eb000 r--p 00001000 08:01 3148647                    /lib/x86_64-linux-gnu/libutil-2.23.so
7ffff75eb000-7ffff75ec000 rw-p 00002000 08:01 3148647                    /lib/x86_64-linux-gnu/libutil-2.23.so
7ffff75ec000-7ffff75ef000 r-xp 00000000 08:01 3148646                    /lib/x86_64-linux-gnu/libdl-2.23.so
7ffff75ef000-7ffff77ee000 ---p 00003000 08:01 3148646                    /lib/x86_64-linux-gnu/libdl-2.23.so
7ffff77ee000-7ffff77ef000 r--p 00002000 08:01 3148646                    /lib/x86_64-linux-gnu/libdl-2.23.so
7ffff77ef000-7ffff77f0000 rw-p 00003000 08:01 3148646                    /lib/x86_64-linux-gnu/libdl-2.23.so
7ffff77f0000-7ffff79b0000 r-xp 00000000 08:01 3148648                    /lib/x86_64-linux-gnu/libc-2.23.so
7ffff79b0000-7ffff7bb0000 ---p 001c0000 08:01 3148648                    /lib/x86_64-linux-gnu/libc-2.23.so
7ffff7bb0000-7ffff7bb4000 r--p 001c0000 08:01 3148648                    /lib/x86_64-linux-gnu/libc-2.23.so
7ffff7bb4000-7ffff7bb6000 rw-p 001c4000 08:01 3148648                    /lib/x86_64-linux-gnu/libc-2.23.so
7ffff7bb6000-7ffff7bba000 rw-p 00000000 00:00 0
7ffff7bba000-7ffff7bd2000 r-xp 00000000 08:01 3148649                    /lib/x86_64-linux-gnu/libpthread-2.23.so
7ffff7bd2000-7ffff7dd1000 ---p 00018000 08:01 3148649                    /lib/x86_64-linux-gnu/libpthread-2.23.so
7ffff7dd1000-7ffff7dd2000 r--p 00017000 08:01 3148649                    /lib/x86_64-linux-gnu/libpthread-2.23.so
7ffff7dd2000-7ffff7dd3000 rw-p 00018000 08:01 3148649                    /lib/x86_64-linux-gnu/libpthread-2.23.so
7ffff7dd3000-7ffff7dd7000 rw-p 00000000 00:00 0
7ffff7dd7000-7ffff7dfd000 r-xp 00000000 08:01 3148659                    /lib/x86_64-linux-gnu/ld-2.23.so
7ffff7dfe000-7ffff7e00000 r--p 00000000 08:06 44439116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/lapack_lite.cpython-37m-x86_64-linux-gnu.so
7ffff7e00000-7ffff7e02000 r-xp 00002000 08:06 44439116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/lapack_lite.cpython-37m-x86_64-linux-gnu.so
7ffff7e02000-7ffff7e03000 r--p 00004000 08:06 44439116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/lapack_lite.cpython-37m-x86_64-linux-gnu.so
7ffff7e03000-7ffff7e04000 r--p 00004000 08:06 44439116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/lapack_lite.cpython-37m-x86_64-linux-gnu.so
7ffff7e04000-7ffff7e05000 rw-p 00005000 08:06 44439116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/linalg/lapack_lite.cpython-37m-x86_64-linux-gnu.so
7ffff7e05000-7ffff7e0d000 r--p 00000000 08:06 44567116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ctypes.cpython-37m-x86_64-linux-gnu.so
7ffff7e0d000-7ffff7e1f000 r-xp 00008000 08:06 44567116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ctypes.cpython-37m-x86_64-linux-gnu.so
7ffff7e1f000-7ffff7e26000 r--p 0001a000 08:06 44567116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ctypes.cpython-37m-x86_64-linux-gnu.so
7ffff7e26000-7ffff7e27000 r--p 00020000 08:06 44567116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ctypes.cpython-37m-x86_64-linux-gnu.so
7ffff7e27000-7ffff7e2b000 rw-p 00021000 08:06 44567116                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_ctypes.cpython-37m-x86_64-linux-gnu.so
7ffff7e2b000-7ffff7fab000 rw-p 00000000 00:00 0
7ffff7fab000-7ffff7fd2000 r--p 00000000 08:01 262388                     /usr/lib/locale/C.UTF-8/LC_CTYPE
7ffff7fd2000-7ffff7fd7000 rw-p 00000000 00:00 0
7ffff7fd9000-7ffff7fda000 rw-p 00000000 00:00 0
7ffff7fda000-7ffff7fdb000 rw-s 00000000 00:15 8                          /dev/shm/24S4pM (deleted)
7ffff7fdb000-7ffff7fdc000 rwxp 00000000 00:00 0
7ffff7fdc000-7ffff7fdf000 r--p 00000000 08:06 44567092                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/math.cpython-37m-x86_64-linux-gnu.so
7ffff7fdf000-7ffff7fe5000 r-xp 00003000 08:06 44567092                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/math.cpython-37m-x86_64-linux-gnu.so
7ffff7fe5000-7ffff7fe7000 r--p 00009000 08:06 44567092                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/math.cpython-37m-x86_64-linux-gnu.so
7ffff7fe7000-7ffff7fe8000 r--p 0000a000 08:06 44567092                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/math.cpython-37m-x86_64-linux-gnu.so
7ffff7fe8000-7ffff7fea000 rw-p 0000b000 08:06 44567092                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/math.cpython-37m-x86_64-linux-gnu.so
7ffff7fea000-7ffff7feb000 r--p 00000000 08:06 44567074                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_heapq.cpython-37m-x86_64-linux-gnu.so
7ffff7feb000-7ffff7fed000 r-xp 00001000 08:06 44567074                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_heapq.cpython-37m-x86_64-linux-gnu.so
7ffff7fed000-7ffff7fee000 r--p 00003000 08:06 44567074                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_heapq.cpython-37m-x86_64-linux-gnu.so
7ffff7fee000-7ffff7fef000 r--p 00003000 08:06 44567074                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_heapq.cpython-37m-x86_64-linux-gnu.so
7ffff7fef000-7ffff7ff1000 rw-p 00004000 08:06 44567074                   /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/lib-dynload/_heapq.cpython-37m-x86_64-linux-gnu.so
7ffff7ff1000-7ffff7ff8000 r--s 00000000 08:01 412119                     /usr/lib/x86_64-linux-gnu/gconv/gconv-modules.cache
7ffff7ff8000-7ffff7ffa000 r--p 00000000 00:00 0                          [vvar]
7ffff7ffa000-7ffff7ffc000 r-xp 00000000 00:00 0                          [vdso]
7ffff7ffc000-7ffff7ffd000 r--p 00025000 08:01 3148659                    /lib/x86_64-linux-gnu/ld-2.23.so
7ffff7ffd000-7ffff7ffe000 rw-p 00026000 08:01 3148659                    /lib/x86_64-linux-gnu/ld-2.23.so
7ffff7ffe000-7ffff7fff000 rw-p 00000000 00:00 0
7ffffffde000-7ffffffff000 rw-p 00000000 00:00 0                          [stack]
ffffffffff600000-ffffffffff601000 r-xp 00000000 00:00 0                  [vsyscall]

Thread 49 "python" received signal SIGSEGV, Segmentation fault.
[Switching to Thread 0x7fff91e74700 (LWP 8763)]
0x00007ffff59b9f9a in dcopy_k_HASWELL () from /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/../../../../libopenblas.so.0
(gdb)
```
and backtrace
```
#0  0x00007ffff59b9f9a in dcopy_k_HASWELL () from /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/../../../../libopenblas.so.0
#1  0x00007ffff59bf13d in dger_k_HASWELL () from /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/numpy/core/../../../../libopenblas.so.0
#2  0x00007fffacd6a31d in dgemv_ () from /home/chadi/miniconda3/envs/ocp-covid/lib/python3.7/site-packages/casadi/../../../././././liblapack.so.3
#3  0x00007fffe4bf4c8a in hsl_ma86_double::factor_solve_block (n=81, p=32, nb=32, s=0, a=..., la=2592, buf=..., d=...,
    perm=<error reading variable: value requires 143448 bytes, which is more than max-value-size>, q=2, tstats=..., control=..., flag=0) at hsl_ma86/hsl_ma86d.f90:4216
#4  0x00007fffe4bfdf52 in hsl_ma86_double::task_dispatch (scale=..., rhs_local=0x555580b996f0, total_threads=-21012, ldr=-19060, rhs=0x55557ca76750, nrhs=-18448, st=-1847116408,
    info=2075617428, tstats=..., control=..., nodes=..., blocks=..., stack=<error reading variable: value requires 524937920 bytes, which is more than max-value-size>,
    pos=<error reading variable: Cannot access memory at address 0x841868>, map=<error reading variable: Cannot access memory at address 0x841868>, lmap=0x55557cb948c0,
    lfact=0x55557cf320a0, maxn=<optimized out>, maxm=<optimized out>, nbcol=<optimized out>, invp=..., val=...) at hsl_ma86/hsl_ma86d.f90:2922
#5  __hsl_ma86_double_MOD_factorize_indef._omp_fn.1 () at hsl_ma86/hsl_ma86d.f90:2280
#6  0x00007ffff1a886d5 in gomp_thread_start (xdata=<optimized out>)
    at /home/nwani/m3/conda-bld/compilers_linux-64_1560109574129/work/.build/x86_64-conda_cos6-linux-gnu/src/gcc/libgomp/team.c:123
#7  0x00007ffff7bc16ba in start_thread (arg=0x7fff91e74700) at pthread_create.c:333
#8  0x00007ffff78f74dd in clone () at ../sysdeps/unix/sysv/linux/x86_64/clone.S:109
```