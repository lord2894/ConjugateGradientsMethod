mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 1 -w 00:15:00 --stdout ./resultsMP/rezultMP-1-1000.out ./sequential -- -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsMP/rezult_Par_1_1000.txt
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 128 -w 00:15:00 --stdout ./resultsMP/rezultMP-128-1000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsMP/rezult_Par_128_1000.txt
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 256 -w 00:10:00 --stdout ./resultsMP/rezultMP-256-1000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsMP/rezult_Par_256_1000.txt
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 512 -w 00:05:00 --stdout ./resultsMP/rezultMP-512-1000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsMP/rezult_Par_512_1000.txt
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 1 -w 00:15:00 --stdout ./resultsMP/rezultMP-1-2000.out ./sequential -- -x 0 -y 0 -a 2 -b 2 -r 2000 -c 2000 -f ./resultsMP/rezult_Par_1_2000.txt
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 128 -w 00:15:00 --stdout ./resultsMP/rezultMP-128-2000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 2000 -c 2000 -f ./resultsMP/rezult_Par_128_2000.txt
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 256 -w 00:10:00 --stdout ./resultsMP/rezultMP-256-2000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 2000 -c 2000 -f ./resultsMP/rezult_Par_256_2000.txt
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 512 -w 00:05:00 --stdout ./resultsMP/rezultMP-512-2000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 2000 -c 2000 -f ./resultsMP/rezult_Par_512_2000.txt

mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 2 -w 02:00:00 --stdout ./resultsMP/rezultMP-2-2000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 2000 -c 2000 -f ./resultsMP/rezult_Par_2_2000.txt
