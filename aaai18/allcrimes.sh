#!/usr/bin/sh

./runcrime --config=cfg-crime-nohidden.txt --out=crimedata/nohidden1.txt
./runcrime --config=cfg-crime-hidden.txt --out=crimedata/hidden1.txt

./runcrime --config=cfg-crime-nohidden.txt --out=crimedata/nohidden2.txt
./runcrime --config=cfg-crime-hidden.txt --out=crimedata/hidden2.txt

./runcrime --config=cfg-crime-nohidden.txt --out=crimedata/nohidden3.txt
./runcrime --config=cfg-crime-hidden.txt --out=crimedata/hidden3.txt

./runcrime --config=cfg-crime-nohidden.txt --out=crimedata/nohidden4.txt
./runcrime --config=cfg-crime-hidden.txt --out=crimedata/hidden4.txt

./runcrime --config=cfg-crime-nohidden.txt --out=crimedata/nohidden5.txt
./runcrime --config=cfg-crime-hidden.txt --out=crimedata/hidden5.txt

./runcrime --config=cfg-crime-hidden5.txt --out=crimedata/hidden6.txt
./runcrime --config=cfg-crime-hidden5.txt --out=crimedata/hidden7.txt
./runcrime --config=cfg-crime-hidden5.txt --out=crimedata/hidden8.txt
./runcrime --config=cfg-crime-hidden5.txt --out=crimedata/hidden9.txt
./runcrime --config=cfg-crime-hidden5.txt --out=crimedata/hidden10.txt

./runcrime --config=cfg-crime-nohidden.txt --lambda=0.01 --out=crimedata/h0l001.txt
./runcrime --config=cfg-crime-hidden5.txt --lambda=0.01 --out=crimedata/h5l001.txt

./runcrime --config=cfg-crime-nohidden.txt --lambda=1 --out=crimedata/h0-l100-1.txt

./runcrime --config=cfg-crime-nohidden.txt --lambda=10000000 --out=crimedata/h0-l7-0
./runcrime --config=cfg-crime-hidden5.txt --lambda=10000000 --out=crimedata/h5-l7-0
./runcrime --config=cfg-crime-hidden5.txt --lambda=10000000 --out=crimedata/h5-l7-1

./runcrime --config=cfg-crime-nohidden.txt --lambda=1000000 --out=crimedata/h0-l6-0
./runcrime --config=cfg-crime-hidden5.txt --lambda=1000000 --out=crimedata/h5-l6-0
./runcrime --config=cfg-crime-hidden5.txt --lambda=1000000 --out=crimedata/h5-l6-1

./runcrime --config=cfg-crime-nohidden.txt --lambda=100000 --out=crimedata/h0-l5-0
./runcrime --config=cfg-crime-hidden5.txt --lambda=100000 --out=crimedata/h5-l5-0
./runcrime --config=cfg-crime-hidden5.txt --lambda=100000 --out=crimedata/h5-l5-1 --nem=200

./runcrime --config=cfg-crime-nohidden.txt --lambda=10000 --out=crimedata/h0-l4-0 --nem=200
./runcrime --config=cfg-crime-hidden5.txt --lambda=10000 --out=crimedata/h5-l4-0 --nem=200
./runcrime --config=cfg-crime-hidden5.txt --lambda=10000 --out=crimedata/h5-l4-1 --nem=200

./runcrime --config=cfg-crime-hidden5.txt --lambda=100000 --clampWitt=50 --nem=100 --out=crimedata/h5-l5-c50-0
