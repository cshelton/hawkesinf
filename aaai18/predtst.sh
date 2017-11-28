./runcrime --config=cfg-test-hidden0.txt --out=h0lastyr.model
./runcrime --config=cfg-test-hidden5.txt --out=h5lastyr.model

./tstcrime --in=h0lastyr.model --out=h0lastyrpred.txt
./tstcrime --in=h5lastyr.model --out=h5lastyrpred.txt
