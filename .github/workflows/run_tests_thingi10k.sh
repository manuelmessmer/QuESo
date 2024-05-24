#!/bin/bash

# Ids of models in thingi10k-database
ids=(36082 36069 36086 36372 36090 36373 37093 37179 37274 37276 37012 37278 37280 34785 37272 37282 32770 37275 37266 37283 37285 35269 37284 37287 37288 37322 37289 37323 37384 37506 37620 37622 37627 37743 37745 37750 37825 37841 37865 37866 37880 37881 37883 37886 37888 37928 37964 37967 37968 37972 37991 38290 38261 38291 38292 38293 38294 38295 38296 38297 38464 38497 38498 38636 38637 38638 38562 38639 38640 38641 38643 38644 38741 39010 39011 39012 39025 39026 39050 39158 39108 39165 39159 39166 39180 39182 39208 39245 39295 39349 39345 39344 39353 39355 39358 39495 39496 39498 39505 39499 39572)

count=0
for id in "${ids[@]}"
do
  found=1
  curl -o "TestExecutables/test.stl" https://www.thingiverse.com/download:$id || found=0

  if [ "$found" == "1" ]; then
    echo "File: $id"
    let "count+=1"
    ./TestExecutables/test_queso_thingi10k -- single TestExecutables/test.cstl 2 10 || exit 1
    ./TestExecutables/test_queso_thingi10k -- single TestExecutables/test.stl 3 30 || exit 1
  else
    echo File: $id is no longer available on Thingiverse.
  fi
  rm TestExecutables/test.stl
done
echo \#$count STLs are tested successfully.

