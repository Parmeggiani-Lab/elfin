#!/bin/bash

echo 'trap "exit" INT'
ls ./res/preprocessed/*/*.pdb | xargs -I{} echo './scripts/Shell/relax.sh {}'