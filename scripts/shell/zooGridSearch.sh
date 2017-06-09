#!/bin/bash

export gsVer=2
export queueName=batch
export ppn=256:knl
export startGridId=0
export endGridId=80
./doGridSearch
