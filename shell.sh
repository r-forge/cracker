#!/bin/sh

echo $0 | sed 's|\(.*\)/.*|\1|'
ls $0 | sed 's|\(.*\)/.*|\1|'
Rscript ./cracker/start-app.R "$0"