#!/bin/bash
if [[ -n "$1" ]]; then
	dt=$1
else
	dt=$(date '+%d/%m/%Y %H:%M:%S')
fi

git pull
git add -A :/
git commit -m "$dt" 
git push -u origin master
echo "Committed/pushed changes to master branch on GitHub"
