#!/usr/bin/bash

for i in *pdf; do
	sips -s format png ${i} ${i%.pdf}.png
done
