#!/bin/sh

echo this is $0
for arg in "$@"
do
	echo argument is $arg
done
echo bye from $0
