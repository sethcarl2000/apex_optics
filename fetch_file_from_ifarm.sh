#!/bin/bash

# Usage: 
#   
#   $ ./fetch_file_from_ifarm.sh    full/path/to/file/on/ifarm     relative/path/to/destination/dir

username="sethhall"

scp -J "${username}@scilogin.jlab.org" "${username}@ifarm:${1}" "./${2}"