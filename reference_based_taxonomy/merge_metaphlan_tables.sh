#!/bin/bash

output_file=$1 #e.g.AZ107

python /mnt/pkgs/metaphlan2/utils/merge_metaphlan_tables.py *metaphlan > "${output_file}"
