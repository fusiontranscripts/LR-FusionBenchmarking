#!/bin/bash

rm -rf ./_*
rm -f ./preds*
rm -f ./all.auc.*
rm -f ./okPara_ignoreUnsure*
rm -f ./ignoreUnsure*
rm -f auc_files.list  fusion_result_file_listing.dat  pipe.log

cd validated_fusions && ./cleanMe.sh

