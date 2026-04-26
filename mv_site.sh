#!/bin/bash
cp ./report/reports_main.html ./docs/index.html
rm -rf ./docs/reports_main_files
cp -R ./report/reports_main_files ./docs
