#!/usr/bin/bash

java -Xmx8g -jar snpEff.jar -ud 20000 hg19 kunkle.vcf > kunkle.ann.vcf
grep -v intergenic kunkle.ann.vcf > kunkle.final.ann.vcf
