ls *.coverage | parallel --keep-order -j 60 --bar 'cat {} | awk '"'"'{if($3>=5)print}'"'"' | wc -l' > coverage_5
ls *.coverage | parallel --keep-order -j 60 --bar 'cat {} | awk '"'"'{if($3>=10)print}'"'"' | wc -l' > coverage_10
ls *.coverage | parallel --keep-order -j 60 --bar 'cat {} | awk '"'"'{if($3>=20)print}'"'"' | wc -l' > coverage_20