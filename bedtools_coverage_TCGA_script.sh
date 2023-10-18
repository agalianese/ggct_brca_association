#!

#BEDTOOLS COVERAGE
file2_path="/Users/agalianese/jiayu/myIds"
token_path="/Users/agalianese/jiayu/ref/anqiToken"

# Read the contents of file2
myIds=$(<"$file2_path")

# Read the contents of file3
token=$(<"$token_path")

for myId in $myIds; do
        curl -H "X-Auth-Token: $token" "https://api.gdc.cancer.gov/slicing/view/$myId?region=chr7:30496621-30504829" -o slice.bam

        #run bedtools
        bedtools coverage -s -a /Users/agalianese/jiayu/ref/ggct_full_exon_hg38.bed -b slice.bam \
        -bed > $myId.ggct.bed

done
