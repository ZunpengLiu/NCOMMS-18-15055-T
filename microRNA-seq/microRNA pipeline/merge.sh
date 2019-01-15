ZunpengdeMacBook-Pro:mirdeep2 zunpengliu$ awk  'FS="," {arr[$1] = arr[$1] + $2}END{for (a in arr) print a "," arr[a]}' DN1_mature_count.csv >DN1_mature_count_merge.csv
ZunpengdeMacBook-Pro:mirdeep2 zunpengliu$ awk  'FS="," {arr[$1] = arr[$1] + $2}END{for (a in arr) print a "," arr[a]}' DN2_miRNA_count.csv >DN2_miRNA_count_merge.csv
ZunpengdeMacBook-Pro:mirdeep2 zunpengliu$ awk  'FS="," {arr[$1] = arr[$1] + $2}END{for (a in arr) print a "," arr[a]}' DT1_miRNA_count.csv > DT1_miRNA_count_merge.csv
ZunpengdeMacBook-Pro:mirdeep2 zunpengliu$ awk  'FS="," {arr[$1] = arr[$1] + $2}END{for (a in arr) print a "," arr[a]}' DT2_miRNA_count.csv >DT2_miRNA_count_merge.csv

