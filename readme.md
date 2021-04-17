 # summary_col_data file content format
 first is header.
 col_name: column name will be used to get data
 conver: func name, will used to get func convert value data
 rename: after covert, if need rename that column, this field will have data
 file content format as csv file, sep = ","

 btw explore_ouput are lib, it main func dont have much summary. To view what that lib can do, please check explore_output.ipynb
 # where we can find data was used
 hotspot 22A: https://www.ncbi.nlm.nih.gov/gene/107832855
 click downdataset -> gen sequence (fasta) -> download
 hotspot statistic: https://raw.githubusercontent.com/auton1/Campbell_et_al/master/hotspots-coldspots-b37-filtered.txt
 # how to run explore_output main func

 python explore_output.py -data <đường dẫn đến thư mục chứa kết quả của ldhat/rhomap> -compare <file định dạng nguồn dữ liệu cần compare>

 - compare: file đưa vô đây là file csv gồm các trường:
    'path':đường dẫn tới file compare,
    'rates': tên của cột để lấy giá trị recombination rate,
    'poss': tên cột để lấy position,
    'name': tên của biểu đồ được vẽ lên,
    'ylabel': tên của trục y (mặc định lấy theo rates),
    'xlabel': tên của trục x (mặc định lấy theo poss),
    'query': câu query để filter dữ liệu pandas DataFrame

>template_compare = pd.DataFrame({
    'path':['../HapmapII_GRCh37_RecombinationHotspots/genetic_map_GRCh37_chr22.txt','./genetic_map_GRCh38_merged.tab'],
    'rates': ['Rate(cM/Mb)','recomb_rate'],
    'poss': ['Position(bp)','pos'],
    'name': ['genetic_map_GRCh38_chr22','genetic_map_GRCh38_chr22_merged'],
    'ylabel':['Recombination Rate (cM/Mb)','Recombination Rate (cM/Mb)'],
    'xlabel': ['Position (Mb)','Position (Mb)'],
    'query':['',"chrom == 'chr22'"]
})
template_compare.to_csv('template_compare.csv')

Tạo requirements file
> pipreqs ./ --force