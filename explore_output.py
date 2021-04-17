import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import re
from itertools import groupby
import argparse

from .config_genhelper import convert_values

# rename ldhat data file
def rename_ldhat_file(fdir,custume_replace = []):
    '''
    Input:
        fdir: datafolder
        custume_replace: use to define other replace we want. Ex change "hello_work" to "hw", custume_replace =[['hello_work','hw']]
    '''
    fnames = os.listdir(fdir)
    for fname in fnames:
        new_name = fname
        for creplace in custume_replace:
            new_name = new_name.replace(creplace[0],creplace[1])
        new_name = new_name.replace('_acceptance_rates.txt','.acrates')
        new_name = new_name.replace('_hotspots.txt','.hotspots')
        new_name = new_name.replace('_summary.txt','.summary')
        new_name = new_name.replace('_rates.txt','.rates')
        old_path = os.path.join(fdir,fname)
        new_path = os.path.join(fdir,new_name)
        # print(old_path)
        # print(new_path)
        os.rename(r'{}'.format(old_path),r'{}'.format(new_path))

# item ccủa batch phải được đặt tên theo định dạng:
# <..phần tên trước định dạng..>_<tên batch(số hệu)>_<tên item của batch>.<loại file>
def get_batch_item_name(fdir,ftype,get_locs=False):
    assert ftype in ['summary'], 'ftype must in [summary]'
    ftype = '.'+ftype
    fnames = [os.path.join(fdir,name) for name in os.listdir(fdir) if os.path.splitext(name)[-1] == ftype]    
    gfnames = []
    for key,group in groupby(fnames, lambda x: x.split('_')[:-1]):
        gfnames.append(list(group))
    gflnames = []
    if get_locs == False:
        return gfnames, gflnames
    flnames = [os.path.join(fdir,name) for name in os.listdir(fdir) if os.path.splitext(name)[-1] == '.locs']
    for key,group in groupby(flnames, lambda x: x.split('_')[:-1]):
        gflnames.append(list(group))
    return gfnames, gflnames

def stack_data(data,poss):
    nb = len(data)
    redata = data[0].copy()
    reposs = poss[0].copy()
    for i in range(1,nb):
        lindexs = np.argwhere(reposs >= poss[i][0])
        rindexs = np.argwhere(poss[i] <= reposs[-1])
        for key in redata.keys():
            redata[key][lindexs] = (redata[key][lindexs]+data[i][key][rindexs])/2
            get_index = rindexs[-1][0]+1
            redata[key] = np.hstack([redata[key],data[i][key][get_index:]])
        reposs = np.hstack([reposs,poss[i][get_index:]])
    return redata, reposs

def get_data_from_batchs(gfnames,gflnames,col_data,sep='\t'):
    '''
    Input:
        gfnames: Group file name, danh sách file .summary được group lại theo nhóm
        gflnames: Group file name, danh sách file .locs được group lại theo nhóm
        col_data: danh sách các object, mỗi object có trường
                'col_name' : tên của column cần lấy dữ liệu
                'convert' : là func, đầu vô là data được lấy dựa trên col_name, đầu ra trả về giá trị tương ứng
    Ouput:
        data được lấy theo col_data và position tương ứng
    '''
    assert len(gfnames) == len(gflnames), 'gflnames and gfnames must same number element'
    assert [len(g) for g in gfnames] == [len(g) for g in gflnames], 'element each gfnames and gflnames must same length'
    nb = len(gfnames)
    gsnp_data = []
    gsnp_poss = []
    for i in range(nb):
        nbp = len(gfnames[i])
        snp_data = []
        snp_poss = []
        for j in range(nbp):
            df = pd.read_csv(gfnames[i][j],sep=sep)
            fposs = open(gflnames[i][j])
            poss = np.array(fposs.readline().split(),dtype=int)
            snp_poss.append(poss[:-1])
            fposs.close()
            data = {}
            for col in col_data:
                data[col['col_name']] = df[col['col_name']].values
                if 'convert' in col and col['convert'] is not None:
                    data[col['col_name']] = col['convert'](data[col['col_name']])
            snp_data.append(data)
        tdata, tposs = stack_data(snp_data,snp_poss)
        gsnp_data.append(tdata)
        gsnp_poss.append(tposs)
    return gsnp_data, gsnp_poss

def view_small_region(data,poss,from_location = 0, to_location=-1):
    '''
    Ouput:
        Trả về dữ liệu vùng region muốn xem
    '''
    if to_location == -1:
        to_location = poss[-1]
    indexs = np.argwhere((poss >= from_location) & (poss <= to_location))
    vposs = np.hstack([poss[indexs],poss[indexs]]).flatten()
    vdata = np.hstack([data[indexs],data[indexs]]).flatten()
    return vdata[1:], vposs[:-1]

def post_batch_summary_avg(fdir,col_data = [{'col_name': '4Ner/kb','rename':'Recomb_rate(cM/Mb)','convert': lambda x : x*convert_values._4Ner_per_kb_to_cM_per_Mb}],need_save = True):
    '''
    Input:
        fdir: thư mục cần lấy dữ liệu
        col_data: danh sách các object, mỗi object có trường
                'col_name' : tên của column cần lấy dữ liệu
                'convert' : là func, đầu vô là data được lấy dựa trên col_name, đầu ra trả về giá trị tương ứng
                'rename' : tên mới của column
        need_save: Có cần lưu lại hay không
    Output:
        dataframe được lấy dữ liệu theo col_name và position
    '''
    gfnames, gflnames = get_batch_item_name(fdir,ftype='summary',get_locs=True)
    group_data, group_poss = get_data_from_batchs(gfnames,gflnames,col_data=col_data)
    avg_data = None
    for e in group_data:
        if avg_data is None:
            avg_data = pd.DataFrame(e)
        else:
            avg_data = avg_data + pd.DataFrame(e)
    avg_data = avg_data/len(group_data)
    for col in col_data:
        if 'rename' in col and col['rename'] is not None:
            avg_data.rename(columns={col['col_name']: col['rename']},inplace=True)
    poss = group_poss[0]    
    if need_save:
        data = avg_data.copy()
        data.insert(0,'position(Mb)',poss/1e6)
        csv_name = os.path.basename(os.path.dirname(fdir))
        data.to_csv(csv_name+'_summary.csv')
    return avg_data, poss

def get_data_compare(fpath,rates_name,poss_name,sep='\t',query=None):
    compare_data = pd.read_csv(fpath,sep=sep)
    if query is not None:
       compare_data = compare_data.query(query)
    compare_rates = compare_data[rates_name].values
    compare_poss = compare_data[poss_name].values
    return compare_rates, compare_poss

def get_fail_data_liftovered(fpath):
    err_new_compare_file = open(fpath)
    flag = False
    err_poss = []
    for line in err_new_compare_file.readlines():
        if flag:
            err_poss.append(int(line.split('\t')[1]))
        if line.find('Deleted') != -1:
            flag = True
        else:
            flag = False
    err_new_compare_file.close()
    return err_poss

def get_index_liftovered(convered_path,fail_path,old_poss):
    new_compare_data = pd.read_csv(convered_path,sep='\t',header=None)
    new_compare_poss = new_compare_data.values[:,1]
    remove_indexs = []
    pmax = -1
    for i in range(len(new_compare_poss)):
        if pmax > new_compare_poss[i]:
            print(new_compare_poss[i])
            remove_indexs.append(i)
        else:
            pmax = new_compare_poss[i]
    remove_indexs = np.array(remove_indexs)[::-1]
    err_poss = get_fail_data_liftovered(fail_path)
    hold_indexs = np.array([i for i in range(len(old_poss)) if old_poss[i] not in err_poss])
    print(len(hold_indexs))
    for ri in remove_indexs:
        hold_indexs = np.delete(hold_indexs,ri)
        new_compare_poss = np.delete(new_compare_poss,ri)
    print(len(hold_indexs))
    return hold_indexs, new_compare_poss

def view_plot_data(view_data,name,figsize=(14,4),from_location = 0, to_location=-1,view_same=False):
    nbrow = len(view_data)
    nbcol = 1
    plt.figure(figsize=(figsize[0],figsize[1]*nbrow))
    axs = []
    max_value = 0
    max_index = 0
    for i in range(nbrow):
        if view_same is not True:
            ax = plt.subplot(nbrow,nbcol,i+1)
            axs.append(ax)
        rates = view_data[i]['rates']
        poss = view_data[i]['poss']
        rates, poss = view_small_region(rates,poss,from_location=from_location,to_location=to_location)
        plt.plot(poss,rates,label=view_data[i]['label'])
        if len(rates) == 0:
            continue
        max_r = max(rates)
        if max_value < max_r:
            max_value = max_r
            max_index = i        
        if view_same is not True:
            plt.ylabel(view_data[i]['ylabel'])
            plt.xlabel(view_data[i]['xlabel'])
            plt.title(view_data[i]['title'])

    if view_same is not True:
        for i in range(nbrow):
            if i != max_index:
                axs[i].sharex(axs[max_index])
                axs[i].sharey(axs[max_index])
    else:
        plt.ylabel(view_data[0]['ylabel'])
        plt.xlabel(view_data[0]['xlabel'])
        plt.title(view_data[0]['title'])
    plt.legend()
    plt.savefig(name,facecolor='white')
    plt.tight_layout()
    plt.show()

def summary(data_dir,compare_path,plot_name='No Name'):
    avg_data, poss = post_batch_summary_avg(data_dir)
    view_data = []
    view_data.append({
        'title': 'VN Recombination Rate chr22',
        'ylabel': 'Recombination Rate (cM/Mb)',
        'xlabel': 'Position (Mb)',
        'rates': avg_data.values.T[0],
        'poss':poss/1e6
    })
    compares = pd.read_csv(compare_path)
    for i in range(len(compares)):
        query = compares['query'][i]
        if pd.isnull(query):
            query = None
        compare_rates, compare_poss = get_data_compare(compares['path'][i],compares['rates'][i],compares['poss'][i],query=query)
        xlabel = compares['xlabel'][i] if pd.isnull(compares['xlabel'][i]) is not True else compares['poss'][i]
        ylabel = compares['ylabel'][i] if pd.isnull(compares['ylabel'][i]) is not True else compares['rates'][i]
        title = compares['name'][i]
        view_data.append({
            'title': title,
            'ylabel': ylabel,
            'xlabel': xlabel,
            'rates': compare_rates,
            'poss':compare_poss/1e6
        })
    view_plot_data(view_data,plot_name)

def get_data_summary(data_dirs,col_data):
    data_avg = None
    poss = None
    for data_dir in data_dirs:
        data, tposs = post_batch_summary_avg(data_dir,col_data=col_data,need_save=False)
        if data_avg is None:
            data_avg = data
        else:
            data_avg = data_avg + data
        poss = tposs
    data_avg = data_avg/len(data_dirs)
    data_avg.insert(0,'Position(bp)',poss)
    return data_avg

def map_lossed_data(data,cdata,lossed_name,value_name):
    '''
    Input:
    '''
    ldata = data.copy()
    lcdata = cdata.copy()
    temp = ldata.append(lcdata,ignore_index=True).sort_values(lossed_name,ignore_index=True)
    temp.drop_duplicates(subset=[lossed_name],inplace=True, ignore_index=True)
    if np.isnan(temp[value_name][0]):
        vmin = np.min(temp[value_name])
        if np.isnan(vmin):
            temp[value_name][0] = 0
        else:
            temp[value_name][0] = vmin
    temp.fillna(method='pad',inplace=True)
    return temp

def convert_4NerKb_to_cMMb(x):
    return x*convert_values._4Ner_per_kb_to_cM_per_Mb

def convert_4Ner_to_cM(x):
    return x*convert_values._4Ner_to_cM

def get_convert_func(name):
    if name == '4Ner/kb_to_cM/Mb':
        return convert_4NerKb_to_cMMb
    elif name == '4Ner_to_cM':
        return convert_4Ner_to_cM
    else:
        return lambda x: x

def get_col_data(fpath):
    col_name = 'col_name'
    convert = 'convert'
    rename = 'rename'
    csv = pd.read_csv(fpath)
    col_data = []
    for i, row in csv.iterrows():
        temp = {col_name: row[col_name]}
        if row[convert] is not np.NaN and row[convert] is not None:
            temp[convert] = get_convert_func(row[convert])
        if row[rename] is not np.NaN and row[rename] is not None:
            temp[rename] = row[rename]
        col_data.append(temp)
    return col_data

def get_indexs_position(poss,f_poss,t_poss):
    return np.argwhere((poss>=f_poss)&(poss <= t_poss)).flatten()

def draw_hist(data,dlabels,n,step,xshow=True,title=None,xlabel=None,ylabel=None,need_save=False,save_dir='./',f_value=None,t_value=None):
    '''
    input:
        data: dữ liệu cần vx hist
        n: làm tròn sau n chữ số sau dấu "."
        step: độ rộng mỗi bin
        xshow: hiện giá trị trục x
        xlabel: hiện label cho trục x
        ylabel: hiện label cho trục y
        need_save: lưu lại hình hay không
        save_dir: thư mục lưu
    '''
    center = 0.5*(10**(-n)) # Biến hỗ trợ làm tròn lên hay xuống cho hàm round    
    plt.figure(figsize=(20,7))
    lbound = -1
    rbound = -1
    for d in data:
        left_bound = np.round(d.min() - center,n)
        right_bound  = np.round(d.max() + center,n)
        right_bound = left_bound + ((right_bound-left_bound)//step + 1)*step+step/2
        if lbound > left_bound or lbound == -1:
            lbound = left_bound
        if rbound < right_bound or right_bound == -1:
            rbound = right_bound
    if f_value is not None:
        lbound = f_value
    if t_value is not None:
        rbound = t_value
    for i in range(len(data)):
        vindexs = get_indexs_position(data[i],lbound,rbound)
        plt.hist(data[i][vindexs],bins=np.arange(lbound,rbound,step=step),histtype='step',label=dlabels[i])
    plt.legend()
    if xshow:
        plt.xticks(ticks=np.arange(lbound,rbound,step=step))
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if need_save:
        fig_name = title
        fig_path = os.path.join(save_dir,fig_name)
        plt.savefig(fig_path,facecolor='white')
    plt.show()

def view_same_plot(data,rate_name,poss_name,title,from_location,to_location,lline=None,rline=None,need_save=False,save_dir=None):
    plt.figure(figsize=(14,7))
    ymax = -1
    for d in data:
        rates, poss = view_small_region(d[rate_name].values,d[poss_name].values,from_location=from_location,to_location=to_location)
        if rates.max() > ymax:
            ymax = rates.max()
        plt.plot(poss,rates)
    if lline is not None:
        plt.vlines(lline,0,ymax,linestyles='dashed')
    if rline is not None:
        plt.vlines(rline,0,ymax,linestyles='dashed')
    plt.xlabel(poss_name)
    plt.ylabel(rate_name)
    plt.title(title)
    if need_save and save_dir is not None:
        plt.savefig(os.path.join(save_dir,title),facecolor='white')
    plt.show()

def get_hotspot_position_to_bed(path,chr,chr_col,start_col,end_col,sep,skiprows,save_path=None):
    temp = pd.read_csv(path,sep=sep,skiprows=skiprows)
    temp = temp.query("{0} == '{1}'".format(chr_col,chr))
    temp = temp[[chr_col,start_col,end_col]]
    if save_path is not None:
        temp.to_csv(save_path,sep='\t',index=False,header=['#chrom','chromStart','chromEnd'])
    return temp

def get_data_hotspot(data,rate_name,poss_name,from_location,to_location,around_distance,scale=None):
    center = (from_location+to_location)/2
    indexs = get_indexs_position(data[poss_name].values,center-around_distance,center+around_distance)
    rates = data[rate_name].values[indexs]
    poss = data[poss_name].values[indexs]
    if len(poss) == 0:
        poss = np.array([center])
    poss = np.hstack([[center-around_distance],poss,[center+around_distance]])
    m = data[rate_name][data[poss_name] < poss[0]].values
    if (len(rates) == 0) and (len(m) > 0):
        rates = np.array([m[-1]])
    elif len(rates) == 0:
        rates = np.array([0])
    if len(m) > 0:
        rates = np.hstack([[m[-1]],rates,[rates[-1]]])
    else:
        rates = np.hstack([[0],rates,[rates[-1]]])
    poss = poss - center
    if scale is not None:
        poss = np.round(poss*scale,6)
    temp = pd.DataFrame({'0':poss,'1':rates})
    temp.columns = [poss_name,rate_name]
    return temp

def get_slide_hotspot_data(slidesdf,start_name,end_name,vrange,data,rate_name,poss_name,scale,view_range=None):
    nb_slides = len(slidesdf)
    if view_range is None:
        view_range = range(nb_slides)
    elif view_range[-1] > nb_slides:
        view_range = range(view_range[0],nb_slides)
    data_hotspot = pd.DataFrame({poss_name:[],rate_name:[]})
    for i in view_range:
        row = slidesdf.iloc[i]
        temp = get_data_hotspot(data,rate_name,poss_name,row[start_name],row[end_name],vrange,scale=scale)
        # Map position của 2 dữ liệu vô nhau
        # Lấy position của hotspot
        clone = data_hotspot.copy()
        clone[rate_name] = np.NaN
        # Gọi hàm map vô dữ liệu hiện tại và fill dữ liệu
        temp = map_lossed_data(temp,clone,poss_name,rate_name)
        # Lấy position của dữ liệu hiện tại
        clone = temp.copy()
        clone[rate_name] = np.NaN
        # Map vô dữ liệu hotspot
        data_hotspot = map_lossed_data(data_hotspot,clone,poss_name,rate_name)
        # Cộng dồn giá trị rate
        data_hotspot[rate_name] = data_hotspot[rate_name] + temp[rate_name]
    return data_hotspot

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-data', type=str)
    parser.add_argument('-compare', type=str)
    parser.add_argument('-plot',default='No Name', type=str)
    args = parser.parse_args()    
    summary(args.data,args.compare,args.plot)