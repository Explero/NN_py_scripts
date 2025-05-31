import os
import subprocess

def extract_tar_gz(folder_path):
    #folder_path = './test'
    file_names = os.listdir(folder_path)
    tar_files = [f for f in file_names if f.endswith('.tar.gz')]

    count_file=0
    file_no_suffixs=[]
    file_work_dirs=[]   # 解压所有 .tar.gz 文件并重命名
    for file in tar_files:
        count_file+=1
        file_no_suffix = file.replace('.tar.gz', '', 1)
        file_no_suffixs.append(file_no_suffix)
        file_work_dirs.append(f'{folder_path}/{file_no_suffix}_work')
        #print(file_work_dirs[-1])
        os.makedirs(file_work_dirs[-1], exist_ok=True)
        print(f'处理第{count_file}个文件：{file}')
        try:
            subprocess.run(['tar', '-xzf', f'{folder_path}/{file}','-C',folder_path], check=True)
            subprocess.run(['cp','-r',f'{folder_path}/{file_no_suffixs[-1]}',file_work_dirs[-1] ], check=True)
            subprocess.run(['rm','-r',f'{folder_path}/{file_no_suffixs[-1]}'], check=True)
        except subprocess.CalledProcessError as e:
            print(f"处理 {file} 失败：{e}")
    return file_work_dirs,file_no_suffixs,count_file

#将XYZ文件按元素组成分类到CHON/other目录
import os
import shutil
from tqdm import tqdm  # 进度条支持

def classify_xyz_files(source_dir, chon_dir, other_dir, allowed_elements={'C', 'H', 'O', 'N'}):
    """
    将XYZ文件按元素组成分类到CHON/other目录
    :param source_dir: XYZ文件源目录
    :param chon_dir: 仅含C/H/O/N的文件存储目录
    :param other_dir: 其他文件存储目录
    :param allowed_elements: 允许的元素集合
    """
    # 创建目标目录
    os.makedirs(chon_dir, exist_ok=True)
    os.makedirs(other_dir, exist_ok=True)
    
    # 遍历源目录[4,5](@ref)
    for root, _, files in os.walk(source_dir):
        # 使用进度条可视化处理过程
        for filename in tqdm(files, desc=f"Processing {os.path.basename(root)}"):
            if not filename.lower().endswith('.xyz'):
                continue
            
            src_path = os.path.join(root, filename)
            try:
                with open(src_path, 'r') as f:
                    # 跳过第一行（原子数）
                    f.readline()
                    # 读取第二行（元素行）
                    elements_line = f.readline().strip()
                
                # 解析元素列表[1,2](@ref)
                elements_str = elements_line.split('[')[-1].split(']')[0]
                elements = {e.strip(" '") for e in elements_str.split(',')}
                
                # 判断元素是否符合条件
                if elements.issubset(allowed_elements):
                    dest_dir = chon_dir
                else:
                    dest_dir = other_dir
                
                # 构建目标路径并复制文件[7,8](@ref)
                dest_path = os.path.join(dest_dir, filename)
                shutil.copy2(src_path, dest_path)
                
            except (IndexError, FileNotFoundError) as e:
                print(f"格式错误或文件不存在: {filename} → {str(e)}")
            except Exception as e:
                print(f"未知错误处理 {filename}: {str(e)}")

#用于将CHON划分为CHON、CHO、CHN、CH、CO、CN、HON、CON等8个子类型
import os
import shutil
from tqdm import tqdm  # 新增进度条库

def classify_chon_subtypes(source_dir, target_base, categories):
    """ 文件分类主函数（含进度条） """
    # 预计算总文件数（用于进度条初始化）[1,4](@ref)
    total_files = sum(
        len([f for f in files if f.endswith('.xyz')])
        for _, _, files in os.walk(source_dir)
    )

    # 初始化进度条[2,4](@ref)
    with tqdm(total=total_files, desc="分类进度", unit="file", dynamic_ncols=True) as pbar:
        for root, _, files in os.walk(source_dir):
            for filename in files:
                if not filename.endswith('.xyz'):
                    continue
                
                src_path = os.path.join(root, filename)
                try:
                    # 读取元素信息
                    with open(src_path, 'r') as f:
                        f.readline()  # 跳过原子数行
                        second_line = f.readline().strip()
                    
                    # 解析元素集合（兼容多种格式）[3](@ref)
                    elements_str = second_line.split('[')[-1].split(']')[0]
                    elements = {e.strip().strip("'\" ").upper() for e in elements_str.split(',')}
                    
                    # 分类匹配（优化循环逻辑）
                    category = next(
                        (cat for cat, req in categories.items() if elements == req),
                        'Other_CHON'
                    )

                    # 构建目标路径
                    relative_path = os.path.relpath(root, source_dir)
                    dest_dir = os.path.join(target_base, category, relative_path)
                    os.makedirs(dest_dir, exist_ok=True)
                    
                    # 复制文件
                    shutil.copy2(src_path, os.path.join(dest_dir, filename))
                    
                except IndexError:
                    tqdm.write(f"格式异常: {filename}")  # 使用tqdm.write避免破坏进度条[4](@ref)
                except Exception as e:
                    tqdm.write(f"处理失败: {filename} → {str(e)}")
                finally:
                    pbar.update(1)  # 更新进度条[1](@ref)
                    pbar.set_postfix_str(f"当前: {filename[:15]}")  # 显示短文件名

import os
import shutil
from tqdm import tqdm
import ase.io

def regroup_by_atomic_sequence(base_dir):
    """
    根据原子序列模式对已分类分子进行二次分组
    :param base_dir: 初始分类目录（如LT10_CHON_SUBTYPES）
    """
    # 遍历所有元素类别目录
    for element_dir in os.listdir(base_dir):
        element_path = os.path.join(base_dir, element_dir)
        if not os.path.isdir(element_path):
            continue
        
        # 初始化原子序列类型计数器
        type_counter = 1
        sequence_map = {}
        
        # 获取所有XYZ文件路径
        all_files = []
        for root, _, files in os.walk(element_path):
            all_files.extend([os.path.join(root, f) for f in files if f.endswith('.xyz')])
        
        # 进度条配置
        with tqdm(total=len(all_files), desc=f"处理 {element_dir}", unit="file") as pbar:
            for src_path in all_files:
                try:
                    # 读取原子序列
                    atoms = ase.io.read(src_path, format='xyz')
                    numbers = atoms.numbers
                    sequence = "".join(map(str, numbers))  # 生成原子序列指纹
                    
                    # 动态创建类型子目录
                    if sequence not in sequence_map:
                        sequence_map[sequence] = f"type{type_counter}"
                        type_counter += 1
                    type_name = sequence_map[sequence]
                    
                    # 构建目标路径
                    dest_dir = os.path.join(element_path, type_name)
                    os.makedirs(dest_dir, exist_ok=True)
                    
                    # 移动文件并保留原始文件名
                    shutil.move(src_path, os.path.join(dest_dir, os.path.basename(src_path)))
                    
                except Exception as e:
                    print(f"处理失败: {os.path.basename(src_path)} → {str(e)}")
                finally:
                    pbar.update(1)

import os
#获得一级子目录的名称列表
def find_subdirectories(directory):
    """
    查找指定目录下的一级子目录

    参数:
    directory (str): 要查找的目录路径

    返回:
    list: 一级子目录的名称列表
    """
    # 使用 os.listdir 获取目录下的所有条目
    all_entries = os.listdir(directory)
    
    # 使用列表推导式和 os.path.isdir 筛选出子目录
    subdirectories = [entry for entry in all_entries if os.path.isdir(os.path.join(directory, entry))]
    for i in range(len(subdirectories)):
                subdirectories[i]=f'{directory}/{subdirectories[i]}'
    return subdirectories

import os

def find_xyz_files(root_dir):
    xyz_files = []
    for root, dirs, files in os.walk(root_dir):     #root为当前目录路径，dirs为该目录下的所有子目录，files为该目录下的所有非目录子文件
        for filename in files:
            if filename.endswith('.xyz'):
                full_path = os.path.join(root, filename)
                xyz_files.append(full_path)
                #print(f"发现文件: {full_path}")  # 实时打印路径
    return xyz_files

import os
import subprocess
import numpy as np
def raw_to_npy(sub_folder_path):
    set000_dir=f'{sub_folder_path}/set.000'
    os.makedirs(set000_dir, exist_ok=True)
    # 读取文本格式的 box.raw
    box_data = np.loadtxt(f'{sub_folder_path}/box.raw')
    np.save(f'{set000_dir}/box.npy', box_data)

    # 读取文本格式的 coord.raw
    coord_data = np.loadtxt(f'{sub_folder_path}/coord.raw')
    np.save(f'{set000_dir}/coord.npy', coord_data)
    try:
        #subprocess.run(['cp','-r',f'{sub_folder_path}/coord.npy',f'{set000_dir}'], check=True)
        #subprocess.run(['cp','-r',f'{sub_folder_path}/box.npy',f'{set000_dir}'], check=True)
        #subprocess.run(['rm','-r',f'{sub_folder_path}/coord.npy'], check=True)
        #subprocess.run(['rm','-r',f'{sub_folder_path}/box.npy'], check=True)
        subprocess.run(['rm','-r',f'{sub_folder_path}/coord.raw'], check=True)
        subprocess.run(['rm','-r',f'{sub_folder_path}/box.raw'], check=True)
    except subprocess.CalledProcessError as e:
            print(f"处理 {sub_folder_path} 失败：{e}")

import os, sys
def extxyz_to_deepmd(Data,dirs,count):
    os.makedirs(f'{dirs}/species_{count:07d}',exist_ok=True)
    deep_sub_dir=f'{dirs}/species_{count:07d}'
    type_map = {1:0, 6:1, 7:2, 8:3}
    element_map = {1:'H', 6:'C', 7:'N', 8:'O'}
    numbers=Data[0][0].numbers
    NATOMS=len(numbers)
    os.system(f'touch {deep_sub_dir}/type_map.raw')
    with open(f'{deep_sub_dir}/type_map.raw', 'w') as f:
        f.write(f'H C N O\n')
    os.system(f'touch {deep_sub_dir}/type.raw')
    with open(f'{deep_sub_dir}/type.raw', 'w') as f:
        for i in numbers:
            f.write(f'{type_map[i]} ')
        f.write('\n')
    os.system(f'touch {deep_sub_dir}/nopbc')
    #with open(f'{deep_sub_dir}/box.raw', 'w') as f:
    #    f.write(f'0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n')

    p_count=0
    os.system(f'touch {deep_sub_dir}/coord.raw')
    for p, mol in enumerate(Data):
        if p > 0 and p_count == 0:
            with open('/home/shuhua01/suchaoxu/work/mdcdnn/data3/test.txt','a') as f:    #这个文件需要提前创建,懒得再写代码(虽然不难实现)
                f.write(f'{count}\n')
            p_count+=1
        numbers = mol[0].numbers
        positions = mol[0].positions
        NAtoms = len(numbers)
        if NAtoms != NATOMS:
            raise ValueError(f"第 {p} 个分子包含 {NAtoms} 个原子，与第一个分子 {NATOMS} 个原子不符")
            #return p
            continue
        else:
            with open(f'{deep_sub_dir}/box.raw', 'a') as f:
                f.write(f'0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n')
            with open(f'{deep_sub_dir}/coord.raw', 'a') as f:
                for i in range(NAtoms):
                    f.write(f'{positions[i][0]:12.6f} {positions[i][1]:12.6f} {positions[i][2]:12.6f} ')
                f.write('\n')
        #sys.stdout.write(f'{p}\r')
    raw_to_npy(deep_sub_dir)

# -*- coding: utf-8 -*-

import os, sys
import numpy as NP
import ase, ase.io, ase.units
from functools import reduce
from tqdm import tqdm 
def deal_xyz_files_main_tqdm(target_folder, deepmd_data_dir,count):
    xyz_file_list=[]
    sub_dirs_2=[]
    #target_folder = "./LT10/LT10_CHON_SUBTYPES"
    if not os.path.exists(target_folder):
        raise FileNotFoundError(f"目录 {target_folder} 不存在")
    sub_dirs_1=find_subdirectories(target_folder)
    for i in range(len(sub_dirs_1)):
        sub_dirs_2.append(find_subdirectories(sub_dirs_1[i]))  
    for i in range(len(sub_dirs_2)):
        xyz_file_list.append([])
        for j in range(len(sub_dirs_2[i])):
            xyz_file_list[i].append([])
            xyz_file_list[i][j]=find_xyz_files(sub_dirs_2[i][j])
    #print(xyz_file_list)
    Data=[]
    total_files=0
    for i in range(len(xyz_file_list)):
        Data.append([])
        for j in range(len(xyz_file_list[i])):
            Data[i].append([])
            total_files+=1
            for k in range(len(xyz_file_list[i][j])):
                #Data[i][j].append([])
                Data[i][j].append(ase.io.read(xyz_file_list[i][j][k], index = ':', format = 'extxyz'))
            #print(Data[i][j])
    err=[]
    #processed_files = 0
    #pbar = tqdm(total=count_files, desc="Processing files", dynamic_ncols=True)  
    #deepmd_data_dir='./LT10/deepmd_data'
    #if not os.path.exists(deepmd_data_dir):
    #    os.makedirs(deepmd_data_dir, exist_ok=True)
    #count=0
    #for i in range(len(Data)):
    #    for j in range(len(Data[i])):
    #        if len(Data[i][j]) !=1:
    #            print(i,j)
    #print(Data[1][1][-1])
    #ers=extxyz_to_deepmd(Data[1][1],deepmd_data_dir,count)
    with tqdm(total=total_files, desc="转换数据为 DeepMD 格式", ncols=100,dynamic_ncols=True) as pbar:
        for i in range(len(Data)):
            #err.append([])
            for j in range(len(Data[i])):
                count+=1
                #pbar.set_postfix_str(f"当前文件: {count}")
                #processed_files+=1
                #print(f"正在处理第 {count} 个文件...")
                #err[i].append([])
                extxyz_to_deepmd(Data[i][j],deepmd_data_dir,count)
                #pbar.update(1)
                pbar.set_postfix_str(f"当前文件: {count}")
                pbar.update(1)
        #pbar.close()
    return count

import os
import subprocess

def process_mian(folder_path):
    #folder_path = './test'
    CATEGORIES = {
            'CHON': {'C', 'H', 'O', 'N'},  # 四元素组合
            'CHO': {'C', 'H', 'O'},         # 碳氢氧组合
            'CHN': {'C', 'H', 'N'},
            'CH': {'C', 'H'},
            'CO': {'C', 'O'},
            'CN': {'C', 'N'},
            'HON': {'H', 'O', 'N'},
            'CON': {'C', 'O', 'N'}
    }
    deepmd_data_dir = f'{folder_path}/deepmd_data'
    if not os.path.exists(deepmd_data_dir):
        os.makedirs(deepmd_data_dir, exist_ok=True)
    count=0
    err=[]
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"目录 {folder_path} 不存在")
    else:
        print(f"目录 {folder_path} 存在,执行解压操作")
        work_dir,file_no_suffixs,count_file=extract_tar_gz(folder_path)
        print(f"解压完成，工作目录为{work_dir}，文件总数为{count_file}")
        for i in range(len(work_dir)):
            print(f"Processing {work_dir[i]}")
            SOURCE_DIR = f'{work_dir[i]}/{file_no_suffixs[i]}'
            CHON_DIR = f'{work_dir[i]}/{file_no_suffixs[i]}_CHON'
            OTHER_DIR = f'{work_dir[i]}/{file_no_suffixs[i]}_OTHER'
            TARGET_BASE = f'{work_dir[i]}/{file_no_suffixs[i]}_CHON_SUBTYPES'      
            classify_xyz_files(SOURCE_DIR, CHON_DIR, OTHER_DIR)
            classify_chon_subtypes(CHON_DIR, TARGET_BASE, CATEGORIES)
            regroup_by_atomic_sequence(TARGET_BASE)
            #print(TARGET_BASE)
            count=deal_xyz_files_main_tqdm(TARGET_BASE,deepmd_data_dir,count)
            i+=1
        
folder_path = '/home/shuhua01/suchaoxu/work/mdcdnn/data3/CHEMBL_selected_data'
process_mian(folder_path)    


