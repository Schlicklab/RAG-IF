# -*- coding: UTF-8 -*-

# This function is designed to check the BPSEQ file 
# and it can be called within "runRNAinverse.py" script. 
#
# Sometimes, the RNAfold program will fail to generate the centroid   
# structure (secondary structure), it simply says "...." as the output.
#
# When above output is converted into BPSEQ file, the third column  
# for all NTs will be zero, indicating no pairing. However, the TreeGraph   
# code will crash in case of such BPSEQ files. Thus, we would like to   
# avoid feeding TreeGraph code such BPSEQ file as an input to 
# eliminate errors.
#
# This function read in a BPSEQ file as the input, it returns a signal to   
# tell the main program whether to skip running TreeGraph or not.   
#
# Y. Tao - 2018/07/16
#
import sys,time



class ShowProcess():
    """ 
    显示处理进度的类
    调用该类相关函数即可实现处理进度的显示
    """
    i = 0 # 当前的处理进度
    max_steps = 0 # 总共需要处理的次数
    max_arrow = 50 #进度条的长度
    infoDone = 'done'

    # 初始化函数，需要知道总共的处理次数
    def __init__(self, max_steps, infoDone = 'Done'):
        self.max_steps = max_steps
        self.i = 0 
        self.infoDone = infoDone

    # 显示函数，根据当前的处理进度i显示进度
    # 效果为[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]100.00%
    def show_process(self, i=None):
        if i is not None:
            self.i = i 
        else:
            self.i += 1
        num_arrow = int(self.i * self.max_arrow / self.max_steps) #计算显示多少个'>'
        num_line = self.max_arrow - num_arrow #计算显示多少个'-'
        percent = self.i * 100.0 / self.max_steps #计算完成进度，格式为xx.xx%
        process_bar = '[' + '>' * num_arrow + '-' * num_line + ']'\
                      + '%.2f' % percent + '%' + '\r' #带输出的字符串，'\r'表示不换行回到最左边
        sys.stdout.write(process_bar) #这两句打印字符到终端
        sys.stdout.flush()
        if self.i >= self.max_steps:
            self.close()

    def close(self):
        print('')
        print(self.infoDone)
        self.i = 0







def chkBPSEQ( inpf ): 

    toskip = 'y'

    with open( inpf ) as pointer:
         for line in pointer:
             if len(line) > 1:
                thirdC = line.split()[2] 
                if thirdC != '0':
                   toskip = 'n'

    #
    return toskip


# This function is designed to add another verification for candidate  
# sequence, besides comparing the Graph_ID only.   
# 
# Such verification is to check whether the starting and ending sequence    
# pairs or not, for at least three residues. 
# 
# If all three residues are not paired, then this candidate might not be  
# a good one, even though its Graph_ID is correct.
#
# Y. Tao - 2018/07/18
#
def chkTermi( inpf ):

    toskip = 'n'

    nline = 0 
    with open( inpf ) as pointer:
         for line in pointer:
             if len(line) > 3:
                nline = nline + 1   
   
    ifpair1 = 'n'
    ifpair2 = 'n'
    ifpair3 = 'n'

    counter = 0
    with open( inpf ) as pointer:
         for line in pointer:
             if len(line) > 3:
                counter = counter + 1
                if counter == 1:
                   partner = int(line.split()[2])
                   if partner == nline:
                      ifpair1 = 'y'
                if counter == 2:
                   partner = int(line.split()[2])
                   if partner == (nline-1):
                      ifpair2 = 'y'
                if counter == 3:
                   partner = int(line.split()[2])
                   if partner == (nline-2):
                      ifpair3 = 'y'
    #
    if (ifpair1 == 'n') and (ifpair2 == 'n') and (ifpair3 == 'n'):
       #
       toskip = 'y'

    #
    return toskip










