"""
带有循环的化工系统序贯求解过程: 子系统识别、撕裂物流切割、确定求解序列
给定表示化工流程的邻接矩阵，返回撕裂物流和求解顺序
使用方法可查看测试文件test_sequence1.ipynb

参考文献:
[1] 邱彤,戴一阳,赵永臣.化工过程模拟——理论与实践[M]. 北京:化学工业出版社, 2019.

"""

import copy

def calculate_sequence(blocks, matrix):
    """ 
    输入各单元模块编号blocks和邻接矩阵matrix, 返回各单元求解顺序和撕裂物流编号
    """

    block_index, stream_index = calculate_sequence_index(matrix)

    for i in range(len(block_index)):
        if isinstance(block_index[i], list):
            for j in range(len(block_index[i])):
                block_index[i][j] = blocks[block_index[i][j]]
        else:
            block_index[i] = blocks[block_index[i]]
    
    for subsystem in stream_index:
        for stream in subsystem:
            stream[0] = blocks[stream[0]]
            stream[1] = blocks[stream[1]]

    return block_index, stream_index


def calculate_sequence_index(matrix):
    """ 
    输入邻接矩阵matrix, 先获取各单元的计算顺序的索引,其中不可分割子系统视为一个整体
    再对各个不可分割子系统进行切割, 获取撕裂物流的索引以及切割后子系统的计算顺序的索引
    """

    # 各单元的计算顺序(索引形式)
    sequence = index_sequence(matrix)   

    tear_streams = []
    for i in range(len(sequence)):
        blocks = sequence[i]
        if isinstance(blocks, list):
            # 构建不可分割子系统的邻接矩阵
            submatrix = [[0 for i in range(len(blocks))] for i in range(len(blocks))]
            for row in range(len(blocks)):
                for col in range(len(blocks)):
                    submatrix[row][col] = matrix[blocks[row]][blocks[col]]
            
            # 获取撕裂物流, 并添加到集合中
            tear_stream = greedy_tear_streams(blocks, submatrix)     
            tear_streams.append(tear_stream) 
            
            # 获取撕裂物流的索引[start, end], 将撕裂物流断开, 即submatrix[start, end]=0
            tear_stream_index = greedy_tear_stream_index(submatrix) 
            for index in tear_stream_index:
                start = index[0]
                end = index[1]
                submatrix[start][end] = 0
            # 获取子系统分割后的计算顺序
            sequence[i] = block_sequence(blocks, submatrix)

    return sequence, tear_streams
    

def greedy_tear_streams(blocks, matrix):
    """
    使用贪心算法切割不可分割子系统, matrix为不可分割子系统的邻接矩阵
    返回撕裂物流的起始单元编号[[start, end], [start, end]]
    start为开始单元编号, end为终止单元编号
    """

    tear_stream_index = greedy_tear_stream_index(matrix)
    for index in tear_stream_index:
        index[0] = blocks[index[0]]
        index[1] = blocks[index[1]]
    
    return tear_stream_index


def greedy_tear_stream_index(matrix):
    """
    使用贪心算法切割不可分割子系统, matrix为不可分割子系统的邻接矩阵
    返回撕裂物流的起始单元索引[[start, end], [start, end]]
    start为开始单元索引, end为终止单元索引
    """

    cycles = find_cycles(matrix)    # 所有的回路
    number_of_cycles = len(cycles)  # 回路数量

    # 获取所有流股, 流股以[start, end]表示
    streams = []   
    for cycle in cycles:
        for start in range(len(cycle)):
            end = (start+1) % len(cycle)
            stream = [cycle[start],cycle[end]]
            if stream not in streams:
                streams.append(stream)
    number_of_streams = len(streams)  # 流股的数量

    # 构建回路矩阵
    cycles_matrix = [[0 for i in range(number_of_streams)] for i in range(number_of_cycles)]
    # 遍历所有回路和所有流股, 若流股在回路中则赋值为1
    for row in range(number_of_cycles):
        cycle = cycles[row]
        cycle_streams = []    # 该回路的所有流股
        # 遍历回路
        for start in range(len(cycle)):
            end = (start+1) % len(cycle)
            stream = [cycle[start],cycle[end]]
            cycle_streams.append(stream)
        # 遍历流股
        for col in range(number_of_streams):
            stream = streams[col]
            if stream in cycle_streams:
                cycles_matrix[row][col] = 1
    # 使用贪心算法获取撕裂物流集合
    tear_streams = []   # 
    while len(cycles) > 0:
        max_cut = 0   # 流股所切割回路的数量的最大值
        max_col = 0   # 上述流股所对应的索引
        max_cut_row = []  # 上述流股所切割的回路的索引
        for col in range(len(streams)):
            cut = 0
            cut_row = []
            for row in range(len(cycles)):
                if cycles_matrix[row][col] != 0:
                    cut += 1 
                    cut_row.append(row)
            if cut > max_cut :
                max_cut = cut
                max_col = col
                max_cut_row = cut_row
        # 将撕裂物流从流股集合中移除，并添加到撕裂物流集合中
        tear_stream = streams.pop(max_col)    
        tear_streams.append(tear_stream)
        # 从回路矩阵中移除撕裂物流所在列，并移除被切割的回路所在行
        for row in range(len(cycles)): 
            cycles_matrix[row].pop(max_col)
        max_cut_row.sort(reverse=True)
        for row in max_cut_row:
            cycles.pop(row)
            cycles_matrix.pop(row)    
    return tear_streams

def block_sequence(blocks, matrix):
    """ 
    输入模块列表blocks和邻接矩阵matrix, 返回各模块的求解顺序
    """

    index = index_sequence(matrix)
    for i in range(len(index)):
        if isinstance(index[i], list):
            for j in range(len(index[i])):
                index[i][j] = blocks[index[i][j]]
        else:
            index[i] = blocks[index[i]]
    return index


def index_sequence(matrix):
    """ 
    根据回路搜索法, 输入一个邻接矩阵, 返回流程的求解序列, 回路视为一个不可分割子系统
    结果返回各单元的索引
    """
    
    copy_matrix = copy.deepcopy(matrix)  # 复制一个矩阵进行操作, 避免破坏原矩阵
    n = len(copy_matrix)     # 单元数量

    before = []  # 先求解的单元
    after = []   # 后求解的单元
    # 剩余未处理的单元，该变量用于保存各单元在原始矩阵中的索引。
    # 因为以下操作涉及矩阵的行和列的剔除操作，导致各单元在矩阵中的索引发生变化
    remains = list(range(n))  
    
    while len(copy_matrix) > 0:
        delete_blocks = []  # 需要剔除的单元
        n = len(copy_matrix)     # 单元数量
        # 查找全为0的列
        for col in range(n):
            flag = True
            for row in range(n):
                if copy_matrix[row][col] != 0:
                    flag = False
                    break
            if flag:
                before.append(remains[col])
                delete_blocks.append(col)
              
        # 查找全为0的行
        for row in range(n):
            flag = True
            for col in range(n):
                if copy_matrix[row][col] != 0:
                    flag = False
                    break
            if flag:
                if row not in delete_blocks:
                    after.insert(0, remains[row])
                    delete_blocks.append(row)

        # 根据索引从大到小剔除上述的单元所在的行和列
        delete_blocks.sort(reverse=True)
        for i in delete_blocks:
                delete(copy_matrix, remains, i)

        # 根据邻接矩阵, 找到一条回路
        cycle =find_cycle(copy_matrix)

        # 合并回路的节点, 构造新的邻接矩阵
        if len(cycle) != 0:
            merge(copy_matrix, remains, cycle)

    before.extend(after)  # 合并先后两个求解序列
    return before   


def delete(matrix, blocks, i):
    """ 从邻接矩阵中剔除索引为i的行和列, 并且从单元模块数组中剔除索引为i的单元"""

    blocks.pop(i)       # 从单元模块数组中剔除索引为i的单元
    length = len(matrix)   # 邻接矩阵的长度
    for row in range(length):
        matrix[row].pop(i)   #  从邻接矩阵中剔除索引为i的列
    matrix.pop(i)       #  从邻接矩阵中剔除索引为i的行


def merge(matrix, blocks, cycle):
    """
    cycle为回路中各单元的索引, 将回路中的各单元合并为新的虚拟单元, 构建新的邻接矩阵
    """

    cycle_blocks = []  # 代表回路中的各单元
    for i in cycle:
        if isinstance(blocks[i],list):        # block[i]可能是个回路
            for block in blocks[i]:
                cycle_blocks.append(block)
        else:
            cycle_blocks.append(blocks[i])
    blocks.append(cycle_blocks)   # 添加合并的虚拟单元

    length = len(matrix)   # 原邻接矩阵的长度
    # 在邻接矩阵末尾添加新的一列
    for row in range(length):
        sum = 0
        for i in cycle:
            sum = sum | matrix[row][i]   # 布尔加法(或运算)
        matrix[row].append(sum)
    
    # 在邻接矩阵末尾添加新的一行
    matrix.append([0]*(length+1))
    for col in range(length):
        sum = 0
        for i in cycle:
            sum = sum | matrix[i][col]   # 布尔加法(或运算)
        matrix[length][col] = sum
    
    # 按照索引从大到小，剔除回路各单元所在行和列
    cycle.sort(reverse=True)
    for i in cycle:
        delete(matrix, blocks, i)

    
def find_cycle(matrix):
    """ 输入一个邻接矩阵, 返回一个回路(各单元以索引返回) """

    cycles = find_cycles(matrix)
    if (len(cycles) == 0):
        return []
    else:
        return cycles[0]


def find_cycles_blocks(blocks, matrix):
    """ 输入一个邻接矩阵, 返回所有回路(各单元以编号表示) """
    cycles = find_cycles(matrix)
    for cycle in cycles:
        for i in range(len(cycle)):
            cycle[i] = blocks[cycle[i]]
    return cycles
    
            
def find_cycles(matrix):
    """
    输入一个邻接矩阵, 返回所有回路(各单元以索引表示)
    参考代码: https://www.iteye.com/blog/128kj-1717469
             https://blog.csdn.net/dizhuangrou5770/article/details/101191909
    """

    n = len(matrix)     # 单元数量
    cycles = []         # 储存回路的列表

    # 遍历每个单元进行深度优先搜索
    for i in range(n):
        visited = [0]*n     # 单元访问状态，初始0为未访问
        trace = []          # 从出发点到当前节点的轨迹
        dfs(matrix, n, visited, trace, cycles, i)  

    return cycles


def dfs(matrix, n, visited, trace, cycles, i):
    """ 深度优先搜索 """

    if visited[i] :    
        if i in trace:   
            j = trace.index(i)
            cycle = []
            while j < len(trace): 
                cycle.append(trace[j])
                j += 1  
            add_cycle(cycles, cycle)
            return    
        return    
    visited[i]=1  
    trace.append(i)  
        
    for v in range(n):
        if matrix[i][v] == 1:
            dfs(matrix, n, visited, trace, cycles, v)   
    trace.pop()   


def add_cycle(cycles, cycle):
    """
    遍历结果集cycles, 检查要添加的回路cycle是否在cycles中, 若无, 则添加 
    """
    if len(cycles) == 0:
        cycles.append(cycle)
    else:
        set1 = set(cycle)
        for i in cycles:
            if set1 == set(i):
                return
        cycles.append(cycle)


if __name__ == "__main__":

    matrix = [[0, 1, 0, 1, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0]]
    print(find_cycle(matrix))

    blocks = list(range(6))
    delete(matrix, blocks,1)
    print(blocks)
    print(matrix)
    print(isinstance(blocks,list))

    matrix = [[0,1,1,0,0,0,0],  
            [0,0,0,1,0,0,0],  
            [0,0,0,0,0,1,0],  
            [0,0,0,0,1,0,0],  
            [0,0,1,0,0,0,0],  
            [0,0,0,0,1,0,1],  
            [1,0,1,0,0,0,0]]
    print(find_cycle(matrix))

    matrix = [[0,1,1],
            [0,0,0],
            [0,1,0]]
    print(find_cycle(matrix)) 

    print(50*"*")
    blocks = list(range(1,10))
    matrix = [
            [0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0]]
    print(block_sequence(blocks, matrix))

    blocks = list(range(1,5))
    matrix = [
            [0, 1, 0, 0],
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 1, 0]]
            
    print("求解顺序为: ")  
    print(block_sequence(blocks, matrix))

    print(50*"*")
    blocks = [1,2,3,4]
    matrix = [
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [1, 0, 0, 1],
            [1, 0, 1, 0]]
    print(find_cycles(matrix))
    print(greedy_tear_streams(blocks, matrix))

    print(50*"*")
    matrix = [
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [1, 1, 0, 1],
            [1, 1, 0, 0]]
    print(find_cycles(matrix))
    print(greedy_tear_streams(blocks, matrix))


    blocks = list(range(1,10))
    matrix = [
            [0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0]]
    print("求解顺序为: ")
    print(calculate_sequence(blocks, matrix)[0])
    print("撕裂物流")
    print(calculate_sequence(blocks, matrix)[1])
