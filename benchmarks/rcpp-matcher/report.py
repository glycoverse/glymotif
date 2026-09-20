from pathlib import Path
import csv,statistics,collections
p=Path('benchmarks/rcpp-matcher')
r=list(csv.DictReader((p/'timings.csv').open())); groups=collections.defaultdict(list)
for x in r:groups[(x['workload'],x['output'],x['engine'])].append(x)
cases=[('scalar','have'),('scalar-positive','have'),('512x25','have'),('512x25','count'),('512x25','match'),('duplicates-4096x25','have'),('full-corpus','have')]
labels={'scalar':'单对结构（阴性）','scalar-positive':'单对结构（阳性）','512x25':'512 结构 × 25 motif','duplicates-4096x25':'4,096 行（64 个独立结构）× 25 motif','full-corpus':'8,009 结构 × 13 motif'}
lines=['| 工作负载 | 输出 | 当前 API，ms | C++ 直接入口，ms | 保留 R 参数检查，ms | 直接入口加速 | 保留检查加速 |','|---|---|---:|---:|---:|---:|---:|']
for case in cases:
 med={e:statistics.median(float(x['elapsed_s']) for x in groups[case+(e,)])*1000 for e in ['current','cpp_direct','cpp_validated']}
 lines.append(f'| {labels[case[0]]} | {case[1]} | {med["current"]:.3f} | {med["cpp_direct"]:.3f} | {med["cpp_validated"]:.3f} | {med["current"]/med["cpp_direct"]:.2f}× | {med["current"]/med["cpp_validated"]:.2f}× |')
parity=list(csv.DictReader((p/'parity.csv').open()))
checks=sum(int(x['pairs']) for x in parity)
text=f'''# glycan_structure → C++ → 结果：端到端匹配实验

日期：2026-09-20。结论：**普通、已定位的糖链匹配值得把整条处理路径迁入 C++**。
全量 8,009 × 13 存在性查询中，直接原生入口约快 26.8 倍；保留现有公开 API
参数准备后的保守版本约快 5.0 倍。收益来自移走 R 层图提取、兼容性计算、批量调度
和输出整理；现有 VF2 已经是 C++，本实验保留了同一个 Boost 搜索实现。

**这不是全功能替换已完成的声明**：564 个带浮动信息的结构没有进入性能数据。
浮动结构的候选定位枚举、冲突约束和跨定位汇总尚未移植；原型对此明确报错。
生产 API、NAMESPACE、依赖和现有匹配代码均未修改。

## 范围和计时边界

- 本地 `glydb::glydb_structures()`：8,573 个结构，8,009 个无浮动信息，564 个有浮动信息。
- 随机种子 20260920；从普通结构抽取 512 个，用 13 个边界 motif 加 12 个真实 motif
  组成 25 个 motif；全量测试使用前 13 个 motif。输入键保存在 `inputs.rds`。
- 输入是已经构造好的 `glycan_structure` 对象；三条路线收到相同对象和选项。
- 每次调用包含：对象读取、唯一结构索引、图信息提取、C++ 图与兼容性矩阵构建、
  搜索、映射去重、NA/重复行恢复和 R 结果创建。没有在计时外缓存图 profile 或匹配矩阵。
- 不计：字符串解析、库加载、Rcpp 编译，以及从包命名空间取得静态残基字典。
  每次把该字典转换为 C++ 数据结构的成本已计入。
- `cpp_direct` 是已类型化对象的原生入口，**不等于完整公开 API 替换**；
  `cpp_validated` 每次先调用当前 `prepare_motif_args()`，再进入 C++。
  后者更适合作为现有 API 接入后的收益估计。
- 当前基线用 `pkgload::load_all()` 加载相邻 glymotif checkout，版本 0.19.0.9000，
  commit `3dce7f593da6c655102806490b821c491ef01f81`；不是旧版纯 R 算法。
  glyrepr 使用已安装的 1.0.0，glydb 0.6.0.9000，igraph 2.3.3，Rcpp 1.1.2，BH 1.90.0-1。

## 时间结果

下面是单次调用的墙钟时间中位数，单位 **毫秒**。加速倍数是两条路线中位数之比。
通常 7 轮，全量 5 轮；单对查询每轮执行 100 次后取均值。每轮随机路线顺序，
每项测量前执行 GC；各路线先暖机。CPU 时间和每轮原始结果也已保存。

{chr(10).join(lines)}

这些是本机、上述 workload 的测量，不是跨硬件保证，也不代表所有 motif 复杂度。
常规数据上存在性、计数和完整映射均有明显改善。单对调用保留参数检查后的收益
约为两倍，说明批量入口比反复调用标量入口更适合发挥原生实现的收益。
重复输入会被两种实现去重；4,096 行不应解释为 4,096 个独立结构搜索。

## 一致性证据

- 共 **{checks:,} 次结构对 × 输出类型比较**，{len(parity)} 个 workload/输出检查条目，
  全部使用 `identical()` 通过。数字含同一结构对在不同模式、选项和输出上的重复检查，
  不是独立生物结构数量。
- 27 个边界结构，加入重复项和 NA，遍历 strict/lenient、四种 alignment、
  strict_sub 两值、ignore_linkages 两值；另测 degree mask 全 FALSE/全 TRUE 和空 glycan。
- 包含具体/通用/混合残基、环型、显式构型、内置修饰、模糊及重复取代基、
  未知或多候选连接、端基、分支对称性、alditol、名称及缺失值。
- 512 × 25 真实/边界混合矩阵覆盖两种 mode 和四种 alignment，三种输出逐项一致。
- 全部 8,009 × 13 的 have/count/match 也分别完全一致；完整映射包括列表结构、
  整数节点编号、顺序、同构映射去重、名称和 NA 返回形式。
- 被排除的 564 个浮动结构逐个检查，全部明确报“不支持浮动定位”，没有静默降级。
- 这是针对有效结构对象的差分验证，不是对所有可能输入、异常路径的形式证明。

## 实现和生产化边界

`cpp_structure_match()` 直接读取结构向量的图缓存、节点和边数组，调用期间没有
R 函数回调。残基类别字典来自 glyrepr 包数据，匹配规则本身在 C++ 执行。
原型复用了 `glymotif/src/vf2.cpp` 的 Boost VF2，保持搜索引擎一致。

最重要的两个后续工作是：

1. **图接口稳定性**。原型直接读 igraph 的内部 10 槽布局。setup 明确锁定
   igraph 2.3.3，C++ 有布局和树形基本检查。这适合可行性实验，不能当作稳定 ABI。
   正式实现应建立稳定的原生桥接，或让 glyrepr 持有自己定义的紧凑图数据；
   若改变提取机制，需要重新测量端到端成本。
2. **浮动语义完整性**。需移植候选父节点、位点占用冲突、环/连通约束、256 组合上限、
   所有/任一定位量词以及计数/映射聚合，并对现有浮动测试逐项验证。
   本实验没有给这些功能的加速预测。

建议以“保留 R 公开参数层 + 批量 C++ 匹配后端”推进普通结构路径；
这已实测有约 5–6 倍的批量收益。完整的全 C++ 公共接口和浮动语义仍需要另外验证。

## 可复跑文件

从 glymotif 仓库根目录运行：

```sh
Rscript benchmarks/rcpp-matcher/run.R
Rscript benchmarks/rcpp-matcher/supplement.R
```

`setup.R` 是加载和接口包装，`matcher.cpp` 是独立原型，`timing.R` 是交错计时器。
`parity.csv`、`timings.csv`、`summary.csv` 是核验和计时结果；`inputs.rds` 固定输入键，
`sessionInfo.txt`、`environment.txt` 和 `source-sha256.json` 保存环境与源文件证据。
实验位于 glymotif 的 benchmarks/rcpp-matcher 目录。
'''
(p/'REPORT.md').write_text(text)
# Numeric companion, including ranges, CPU timing and medians.
with (p/'statistics.csv').open('w') as f:
 w=csv.writer(f,lineterminator="\n");w.writerow(['workload','output','engine','reps','median_s','min_s','max_s','median_cpu_s'])
 for k,rows in groups.items():
  v=[float(x['elapsed_s']) for x in rows];cpu=[float(x['cpu_s']) for x in rows]
  w.writerow([*k,len(v),statistics.median(v),min(v),max(v),statistics.median(cpu)])
print('parity',checks,'rows',len(parity),'timing rows',len(r))
