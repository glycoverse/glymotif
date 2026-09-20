from pathlib import Path
import csv, statistics, collections
p=Path('benchmarks/rcpp-floating-matcher')
a=list(csv.DictReader((p/'localization-audit.csv').open()))
f=list(csv.DictReader((p/'fixture-parity.csv').open()))
t=list(csv.DictReader((p/'timings.csv').open()))
groups=collections.defaultdict(list)
for r in t: groups[(r['workload'],r['output'],r['strict_floating'],r['engine'])].append(r)
lines=['| 数据／输出 | 当前 API（秒） | 直接 C++（秒） | 保留参数检查（秒） | 保留检查加速 |', '|---|---:|---:|---:|---:|']
for work,kind,sf in [('32x6','have','TRUE'),('32x6','count','TRUE'),('32x6','match','TRUE'),('32x6-any','have','FALSE'),('503x6','have','TRUE')]:
    med={e:statistics.median(float(x['elapsed_s']) for x in groups[(work,kind,sf,e)]) for e in ['current','cpp_direct','cpp_validated']}
    lines.append(f'| {work} / {kind} / strict_floating={sf} | {med["current"]:.3f} | {med["cpp_direct"]:.3f} | {med["cpp_validated"]:.3f} | {med["current"]/med["cpp_validated"]:.1f}× |')
variants=sum(int(x['variants']) for x in a if x['result']=='identical')
checks=sum(int(x['pairs']) for x in f)+503*6*6
text=f'''# C++ 浮动结构匹配：核验与端到端 benchmark

日期：2026-09-20。在独立原型中，**浮动定位、子结构匹配和跨定位汇总已全部在 C++ 执行**，
从原始 `glycan_structure` 读取图和浮动元数据，返回 R 结果。定位和匹配期间没有 R 函数回调。
当前生产 API 未修改。此前普通结构实验保存在相邻 `rcpp-matcher` 目录。

## 数据范围与一致性

本地 glydb 共 8,573 个结构，其中 564 个有浮动信息：560 个包含浮动糖链组件，
5 个包含浮动取代基，二者有 1 个重叠。

- **503 个**在现有 256 种原始候选组合上限内，产生 **{variants:,} 个有效定位**。
  每个定位的边、连接、取代基、原始节点编号以及父节点选择顺序，与 R 完全一致。
- **61 个**超过上限；C++ 和现有 R API 均报错。本实验不提高该上限，
  也不把这些结构计入性能分母。
- 503 × 6 motif 的存在性、计数和完整映射，在两种浮动模式下全部一致。
  对每个原始结构，R 参考先定位，再调用现有 `.match_motif_single()`，用现有 R 汇总规则
  得到期望结果；为减少验证时间，这些中间结果只在正确性核对中复用。
  另用固定前 12 个结构核对六种输出与公开 API 完全一致。
- 13 个边界结构加重复项、NA（15 行）× 12 motif，测试 strict/lenient、四种 alignment、
  strict_floating 两值、链接与取代基选项、degree mask 全 TRUE/FALSE。
  {len(f)} 个配置／输出检查、{sum(int(x['pairs']) for x in f):,} 次结构对／输出比较全部通过。
- 上述边界加全量匹配共 **{checks:,} 次结构对／输出比较**，使用 `identical()`；
  含参数和输出的重复，不是独立结构数量。不将定位图检查混入此数字。
- 单独验证组件互相成环、糖基／取代基占用相同碳位点时拒绝，
  恰好 256 种组合时保留全部定位，512 种组合时报错。
- 未定位的浮动 motif 仍被拒绝，保持当前公开 API 的边界。

## 端到端结果

以下是三轮交错计时的墙钟中位数，单位 **秒**。输入为同一批已经构造好的浮动结构对象。
每次调用重新读取图、枚举定位、构造匹配图、搜索、汇总、分配并返回 R 结果；
**没有在计时外预定位或缓存图 profile**。不计字符串解析、包加载和编译。
字典转为 C++ 数据结构的成本在计时内。

{chr(10).join(lines)}

32 个结构按有效定位数分层选取，6 个 motif 与全量一致；全量是全部 503 个可处理的
浮动结构 × 6 motif，共 3,018 对查询。严格浮动模式为默认。
各路线先暖机，每轮随机顺序，每次测量前 GC，计时包含调用内 GC；同时记录 CPU 时间。

`cpp_direct` 是有效结构对象的专用入口。`cpp_validated` 每次额外运行当前
`prepare_motif_args()`，更接近保留 R 公开参数层后接入 C++ 的收益。
从代码路径看，收益包括省去 R 层反复生成／校验 igraph、tibble、兼容性计算和逐定位调度；
此外，存在性查询在严格模式遇到首个失败定位、任一模式遇到首个成功定位时短路。
现有 R 路径先匹配所有定位再汇总；计数和映射没有这种跨定位短路优化。
VF2 仍使用原型继承的相同 Boost 实现，不能解释为 VF2 搜索算法本身有同等加速。

这是单机、这些结构和 motif 的测量，不是所有复杂度或硬件上的速度保证。
全部计时在正确性任务完成后独立执行。

## C++ 语义

1. 对每个浮动组件，读取显式候选父节点；未指定时，候选为自身组件以外全部节点。
   浮动取代基未指定父节点时可选全部节点。
2. 按 R `expand.grid()` 的顺序枚举，先检查原始笛卡尔积是否超过 256。
3. 检查组件父链最终连到主树、没有环；将固定连接、固定取代基和本次浮动选择的
   已知碳位点作为候选集合，以二分匹配判断是否有无冲突位点分配。
   未知位置不占用已知位点，多候选位置表示可选择的集合。
4. 在 C++ 中追加连接或取代基，保留原节点编号；定位之间不做结构同一性去重。
5. 存在性按所有／任一定位汇总，并在结果已确定时短路；计数取最小／最大值。
   映射取各定位映射的并集，按原向量去重并保留首次出现顺序。
   公开 `match_motifs()` 不接受 strict_floating；它始终使用映射并集语义。

## 限制与复现

这是实验后端，**不是已合入的生产实现**。仍直接读取 igraph 2.3.3 的内部布局，
setup 有版本检查；正式接入需要稳定图桥接或 glyrepr 自有的紧凑图表示。
输入前提是有效的 glyrepr 对象，未覆盖任意手工伪造对象、所有错误信息样式、
大于当前上限的无界定位，或所有可能的 motif 组合。

基线为 glymotif 0.19.0.9000，checkout `4e239db1263799cc4b2732fd7239cb8249e139d6`；
glyrepr 使用已安装 1.0.0，igraph 2.3.3。完整环境和源文件指纹另存。

从 glymotif 根目录按 [README](README.md) 顺序运行审计、边界测试、全量匹配与计时脚本。
`localization-audit.csv`、`fixture-parity.csv`、`corpus-parity.csv`、`guard-audit.csv`
保存核验结果，`timings.csv` 为逐轮原始计时，`statistics.csv` 还包含范围与 CPU 中位数。
'''
(p/'REPORT.md').write_text(text)
with (p/'statistics.csv').open('w') as fp:
    w=csv.writer(fp,lineterminator='\n');w.writerow(['workload','output','strict_floating','engine','reps','median_s','min_s','max_s','median_cpu_s'])
    for k,v in groups.items():
        times=[float(x['elapsed_s']) for x in v];cpu=[float(x['cpu_s']) for x in v]
        w.writerow([*k,len(v),statistics.median(times),min(times),max(times),statistics.median(cpu)])
print('Report generated; comparisons:',checks,'localizations:',variants)
