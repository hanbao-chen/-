# 近岸 H3 小海域热点识别中的证据与代价

《海洋数据挖掘技术》课程研究作业　**朱宸浩　22634107**

公开仓库：https://github.com/hanbao-chen/-

- [课程研究报告 PDF](report/course_report_22634107.pdf)
- [课程研究报告 Markdown](report/course_report_22634107.md)
- [完整人机对话 Markdown](conversations/full_conversation.md)
- [对话导出说明与条数核查](conversations/export_metadata.json)

## 研究内容

基于实际 H3 r7 网格、人工岸线和已知真值的模拟船位计数，检查三件事：面积收缩减少误报是否伴随小热点漏检；把重复船位当成独立证据会怎样夸大确认率；在公平误报约束下，连续观测与“待观察”需要付出多少等待和未决代价。

本项目的贡献是具体实验设计、机制分析和评价口径，不宣称发明统计收缩、泊松检验或三分类。所有结果均为模拟结果。相关结构在重复观测场景中由生成器直接给出，未实现真实 AIS 中的相关结构识别。

## 工程目录

```text
report/                             署名报告与图表
conversations/                      实际完整对话、原始结构化答复与导出元数据
experiments/h3_coastal_experiment/   第一阶段原始代码、方案、计数矩阵、图表及指纹
experiments/threshold_diagnostic/    追加门槛等价诊断的原始脚本与输出
experiments/h3_continuous_checks/    连续观测、相关性与待观察实验
run_project.py                      无覆盖核查和隔离复现入口
requirements.txt                    固定实验依赖
```

不公开分发老师课件或同学参考作业原件。公开的课程来源只保留文件名称和核对过的页码。对话保留实际的选题、理解困难、追问、修正和工程协作过程。

## 快速核查

Python 3.11 环境下，在本项目根目录执行：

```powershell
python run_project.py verify
```

此核查只依赖标准库，验证原始实验文件指纹及三状态比例守恒。原始结果不会被覆盖。

## 完整复现

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python run_project.py reproduce --output reproductions/my_run
```

macOS/Linux可将`.venv\Scripts\python`替换为`.venv/bin/python`。为避免NumPy 1.24.3的Python版本兼容问题，推荐Python 3.11。新目录必须不存在。程序复制原始源码，在新目录运行两阶段实验，最后将数值摘要与归档结果比较。包含造数据、统计、保存文件和绘图，耗时随机器变化；原始第二阶段墙钟约33秒，不能换算为真实监测时长。

图表使用Microsoft YaHei，其他平台可替换字体；字体差异不影响计数和数值摘要。原始指纹清单采用Windows相对路径，核查入口会规范路径分隔符。

## 数据与评价纪律

- 第一阶段：两个场景各2,000次模拟快照。
- 第二阶段：3种密度场景×2种时间结构×2,000条序列，每条120轮，在7个预定轮次检查。
- 所有方法在同一场景内读取相同数据。第二阶段从首帧开始，不调用第一阶段20帧历史数据；背景0.5为共同已知的模拟参数。
- 截止仍待观察的真实热点计入未检出，不能从分母排除。
- 全图至少一次误报的序列比例，与单格误报率、误报格数分开报告。
- 固定120轮和逐步检查按相同全图错误预算比较；原始密度二分类不具有同样的风险保证。
- 已知分组校正是理想信息条件参照，不是实际AIS自动校正性能。

完整复现由AI协作执行；不将其表述为学生已经独立完成复现。最终报告如实列出人机分工。

## 报告重建

PDF由同目录Markdown生成。安装requirements-report.txt后，Windows可执行`python tools/build_report.py`。其他平台通过`--body-font`和`--bold-font`指定包含中文的TrueType字体；排版有变更时需重新逐页检查。报告中的实验数值来源于归档工程，未根据重建结果另行调参。

本次交付前的隔离复现记录见[reproduction_check.json](verification/reproduction_check.json)，两阶段数值摘要和追加门槛诊断全部一致。
