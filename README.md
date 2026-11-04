# Svanberg, Cui & Latter (2022) Figure 2：Athena++ 模拟

论文：*Wavelike nature of the vertical shear instability in global protoplanetary discs*, MNRAS 514, 4581–4587，DOI 10.1093/mnras/stac1598。

本项目实现 Section 3 的轴对称二维纯流体盘，目标是 Figure 2 的 **t = 500 P0** 快照。当前完整计算在 `results/figure2/`；只有出现 `COMPLETE` 后，才表示目标模拟和绘图均已成功结束。论文中 3000 P0 用于后续波分析，本图只需 500 P0。

## 物理与数值设置

- 球坐标 (r, theta, phi)，轴对称（nx3=1），无磁场、无显式黏性。
- GM=R0=rho0=1，P0=2π，gamma=7/5，h=0.05。
- cs²=0.0025/R，rho_mid=R^-2，三维垂直平衡密度与旋转速度严格采用式 (24)、(25)。内置点质量源项负责引力，避免重复施加。
- r∈[1,30]、theta∈[π/2−0.25,π/2+0.25]，704×448 网格；径向相邻单元宽度比 x1rat=1.005。
- 四个活动边界均在幽灵单元恢复平衡解；不另加缓冲层。
- 每个速度分量独立加入 [-0.01 cs,0.01 cs] 均匀噪声，固定种子 202211；噪声由坐标散列产生，与分块无关。
- 每个完整步之后，对内能应用式 (27) 的精确指数冷却，tau=10^-10 P0。保持动能不变，同时更新守恒能量与原始压力；没有以显式 dt/tau 限制时间步。
- 采用绝热能量方程、VL2、PLM（二阶）、HLLC、CFL=0.3，双精度，4 个 OpenMP 线程。论文未指定这些全部算法细节及种子，属于明确记录的实现选择，不保证与原图逐像素一致。
- 注意论文文字内部存在分辨率估算不一致：10H 的 theta 区间用 448 格意味着约 44.8 格/H，并非文中约 12 格/H；这里优先采用明确写出的网格数与范围。另外，704 格、[1,30] 与严格 Δr∝r 的宽度比应为 30^(1/704)≈1.004843；这里按文中的相邻宽度比 1.005，由 Athena++ 保持端点，所得网格为几何递增宽度。

## 构建与运行

需要 g++、make、Python 3、NumPy、Matplotlib。Athena++ 上游版本固定为 `053bc738d5c85dac717326b4338b34c7db4e715d`。构建脚本把独立问题生成器复制到上游，不修改上游 disk.cpp 或求解器。

```bash
cd /workspace/course_VSI
bash scripts/build.sh
bash scripts/run_figure2.sh
```

长任务后台运行（本次已启动，勿重复启动）：

```bash
nohup bash scripts/run_figure2.sh > results/figure2/run.log 2>&1 < /dev/null &
```

runner 用锁防止同目录重复运行；若中断，重新执行同一命令会从 final 或最近的编号 restart 文件续算。运行中的进程不保证随云环境休眠或重启保留。使用现有 checkout；除非用户明确要求，不创建 Git worktree。

输入文件 `inputs/athinput.figure2` 默认运行至 500P0。每 1P0 输出完整 VTK（当前从约 10P0 的检查点切换，之前时段没有每轨道快照），每 10P0 保存重启，每 0.1P0 输出 history。达到目标后自动生成：

- `results/figure2/figure2.png`：双面板图；左侧带符号对数色标表示 vz/cs，右侧为 (vR,vz)/cs 箭头。
- `results/figure2/diagnostics.json`：目标时间、网格、冷却偏差及内盘 vz/cs RMS。
- `results/figure2/COMPLETE`：模拟及绘图成功标志。

查看进展：

```bash
tail -n 5 results/figure2/run.log
```

日志 time 为代码单位；除以 2π 得 P0。当前机器约 4 CPU 配额，200 步实测耗时 8.52 秒，到 500P0 线性外推约 6.7 小时；饱和后的速度与负载会改变实际时间。

手工重画（仅在目标快照已产生后）：

```bash
python3 scripts/plot_figure2.py results/figure2 --require-target --output results/figure2/figure2.png
```

坐标使用球坐标体积平均中心，vz=vr cos(theta)−vtheta sin(theta)，vR=vr sin(theta)+vtheta cos(theta)。读取及检查所有径向分块；VTK 存储为单精度，不改变求解器双精度计算。

## 已执行的验证

完整 704×448 网格运行了 200 步，无非有限量、负密度或负压；初始平衡密度相对误差约 6×10^-8（输出舍入），冷却温度相对误差约 1.2×10^-7。另对比连续 20 步与 10 步后重启再运行 10 步，输出在单精度容差内一致（最大密度差约 1.5×10^-8）。

```bash
python3 scripts/validate.py
```

验证报告见 `results/validation.json`。`results/benchmark/early_snapshot.png` 仅为 t≈0.16P0 的绘图检查，不是 Figure 2 的复现结果。500P0 的波结构及幅度仍须等完整模拟结束后评估。

复现短程验证数据的命令：

```bash
bin/athena -i inputs/athinput.figure2 -d results/benchmark time/nlim=200 output1/dt=0.5 output2/dt=0.5
bin/athena -i inputs/athinput.figure2 -d results/restart_split time/nlim=10 output1/dt=0.01 output2/dt=0.01
bin/athena -r results/restart_split/vsi.final.rst -d results/restart_split time/nlim=20
bin/athena -i inputs/athinput.figure2 -d results/restart_direct time/nlim=20 output1/dt=0.01 output2/dt=0.01
```

运行这些验证前先创建对应输出目录，保留已有结果，勿同时向相同目录写入。

## Figure 2 样式的连续 MP4

`scripts/run_movie_data.sh` 用相同分辨率、平衡盘和确定性扰动单独补算缺失的 0–11 P0 快照，不改动主模拟数据。常规输出在首次跨过整数轨道时刻的时间步写出，误差小于一个时间步；图中标出实际时间。

```bash
python3 scripts/make_movie.py
```

该命令合并 `results/movie_early/` 和 `results/figure2/` 的快照，严格按 0,1,2,… 选取连续序列，遇到缺帧即停止，不跳过时刻。前期补算与主模拟使用相同初始条件及数值方法，主模拟在约 10P0 有一次重启，两段不是完全相同的逐步浮点轨迹。视频不是时间插值，也不是把论文图像变形得到。

输出为 `results/movie_figure2/figure2_evolution.mp4`；每个快照持续 0.25 秒，H.264、1440×960、24 fps；采用深色背景、青蓝—正红带符号对数色标，所有帧保持同一色标与坐标范围。仅保留流场主面板，标题显示整数轨道时刻，实际输出时间仍记录在元数据中。`movie_metadata.json` 记录每帧实际模拟时间与数据来源。PNG 与单帧 MP4 会缓存，重运行时只渲染新增快照并更新视频。

`python3 scripts/watch_movie.py` 可持续监测新数据并自动更新视频，直到最后一个 500P0 快照合成成功。当前监测进程已经在云环境中启动，进程不会保证随环境休眠保留。

仓库中的 [`media/figure2_evolution.mp4`](media/figure2_evolution.mp4) 是已上传的 **0–25 P0** 视频（26 个时刻，每帧 0.25 秒）；对应真实时间记录见 `media/movie_metadata.json`。云端 `results/movie_figure2/figure2_evolution.mp4` 会随新数据继续更新；GitHub 中的版本需要再次提交才会更新。500P0 的计算尚未完成。

图像和视频的标题、刻度及公式统一使用 Times 风格；当前云环境采用 Times 兼容的 Nimbus Roman。若系统安装了 Times New Roman，字体配置会优先使用它。
