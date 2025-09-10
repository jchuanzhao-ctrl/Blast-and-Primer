#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CRISPRessoPooled Pipeline with GUI improvements.

This script wraps a CRISPRessoPooled analysis workflow and provides a
Tkinter-based interface.  The original logic is preserved while adding
quality-of-life improvements such as better error handling and GUI
progress feedback.
"""

import sys
import shutil
import os
import re
import gzip
import queue
import subprocess
from datetime import datetime
from threading import Thread

import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# ------------------ 常用工具 ------------------
COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def rc(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    return (seq or "").translate(COMPLEMENT)[::-1]


def npath(p: str) -> str:
    """Normalize a filesystem path, returning the input if falsy."""
    return os.path.normpath(p) if p else p


def ensure_dir(d: str):
    """Create directory *d* if it does not already exist."""
    os.makedirs(d, exist_ok=True)


def file_nonempty(p: str) -> bool:
    """Return True if file *p* exists and is non-empty."""
    return os.path.exists(p) and os.path.getsize(p) > 0


def is_dna(s: str) -> bool:
    """Check whether a string consists solely of A/T/C/G characters."""
    if not isinstance(s, str) or not s.strip():
        return False
    s = re.sub(r"\s+", "", s.strip().upper())
    return all(ch in "ATCG" for ch in s)


def to_docker_subpath(host_path: str, mounted_host_root: str, container_root="/DATA") -> str:
    """把宿主机 host_path（位于 mounted_host_root 下）转成容器内路径 /DATA/rel"""
    rel = os.path.relpath(os.path.abspath(host_path), os.path.abspath(mounted_host_root))
    rel = rel.replace("\\", "/")
    return f"{container_root}/{rel}"


# ------------------ 线程安全日志 ------------------
class TSLogger:
    """Thread-safe logger that writes to both a Tk Text widget and a file."""

    def __init__(self, widget: tk.Text, log_path: str):
        self.q = queue.Queue()
        self.widget = widget
        self.log_path = log_path
        ensure_dir(os.path.dirname(log_path))
        with open(log_path, "w", encoding="utf-8") as f:
            f.write(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] 日志开始\n")

    def put(self, msg: str):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        self.q.put(line)
        try:
            with open(self.log_path, "a", encoding="utf-8") as f:
                f.write(line + "\n")
        except Exception:
            pass

    def pump(self):
        try:
            while True:
                s = self.q.get_nowait()
                self.widget.insert(tk.END, s + "\n")
                self.widget.see(tk.END)
        except queue.Empty:
            pass
        finally:
            self.widget.after(60, self.pump)


# ------------------ 外部程序检测 ------------------
def run_cmd(cmd, logger: TSLogger, timeout=None, check=True, cwd=None):
    """Run an external command and log its output."""
    logger.put(f"[CMD] {' '.join(cmd)}")
    try:
        p = subprocess.run(cmd, text=True, capture_output=True, timeout=timeout, cwd=cwd)
    except FileNotFoundError:
        logger.put(f"[ERROR] command not found: {cmd[0]}")
        raise
    if p.stdout:
        for line in p.stdout.splitlines():
            logger.put(line)
    if p.stderr:
        for line in p.stderr.splitlines():
            logger.put(line)
    if check and p.returncode != 0:
        raise RuntimeError(f"命令失败，退出码 {p.returncode}")
    return p


def which(exe):
    from shutil import which as _which
    return _which(exe) is not None


def check_cutadapt(logger: TSLogger):
    py = sys.executable.replace("pythonw.exe", "python.exe")
    run_cmd([py, "-m", "cutadapt", "--version"], logger, timeout=120, check=True)


def check_docker(logger: TSLogger):
    try:
        run_cmd(["docker", "version"], logger, timeout=60, check=True)
    except FileNotFoundError:
        raise RuntimeError("未找到 docker 命令，请安装 Docker 并确保其在 PATH 中")
    except Exception as e:
        raise RuntimeError("Docker 未运行或不可用，请确保 Docker 服务已启动") from e


# ------------------ FASTQ/Excel 检查 ------------------
def check_fastq_gz(p, logger: TSLogger):
    if not file_nonempty(p):
        raise RuntimeError(f"FASTQ 文件不存在或为空：{p}")
    try:
        with gzip.open(p, "rt") as fh:
            head = fh.readline(1 << 15)
            if not head.startswith("@"):
                raise RuntimeError(f"不是有效的 FASTQ：{p}")
    except Exception as e:
        raise RuntimeError(f"读取 FASTQ 失败：{p} -> {e}")
    logger.put(f"[OK] FASTQ 检查通过：{p}")


def read_design(xlsx_path, logger: TSLogger) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path)
    df.columns = [str(c).strip() for c in df.columns]
    need = ["样品名", "扩增子", "spacer", "正向PCR引物(5'-3')", "反向PCR引物(5'-3')"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise RuntimeError(f"Excel 缺少必要列：{missing}")
    # 校验
    for i, r in df.iterrows():
        name = str(r["样品名"]).strip()
        amp = str(r["扩增子"]).strip().upper()
        sp = str(r["spacer"]).strip().upper()
        fwd = str(r["正向PCR引物(5'-3')"]).strip().upper()
        rev = str(r["反向PCR引物(5'-3')"]).strip().upper()
        if not name:
            raise RuntimeError(f"第 {i+2} 行样品名为空")
        if not (is_dna(amp) and is_dna(sp) and is_dna(fwd) and is_dna(rev)):
            raise RuntimeError(f"第 {i+2} 行存在非 DNA 字符：{name}")
        if sp not in amp:
            logger.put(f"[WARN] spacer 不在扩增子内：{name} ({sp[:12]}...)")
    logger.put(f"[OK] 读取 Excel：{xlsx_path}；行数={len(df)}")
    return df


# ------------------ amplicon 内切（去掉两端引物） ------------------
def inner_amplicon(amp_seq, fwd, rev, logger: TSLogger, sample: str) -> str:
    amp = re.sub(r"\s+", "", str(amp_seq or "").upper())
    F = re.sub(r"\s+", "", str(fwd or "").upper())
    Rrc = rc(str(rev or "").upper())
    i = amp.find(F)
    j = amp.rfind(Rrc)
    if i != -1 and j != -1 and j >= i + len(F):
        inner = amp[i + len(F):j]
        if len(inner) < 20:
            logger.put(f"[WARN] {sample} 内切后过短({len(inner)}bp)，改用原扩增子")
            return amp
        logger.put(f"[OK] {sample} amplicon 已内切：{len(amp)}bp -> {len(inner)}bp")
        return inner
    else:
        logger.put(f"[INFO] {sample} 未检测到 F/rc(R) 同向包裹，使用原扩增子")
        return amp


# ------------------ Demux：按 F/R 引物匹配并去除 5' 引物 ------------------
def demux_by_primers(df, r1, r2, out_dir, threads, logger: TSLogger):
    tmp = os.path.join(out_dir, ".__tmp"); ensure_dir(tmp)
    demux_dir = os.path.join(tmp, "demux"); ensure_dir(demux_dir)
    ffa = os.path.join(tmp, "Fprimers.fasta")
    rfa = os.path.join(tmp, "Rprimers.fasta")
    with open(ffa, "w", encoding="ascii") as fF, open(rfa, "w", encoding="ascii") as fR:
        for _, row in df.iterrows():
            name = str(row["样品名"]).strip()
            fwd = str(row["正向PCR引物(5'-3')"]).strip().upper()
            rev = str(row["反向PCR引物(5'-3')"]).strip().upper()
            fF.write(f">{name}\n{fwd}\n")
            fR.write(f">{name}\n{rev}\n")
    py = sys.executable.replace("pythonw.exe", "python.exe")
    cmd = [
        py, "-m", "cutadapt",
        "--pair-adapters",
        "-g", f"file:{npath(ffa)}",
        "-G", f"file:{npath(rfa)}",
        "--discard-untrimmed",
        "--minimum-length", "30",
        "--maximum-length", "1000",
        "-e", "0.1",
        "-j", str(threads),
        "-o", npath(os.path.join(demux_dir, "{name}.R1.fastq.gz")),
        "-p", npath(os.path.join(demux_dir, "{name}.R2.fastq.gz")),
        npath(r1), npath(r2)
    ]
    run_cmd(cmd, logger, timeout=6 * 3600, check=True)
    produced = {}
    for f in os.listdir(demux_dir):
        if f.endswith(".R1.fastq.gz"):
            name = f[:-len(".R1.fastq.gz")]
            r1p = os.path.join(demux_dir, f)
            r2p = os.path.join(demux_dir, f"{name}.R2.fastq.gz")
            if file_nonempty(r1p) and file_nonempty(r2p):
                produced[name] = (r1p, r2p)
    logger.put(f"[OK] demux 完成；样本数={len(produced)}")
    return demux_dir, produced


# ------------------ 样本级 3' 去互补引物（不丢 reads） ------------------
def trim_3prime_per_sample(df, demux_dir, out_dir, threads, logger: TSLogger):
    trimmed_dir = os.path.join(out_dir, ".__tmp", "trimmed"); ensure_dir(trimmed_dir)
    name2row = {str(r["样品名"]).strip(): r for _, r in df.iterrows()}
    kept = {}
    py = sys.executable.replace("pythonw.exe", "python.exe")
    for f in os.listdir(demux_dir):
        if not f.endswith(".R1.fastq.gz"):
            continue
        name = f[:-len(".R1.fastq.gz")]
        r1p = os.path.join(demux_dir, f)
        r2p = os.path.join(demux_dir, f"{name}.R2.fastq.gz")
        if not (file_nonempty(r1p) and file_nonempty(r2p)):
            continue
        row = name2row.get(name)
        if row is None:
            logger.put(f"[WARN] 缺少 {name} 的引物定义，跳过 3' 修剪")
            continue
        F = str(row["正向PCR引物(5'-3')"]).strip().upper()
        R = str(row["反向PCR引物(5'-3')"]).strip().upper()
        # 仅去 3' 端互补引物；不 discard 未匹配，避免损失覆盖不足的 reads
        cmd = [
            py, "-m", "cutadapt",
            "-j", str(threads),
            "--minimum-length", "30",
            "-a", rc(R),  # R1 末端去 rc(R)
            "-A", rc(F),  # R2 末端去 rc(F)
            "-o", os.path.join(trimmed_dir, f"{name}.R1.fastq.gz"),
            "-p", os.path.join(trimmed_dir, f"{name}.R2.fastq.gz"),
            r1p, r2p
        ]
        run_cmd(cmd, logger, timeout=3600, check=True)
        out1 = os.path.join(trimmed_dir, f"{name}.R1.fastq.gz")
        out2 = os.path.join(trimmed_dir, f"{name}.R2.fastq.gz")
        if file_nonempty(out1) and file_nonempty(out2):
            kept[name] = (out1, out2)
    logger.put(f"[OK] 3' 引物精修完成；样本数={len(kept)}")
    return trimmed_dir, kept


# ------------------ 合并若干 fastq.gz ------------------
def cat_files(files, out_path, logger: TSLogger):
    ensure_dir(os.path.dirname(out_path))
    with open(out_path, "wb") as w:
        for f in files:
            with open(f, "rb") as r:
                shutil.copyfileobj(r, w)
    logger.put(f"[OK] 合并 -> {out_path}")


# ------------------ 写 pooled tsv（自动内切扩增子） ------------------
def write_pooled_tsv(df, out_dir, logger: TSLogger) -> str:
    tmp = os.path.join(out_dir, ".__tmp"); ensure_dir(tmp)
    out_path = os.path.join(tmp, "pooled_amplicons.tsv")
    with open(out_path, "w", encoding="ascii") as f:
        f.write("amplicon_nameamplicon_seqguide_seq\n")
        for _, r in df.iterrows():
            name = str(r["样品名"]).strip()
            amp = str(r["扩增子"]).strip().upper()
            sp = str(r["spacer"]).strip().upper()
            F = str(r["正向PCR引物(5'-3')"]).strip().upper()
            R = str(r["反向PCR引物(5'-3')"]).strip().upper()
            amp_inner = inner_amplicon(amp, F, R, logger, name)
            if sp:
                if sp not in amp_inner:
                    rc_amp = rc(amp_inner)
                    if sp in rc_amp:
                        logger.put(f"[OK] {name} spacer 不在扩增子序列中，已改用互补链序列")
                        amp_inner = rc_amp
                    else:
                        logger.put(f"[WARN] {name} spacer 在扩增子序列及其互补链中均未找到，仍使用提供的 spacer")
            f.write(f"{name}\t{amp_inner}\t{sp}\n")
    logger.put(f"[OK] 写 TSV -> {out_path}")
    return out_path


# ------------------ 运行 CRISPRessoPooled（Docker） ------------------
def run_crispresso_pooled(pooled_tsv, combined_r1, combined_r2, out_dir, name, q_center, q_size, threads, mode, logger: TSLogger):
    ensure_dir(out_dir)
    # 宿主机 -> 容器挂载
    out_dir_abs = os.path.abspath(out_dir)
    input_dir = os.path.abspath(os.path.dirname(combined_r1))  # /INPUT 指向 combined 文件所在目录
    mounts = [
        "-v", f"{npath(out_dir_abs)}:/DATA",
        "-v", f"{npath(input_dir)}:/INPUT",
        "-w", "/DATA",
    ]
    # 给 fastp 报告准备输出目录（容器里写入）
    ensure_dir(os.path.join(out_dir, f"CRISPRessoPooled_on_{name}"))

    fastp_opts = (
        f"--merge "
        f"--merged_out /DATA/CRISPRessoPooled_on_{name}/out.extendedFrags.fastq.gz "
        f"--unpaired1 /DATA/CRISPRessoPooled_on_{name}/out.notCombined_1.fastq.gz "
        f"--unpaired2 /DATA/CRISPRessoPooled_on_{name}/out.notCombined_2.fastq.gz "
        f"--overlap_len_require 20 --correction "
        f"--thread {threads} "
        f"--json /DATA/CRISPRessoPooled_on_{name}/fastp_report.json "
        f"--html /DATA/CRISPRessoPooled_on_{name}/fastp_report.html"
    )

    # 把 TSV 转成容器内可见路径
    f_arg = f"/INPUT/{os.path.basename(pooled_tsv)}"
    r1_arg = f"/INPUT/{os.path.basename(combined_r1)}"
    r2_arg = f"/INPUT/{os.path.basename(combined_r2)}"

    cmd = ["docker", "run", "--rm", *mounts,
           "pinellolab/crispresso2", "CRISPRessoPooled",
           "-r1", r1_arg,
           "-r2", r2_arg,
           "-f", f_arg,
           "-o", "/DATA",
           "--name", name,
           "--quantification_window_center", str(q_center),
           "--quantification_window_size", str(q_size),
           "--plot_window_size", "40",
           "--min_paired_end_reads_overlap", "20",
           "--fastp_options_string", fastp_opts,
           "-p", str(threads)]
    if str(mode).upper() == "BE":
        cmd.append("--base_editor_output")

    try:
        run_cmd(cmd, logger, timeout=6 * 3600, check=True)
    except FileNotFoundError:
        raise RuntimeError("未找到 docker 命令，请确认已安装并配置 Docker")
    except RuntimeError as e:
        raise RuntimeError("Docker 运行失败，请检查 Docker 服务是否启动") from e


# ------------------ 摘要 & 清理 ------------------
def write_summary_html(out_dir, name, logger: TSLogger):
    html = os.path.join(out_dir, "run_summary.html")
    folder = os.path.join(out_dir, f"CRISPRessoPooled_on_{name}")
    ok = os.path.isdir(folder)
    with open(html, "w", encoding="utf-8") as f:
        f.write("<!doctype html><meta charset='utf-8'>")
        f.write("<title>CRISPResso 运行摘要</title>")
        f.write("<h1>CRISPResso 运行摘要</h1>")
        f.write(f"<p>输出目录：{folder}</p>")
        f.write(f"<p>状态：{'成功' if ok else '失败/未找到'}</p>")
        if ok:
            f.write("<ul>")
            for fn in sorted(os.listdir(folder)):
                f.write(f"<li>{fn}</li>")
            f.write("</ul>")
    logger.put(f"[OK] 写出摘要：{html}")


def cleanup_temp(out_dir, name, logger: TSLogger):
    keep_files = {"run.log", "run_summary.html"}
    keep_dir = os.path.join(out_dir, f"CRISPRessoPooled_on_{name}")
    for root, dirs, files in os.walk(out_dir, topdown=False):
        for fn in files:
            p = os.path.join(root, fn)
            if p.startswith(keep_dir):
                continue
            if os.path.basename(p) in keep_files:
                continue
            if "__tmp" in p or os.path.basename(p) in {"combined_R1.fastq.gz", "combined_R2.fastq.gz", "pooled_amplicons.tsv"}:
                try:
                    os.remove(p)
                except Exception:
                    pass
        for d in dirs:
            p = os.path.join(root, d)
            if "__tmp" in p:
                shutil.rmtree(p, ignore_errors=True)
    logger.put("[OK] 已清理临时文件，仅保留 CRISPResso 输出与日志/摘要")


# ------------------ 主流程 ------------------
def run_pipeline(xlsx, r1, r2, out_dir, name, q_center, q_size, threads, mode, logger: TSLogger):
    logger.put("开始运行...")
    check_cutadapt(logger)
    check_docker(logger)
    check_fastq_gz(r1, logger); check_fastq_gz(r2, logger)
    df = read_design(xlsx, logger)

    # 1) demux：去 5' 引物，只保留匹配样本
    demux_dir, produced = demux_by_primers(df, r1, r2, out_dir, threads, logger)

    if not produced:
        raise RuntimeError("demux 未得到任何样本")

    # 2) 样本级 3' 去互补引物（不丢 reads）
    trimmed_dir, kept = trim_3prime_per_sample(df, demux_dir, out_dir, threads, logger)
    if not kept:
        raise RuntimeError("3' 精修后无样本")

    # 3) 合并所有样本到 combined_R1/combined_R2（供 pooled 分析）
    tmp = os.path.join(out_dir, ".__tmp"); ensure_dir(tmp)
    combined_r1 = os.path.join(tmp, "combined_R1.fastq.gz")
    combined_r2 = os.path.join(tmp, "combined_R2.fastq.gz")
    cat_files([p[0] for p in kept.values()], combined_r1, logger)
    cat_files([p[1] for p in kept.values()], combined_r2, logger)

    # 4) 生成 pooled tsv（自动内切扩增子）
    pooled_tsv = write_pooled_tsv(df, out_dir, logger)

    # 5) Docker 运行 CRISPRessoPooled（由其内置 fastp 合并）
    run_crispresso_pooled(pooled_tsv, combined_r1, combined_r2, out_dir, name, q_center, q_size, threads, mode, logger)

    # 6) 摘要与清理
    write_summary_html(out_dir, name, logger)
    cleanup_temp(out_dir, name, logger)
    logger.put("✅ 成功完成！仅保留 CRISPResso 输出 + summary + run.log。")


# ------------------ GUI ------------------
class App(ttk.Frame):
    """Tkinter front-end for the pipeline."""

    def __init__(self, master):
        super().__init__(master, padding=12)
        master.title("CRISPRessoPooled Pipeline")
        master.geometry("980x700"); master.minsize(900, 620)

        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure("TButton", padding=(10, 6))
        style.configure("TEntry", padding=4)
        style.configure("TLabel", padding=4)

        frm = ttk.LabelFrame(self, text="输入与参数", padding=12)
        frm.grid(row=0, column=0, sticky="nsew"); frm.columnconfigure(1, weight=1)

        self.e_xlsx = ttk.Entry(frm); self.e_r1 = ttk.Entry(frm); self.e_r2 = ttk.Entry(frm)
        self.e_out = ttk.Entry(frm); self.e_name = ttk.Entry(frm)
        self.e_name.insert(0, "combined_R1_combined_R2")
        self.e_qc = ttk.Entry(frm, width=10); self.e_qc.insert(0, "0")
        self.e_qs = ttk.Entry(frm, width=10); self.e_qs.insert(0, "21")
        self.e_thr = ttk.Entry(frm, width=10); self.e_thr.insert(0, "8")
        self.mode = ttk.Combobox(frm, values=["NHEJ", "BE"], width=10); self.mode.current(0)

        row = 0
        ttk.Label(frm, text="Excel 文件：").grid(row=row, column=0, sticky="e"); self.e_xlsx.grid(row=row, column=1, sticky="ew")
        ttk.Button(frm, text="选择", command=self.sel_xlsx).grid(row=row, column=2, padx=6); row += 1
        ttk.Label(frm, text="R1 FASTQ：").grid(row=row, column=0, sticky="e"); self.e_r1.grid(row=row, column=1, sticky="ew")
        ttk.Button(frm, text="选择", command=self.sel_r1).grid(row=row, column=2, padx=6); row += 1
        ttk.Label(frm, text="R2 FASTQ：").grid(row=row, column=0, sticky="e"); self.e_r2.grid(row=row, column=1, sticky="ew")
        ttk.Button(frm, text="选择", command=self.sel_r2).grid(row=row, column=2, padx=6); row += 1
        ttk.Label(frm, text="输出目录：").grid(row=row, column=0, sticky="e"); self.e_out.grid(row=row, column=1, sticky="ew")
        ttk.Button(frm, text="选择", command=self.sel_out).grid(row=row, column=2, padx=6); row += 1
        ttk.Label(frm, text="任务名称：").grid(row=row, column=0, sticky="e"); self.e_name.grid(row=row, column=1, sticky="w"); row += 1
        ttk.Label(frm, text="Quant center：").grid(row=row, column=0, sticky="e"); self.e_qc.grid(row=row, column=1, sticky="w"); row += 1
        ttk.Label(frm, text="Quant size：").grid(row=row, column=0, sticky="e"); self.e_qs.grid(row=row, column=1, sticky="w"); row += 1
        ttk.Label(frm, text="线程数：").grid(row=row, column=0, sticky="e"); self.e_thr.grid(row=row, column=1, sticky="w"); row += 1
        ttk.Label(frm, text="模式：").grid(row=row, column=0, sticky="e"); self.mode.grid(row=row, column=1, sticky="w")

        ctrl = ttk.Frame(self, padding=(0, 12, 0, 6))
        ctrl.grid(row=1, column=0, sticky="w")
        self.start_btn = ttk.Button(ctrl, text="运行 Pipeline", command=self.start)
        self.start_btn.grid(row=0, column=0, padx=(0, 10))
        ttk.Button(ctrl, text="退出", command=self.master.destroy).grid(row=0, column=1)
        self.pb = ttk.Progressbar(ctrl, mode="indeterminate")
        self.pb.grid(row=1, column=0, columnspan=2, pady=(8, 0), sticky="ew")
        self.pb.grid_remove()

        logf = ttk.LabelFrame(self, text="运行日志", padding=6)
        logf.grid(row=2, column=0, sticky="nsew")
        self.text = tk.Text(logf, height=22, font=("Consolas", 10))
        self.text.grid(row=0, column=0, sticky="nsew")
        sb = ttk.Scrollbar(logf, command=self.text.yview); sb.grid(row=0, column=1, sticky="ns")
        self.text.configure(yscrollcommand=sb.set)
        logf.rowconfigure(0, weight=1); logf.columnconfigure(0, weight=1)

        self.status = ttk.Label(self, anchor="w")
        self.status.grid(row=3, column=0, sticky="ew")

        self.grid(row=0, column=0, sticky="nsew")
        self.master.rowconfigure(0, weight=1); self.master.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1); self.columnconfigure(0, weight=1)

        # 初始日志器
        out_dir = os.getcwd()
        self.logger = TSLogger(self.text, os.path.join(out_dir, "run.log"))
        self.logger.pump()

    def sel_xlsx(self):
        p = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx;*.xls")])
        if p: self.e_xlsx.delete(0, tk.END); self.e_xlsx.insert(0, p)

    def sel_r1(self):
        p = filedialog.askopenfilename(filetypes=[("FASTQ", "*.fq.gz;*.fastq.gz;*.fq;*.fastq")])
        if p: self.e_r1.delete(0, tk.END); self.e_r1.insert(0, p)

    def sel_r2(self):
        p = filedialog.askopenfilename(filetypes=[("FASTQ", "*.fq.gz;*.fastq.gz;*.fq;*.fastq")])
        if p: self.e_r2.delete(0, tk.END); self.e_r2.insert(0, p)

    def sel_out(self):
        p = filedialog.askdirectory()
        if p: self.e_out.delete(0, tk.END); self.e_out.insert(0, p)

    def start(self):
        try:
            xlsx = self.e_xlsx.get().strip()
            r1 = self.e_r1.get().strip()
            r2 = self.e_r2.get().strip()
            outd = self.e_out.get().strip() or os.getcwd()
            name = self.e_name.get().strip() or "combined_R1_combined_R2"
            qc = int(self.e_qc.get().strip())
            qs = int(self.e_qs.get().strip())
            thr = int(self.e_thr.get().strip())
            mode = self.mode.get().strip()
        except Exception:
            messagebox.showerror("错误", "参数格式不正确")
            return
        if not (xlsx and r1 and r2):
            messagebox.showerror("错误", "请提供 Excel、R1、R2")
            return
        for p, lbl in ((xlsx, "Excel"), (r1, "R1 FASTQ"), (r2, "R2 FASTQ")):
            if not os.path.exists(p):
                messagebox.showerror("错误", f"{lbl} 文件不存在")
                return
        ensure_dir(outd)
        # 将日志切换到输出目录
        self.logger = TSLogger(self.text, os.path.join(outd, "run.log"))
        self.logger.pump()
        self.status.config(text="运行中...", foreground="blue")
        self.start_btn.config(state="disabled")
        self.pb.grid()
        self.pb.start()

        def worker():
            err = None
            try:
                run_pipeline(xlsx, r1, r2, outd, name, qc, qs, thr, mode, self.logger)
            except Exception as e:
                err = e
                self.logger.put(f"[ERROR] {e}")
            finally:
                self.after(0, lambda: self.on_finish(err))

        Thread(target=worker, daemon=True).start()

    def on_finish(self, err=None):
        self.pb.stop()
        self.pb.grid_remove()
        self.start_btn.config(state="normal")
        if err is None:
            self.status.config(text="完成", foreground="green")
            messagebox.showinfo("完成", "Pipeline 成功完成")
        else:
            self.status.config(text="失败", foreground="red")
            messagebox.showerror("失败", str(err))


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
