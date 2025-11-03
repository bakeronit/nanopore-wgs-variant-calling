#!/usr/bin/env python3
import argparse
import json
import re
import subprocess
import sys
from datetime import timedelta
from pathlib import Path
from typing import List, Optional

# eg. "whatshap_phasing" from "snakejob.whatshap_phasing.999.sh"
JOB_NAME_REGEX = re.compile(r"snakejob\.(\w+)\.\d+\.sh")
QSTAT_COMMAND = "qstat -fx -F json"
PBS_SERVER = "hpcpbs02"

def format_memory(mem: str) -> float:
    """Convert memory from KB to GB."""
    mem = int(mem.replace("kb", ""))
    return round(mem / 1024 / 1024, 2)

def parse_time(time_str: str) -> timedelta:
    """Parse time string in HH:MM:SS format."""
    h, m, s = map(int, time_str.split(":"))
    return timedelta(hours=h, minutes=m, seconds=s)

def calculate_average_cpu(cpu_time: str, walltime_used: str) -> float:
    """Calculate average CPU utilization."""
    cpu_td = parse_time(cpu_time)
    wall_td = parse_time(walltime_used)
    return round(cpu_td / wall_td, 2)

def get_wildcards(job_id: str, job_name: str, rundir: Path) -> Optional[str]:
    """Extract Snakemake wildcards from PBS error file."""
    job_id = job_id.split(".")[0]
    error_files = list(rundir.glob(f"snakejob.{job_name}.*.sh.e{job_id}"))
    if not error_files:
        return None

    try:
        with open(error_files[0], 'r') as f:
            for line in f:
                if "wildcards:" in line:
                    return line.strip().split(":", 1)[1].strip()
    except (IOError, IndexError):
        return None
    return None

def parse_job_info(job_info: dict) -> dict:
    """Parse raw PBS job info into structured dictionary."""
    job_name = JOB_NAME_REGEX.search(job_info["Job_Name"]).group(1)
    exec_host = job_info['exec_host'].split("/")[0]
    return {
        "job_name": job_name,
        "exec_host": exec_host,
        "cpu_request": int(job_info['Resource_List']['ncpus']),
        "cpu_max": job_info['resources_used']['cpupercent'] / 100,
        "cpu_average": calculate_average_cpu(job_info['resources_used']['cput'], job_info['resources_used']['walltime']),
        "mem_request": int(job_info['Resource_List']['mem'].replace("gb", "")),
        "mem_used": format_memory(job_info['resources_used']['mem']),
        "walltime_request": job_info['Resource_List']['walltime'],
        "walltime_used": job_info['resources_used']['walltime']
    }

def get_jobs_info(job_ids: List[str], rundir: Path) -> List[dict]:
    """Query qstat and retrieve job information."""
    jobs_data = []
    cmd = f"{QSTAT_COMMAND} {' '.join(job_ids)}"
    result = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
    qstat_json = json.loads(result.stdout)
    for job_id in job_ids:
        full_job_id = f"{job_id}.{PBS_SERVER}"
        job_data = parse_job_info(qstat_json["Jobs"][full_job_id])
        job_data["job_id"] = full_job_id
        job_data["wildcards"] = get_wildcards(full_job_id, job_data["job_name"], rundir)
        jobs_data.append(job_data)
    return jobs_data

def print_table_simple(jobs_data: List[dict]) -> None:
    """Print job data as a simple formatted table."""
    if not jobs_data:
        print("No job data to display")
        return
    headers = {
        "job_id": ("Job ID", 25),
        "wildcards": ("Wildcards", 25),
        "job_name": ("Job Name", 20),
        "exec_host": ("Host", 15),
        "cpu_request": ("CPU Req", 8),
        "cpu_max": ("CPU Max", 8),
        "cpu_average": ("CPU Avg", 8),
        "mem_request": ("Mem Req(GB)", 11),
        "mem_used": ("Mem Used(GB)", 12),
        "walltime_request": ("Time Req", 12),
        "walltime_used": ("Time Used", 12)
    }
    
    header_line = ""
    separator_line = ""
    for key, (header, width) in headers.items():
        header_line += f"{header:<{width}} "
        separator_line += "-" * width + " "
    
    print(header_line)
    print(separator_line)
    
    # Print data rows
    for job in jobs_data:
        row = ""
        for key, (_, width) in headers.items():
            value = str(job.get(key, ""))
            if len(value) > width:
                value = value[:width-3] + "..."
            row += f"{value:<{width}} "
        print(row)

def print_table_tabulate(jobs_data: List[dict]) -> None:
    """Print job data using tabulate library."""
    try:
        from tabulate import tabulate
        
        if not jobs_data:
            print("No job data to display")
            return
        headers = [
            "job_id", "job_name", "wildcards","exec_host", 
            "cpu_request", "cpu_max", "cpu_average",
            "mem_request", "mem_used", 
            "walltime_request", "walltime_used"
        ]
        table_data = []
        for job in jobs_data:
            row = [job.get(h, "") for h in headers]
            table_data.append(row)
        pretty_headers = [
            "Job ID", "Job Name", "Wildcards", "Host", 
            "CPU Req", "CPU Max", "CPU Avg",
            "Mem Req(GB)", "Mem Used(GB)", 
            "Time Req", "Time Used"
        ]
        print(tabulate(table_data, headers=pretty_headers, tablefmt="grid"))
        
    except ImportError:
        print("tabulate library not found. Using simple table format instead.")
        print_table_simple(jobs_data)

def print_summary_stats(jobs_data: List[dict]) -> None:
    """Print summary statistics of the jobs."""
    if not jobs_data:
        return
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    total_jobs = len(jobs_data)
    avg_cpu_usage = sum(job['cpu_average'] for job in jobs_data) / total_jobs
    max_cpu_usage = max(job['cpu_max'] for job in jobs_data)
    avg_mem_usage = sum(job['mem_used'] for job in jobs_data) / total_jobs
    max_mem_usage = max(job['mem_used'] for job in jobs_data)
    print(f"Total Jobs: {total_jobs}")
    print(f"Average CPU Usage: {avg_cpu_usage:.2f}")
    print(f"Max CPU Usage: {max_cpu_usage:.2f}")
    print(f"Average Memory Usage: {avg_mem_usage:.2f} GB")
    print(f"Max Memory Usage: {max_mem_usage:.2f} GB")
    job_groups = {}
    for job in jobs_data:
        name = job['job_name']
        if name not in job_groups:
            job_groups[name] = []
        job_groups[name].append(job)
    print(f"\nJobs by Type:")
    for name, jobs in sorted(job_groups.items()):
        print(f"  {name}: {len(jobs)} jobs")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Analyze PBS job statistics from Snakemake workflows"
    )
    parser.add_argument(
        'rundir',
        type=Path,
        help='Directory containing Snakemake PBS job files'
    )
    parser.add_argument(
        '--summary', '-s',
        action='store_true',
        help='Show summary statistics'
    )
    parser.add_argument(
        '--sort',
        choices=['job_name', 'cpu_average', 'mem_used'],
        default='job_name',
        help='Sort jobs by field (default: job_name)'
    )
    args = parser.parse_args()
    if not args.rundir.exists():
        print(f"Error: Directory not found: {args.rundir}", file=sys.stderr)
        sys.exit(1)

    # find job ids from error files from snakemake pipeline
    job_ids = [fn.suffix[2:] for fn in args.rundir.glob("snakejob.*.sh.e*")]
    if not job_ids:
        print(f"Error: No job files found in {args.rundir}", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(job_ids)} jobs in {args.rundir}")
    try:
        jobs_data = get_jobs_info(job_ids, args.rundir)
    except Exception as e:
        print(f"Error querying jobs: {e}", file=sys.stderr)
        sys.exit(1)
    jobs_data.sort(key=lambda x: x.get(args.sort, ''))
    print_table_tabulate(jobs_data)
    if args.summary:
        print_summary_stats(jobs_data)

if __name__ == "__main__":
    main()