from pathlib import Path
import json
import subprocess
import re 
from datetime import timedelta
from typing import Union

#get "whatshap_phasing" from "snakejob.whatshap_phasing.999.sh"
job_name_regex = re.compile(r"snakejob\.(\w+)\.\d+\.sh")
QSTAT_COMMAND = "qstat -fx -F json"
PBS_SERVER = "hpcpbs02"


def format_mem(mem: str) -> float:
    mem = int(mem.replace("kb", ""))
    return round(mem / 1024 / 1024, 2)

def parse_hms(hms: str) -> timedelta:
    h, m, s = map(int, hms.split(":"))
    return timedelta(hours=h, minutes=m, seconds=s)
def cal_average_cpu(cpu_time: str, walltime_used: str) -> float:
    cpu_time = parse_hms(cpu_time)
    walltime_used = parse_hms(walltime_used)
    return round(cpu_time / walltime_used, 2)

def parse_job_info(job_info: dict) -> dict:
    job_name = job_name_regex.search(job_info["Job_Name"]).group(1)
    exec_host = job_info['exec_host'].split("/")[0]
    return {
        "job_name": job_name,
        "exec_host": exec_host,
        "cpu_request": int(job_info['Resource_List']['ncpus']), 
        "cpu_max": job_info['resources_used']['cpupercent'] / 100,
        "cpu_average": cal_average_cpu(job_info['resources_used']['cput'], job_info['resources_used']['walltime']),
        "mem_request": int(job_info['Resource_List']['mem'].replace("gb", "")),
        "mem_used": format_mem(job_info['resources_used']['mem']),
        "walltime_request": job_info['Resource_List']['walltime'],
        "walltime_used": job_info['resources_used']['walltime']
    }
def get_jobs_info(job_ids: list  = None) -> list[dict]:
    jobs_data = []
    qstat_command = QSTAT_COMMAND
    qstat_result = subprocess.run(f"{qstat_command} {' '.join(job_ids)}", shell=True, check=True, text=True, encoding="utf-8", capture_output=True)
    qstat_json = json.loads(qstat_result.stdout)
    for job_id in job_ids:
        job_id = f"{job_id}.{PBS_SERVER}"
        job_data = parse_job_info(qstat_json["Jobs"][job_id])
        job_data["job_id"] = job_id
        jobs_data.append(job_data)
    return jobs_data
from pathlib import Path
import json
import subprocess
import re 
from datetime import timedelta
from typing import Union

#get "whatshap_phasing" from "snakejob.whatshap_phasing.999.sh"
job_name_regex = re.compile(r"snakejob\.(\w+)\.\d+\.sh")
QSTAT_COMMAND = "qstat -fx -F json"
PBS_SERVER = "hpcpbs02"


def format_mem(mem: str) -> float:
    mem = int(mem.replace("kb", ""))
    return round(mem / 1024 / 1024, 2)

def parse_hms(hms: str) -> timedelta:
    h, m, s = map(int, hms.split(":"))
    return timedelta(hours=h, minutes=m, seconds=s)

def cal_average_cpu(cpu_time: str, walltime_used: str) -> float:
    cpu_time = parse_hms(cpu_time)
    walltime_used = parse_hms(walltime_used)
    return round(cpu_time / walltime_used, 2)

def parse_job_info(job_info: dict) -> dict:
    job_name = job_name_regex.search(job_info["Job_Name"]).group(1)
    exec_host = job_info['exec_host'].split("/")[0]
    return {
        "job_name": job_name,
        "exec_host": exec_host,
        "cpu_request": int(job_info['Resource_List']['ncpus']), 
        "cpu_max": job_info['resources_used']['cpupercent'] / 100,
        "cpu_average": cal_average_cpu(job_info['resources_used']['cput'], job_info['resources_used']['walltime']),
        "mem_request": int(job_info['Resource_List']['mem'].replace("gb", "")),
        "mem_used": format_mem(job_info['resources_used']['mem']),
        "walltime_request": job_info['Resource_List']['walltime'],
        "walltime_used": job_info['resources_used']['walltime']
    }

def get_jobs_info(job_ids: list  = None) -> list[dict]:
    jobs_data = []
    qstat_command = QSTAT_COMMAND
    qstat_result = subprocess.run(f"{qstat_command} {' '.join(job_ids)}", shell=True, check=True, text=True, encoding="utf-8", capture_output=True)
    qstat_json = json.loads(qstat_result.stdout)
    for job_id in job_ids:
        job_id = f"{job_id}.{PBS_SERVER}"
        job_data = parse_job_info(qstat_json["Jobs"][job_id])
        job_data["job_id"] = job_id
        jobs_data.append(job_data)
    return jobs_data

def print_table_simple(jobs_data: list[dict]):
    """Print job data as a simple formatted table using string formatting"""
    if not jobs_data:
        print("No job data to display")
        return
    
    # Define column headers and their widths
    headers = {
        "job_id": ("Job ID", 25),
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
    
    # Print header
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

def print_table_tabulate(jobs_data: list[dict]):
    """Print job data using tabulate library (needs to be installed: pip install tabulate)"""
    try:
        from tabulate import tabulate
        
        if not jobs_data:
            print("No job data to display")
            return
        
        # Define headers in display order
        headers = [
            "job_id", "job_name", "exec_host", 
            "cpu_request", "cpu_max", "cpu_average",
            "mem_request", "mem_used", 
            "walltime_request", "walltime_used"
        ]
        
        # Create table data
        table_data = []
        for job in jobs_data:
            row = [job.get(h, "") for h in headers]
            table_data.append(row)
        
        # Define pretty headers
        pretty_headers = [
            "Job ID", "Job Name", "Host", 
            "CPU Req", "CPU Max", "CPU Avg",
            "Mem Req(GB)", "Mem Used(GB)", 
            "Time Req", "Time Used"
        ]
        
        print(tabulate(table_data, headers=pretty_headers, tablefmt="grid"))
        
    except ImportError:
        print("tabulate library not found. Using simple table format instead.")
        print_table_simple(jobs_data)

def print_summary_stats(jobs_data: list[dict]):
    """Print summary statistics of the jobs"""
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
    
    # Group by job name
    job_groups = {}
    for job in jobs_data:
        name = job['job_name']
        if name not in job_groups:
            job_groups[name] = []
        job_groups[name].append(job)
    
    print(f"\nJobs by Type:")
    for name, jobs in sorted(job_groups.items()):
        print(f"  {name}: {len(jobs)} jobs")

if __name__ == "__main__":
    from pprint import pprint
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python qstat_jobs.py <rundir>")
        sys.exit(1)
    
    rundir = Path(sys.argv[1])
    job_ids = [fn.suffix[2:] for fn in rundir.glob("snakejob.*.sh.e*")]
    
    if not job_ids:
        print(f"No job files found in {rundir}")
        sys.exit(1)
    
    print(f"Found {len(job_ids)} jobs in {rundir}")
    jobs_data = get_jobs_info(job_ids)
    
    
    print_table_tabulate(jobs_data)
    
    print_summary_stats(jobs_data)