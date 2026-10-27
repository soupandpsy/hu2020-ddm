import time,json,csv
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]/'Output/hddm/matching/diagnostic_H2H3_20261026_1425'
CG=Path('/sys/fs/cgroup')
def cpu():return int(dict(line.split() for line in (CG/'cpu.stat').read_text().splitlines())['usage_usec'])
start=time.time();prev=start;used=cpu();peak=0.;mempeak=0
with (BASE/'resource_samples.csv').open('a',buffering=1) as f:
    w=csv.writer(f)
    if f.tell()==0:w.writerow(['epoch','cpu_cores','cgroup_memory_bytes'])
    for _ in range(720):
        time.sleep(10);now=time.time();newcpu=cpu();cores=(newcpu-used)/(now-prev)/1e6
        memory=int((CG/'memory.current').read_text());w.writerow([now,cores,memory]);prev=now;used=newcpu
        peak=max(peak,cores);mempeak=max(mempeak,memory)
        obj=dict(measurement_started_epoch=start,last_epoch=now,peak_cpu_cores=peak,peak_cgroup_memory_bytes=mempeak,
            latest_cpu_cores=cores,latest_memory_bytes=memory,memory_events=dict(line.split() for line in (CG/'memory.events').read_text().splitlines()),
            note='Container cgroup memory includes cache; process RSS would double-count shared pages. CPU quota=8, memory cap=4 GiB.')
        (BASE/'resource_summary.json').write_text(json.dumps(obj,indent=2))
