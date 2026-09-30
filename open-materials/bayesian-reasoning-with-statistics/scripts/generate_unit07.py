#!/usr/bin/env python3
"""Author-only deterministic fictional volunteers; not a sample from a campus."""
import argparse
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def rows():
    # Explicit outcome positions avoid a false suggestion of random assignment.
    desk_yes={2,4,7,9,13,16,20,23}
    for route,n,yes,prefix in [('desk',24,desk_yes,'D'),('poster',12,set(),'P')]:
        for i in range(1,n+1):
            yield (f'{prefix}{i:03d}',route,int(i in yes))
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'data/unit07_help_seeking_routes.csv')
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(('volunteer_id','recruitment_route','reported_contact_7d'));w.writerows(rows())
    print('UNIT07_DATA_OK n=36 desk=8/24 poster=0/12 '+str(a.output))
if __name__=='__main__':main()
