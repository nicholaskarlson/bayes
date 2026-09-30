#!/usr/bin/env python3
"""Deterministic author-only construction of fictional Unit 5 rows.

Parameter values in this file specify a synthetic teaching artifact. They are
not reader-facing estimates, a revealed causal structure, or real participants.
"""
import csv
import argparse
from pathlib import Path
from random import Random
ROOT=Path(__file__).resolve().parents[1]
TARGET=ROOT/'data/unit05_belonging_study.csv'
def make():
    rng=Random(20260905)
    rows=[]
    for i in range(1,73):
        long=int((17*(i-1))%72>=36)
        hours=round(max(1.0,min(10.0,5.7-1.4*long+((11*i)%23-11)*0.235+rng.uniform(-.23,.23))),1)
        belonging=round(max(0.0,min(100.0,55.0+1.45*(hours-5.0)-6.5*long+rng.gauss(0,7.0))),1)
        rows.append((f'B{i:03d}',f'{hours:.1f}',str(long),f'{belonging:.1f}'))
    return rows
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=TARGET)
    target=parser.parse_args().output
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('w',newline='') as out:
        writer=csv.writer(out,lineterminator='\n')
        writer.writerow(['participant_id','study_hours_per_week','long_commute','belonging_score'])
        writer.writerows(make())
    print(f'UNIT05_GENERATED n=72 path={target}')
